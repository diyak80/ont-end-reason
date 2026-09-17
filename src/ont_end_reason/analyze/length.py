"""Length distribution analysis per end_reason category (TOOL_SPEC type 4).

For each end_reason class, computes a structured summary of the read-length
distribution: n, mean, median, percentiles (25/50/75/95/99), N50, std, min,
max. Reads from sequencing_summary.txt (streaming) or any iterable of
ReadRecord with a `.length` field.

Implementation is purely numpy/pandas — no external paper-script
dependencies. The result dataclass `LengthResult` is JSON-serialisable for
downstream report generators.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ..errors import AnalysisError
from ..errors import IOError as OntIOError
from ..io.manifest import ReadRecord
from ..io.readers import detect_format, extract_from_summary


@dataclass
class LengthStats:
    """Length distribution summary for one end_reason class."""

    n: int
    mean: float
    median: float
    std: float
    min: int
    max: int
    p25: float
    p50: float
    p75: float
    p95: float
    p99: float
    n50: int  # length at which 50% of cumulative-length lies above


@dataclass
class LengthResult:
    """Per-end_reason length distributions + raw values for downstream viz."""

    total_reads: int = 0
    per_class: dict[str, LengthStats] = field(default_factory=dict)
    # raw_lengths_by_class is kept for visualisation (histograms, violins).
    # Full per-class population, uncapped and unsampled: plot_length_distribution
    # (viz/static.py) renders it at exact 1 bp resolution via exact_counts.py, and
    # capping or subsampling here would make the plotted population diverge from
    # the per_class stats and the JSON output, which are always exact.
    raw_lengths_by_class: dict[str, list[int]] = field(default_factory=dict)
    source: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "total_reads": self.total_reads,
            "per_class": {
                k: {
                    "n": s.n,
                    "mean": round(s.mean, 2),
                    "median": round(s.median, 2),
                    "std": round(s.std, 2),
                    "min": s.min,
                    "max": s.max,
                    "p25": round(s.p25, 2),
                    "p50": round(s.p50, 2),
                    "p75": round(s.p75, 2),
                    "p95": round(s.p95, 2),
                    "p99": round(s.p99, 2),
                    "n50": s.n50,
                }
                for k, s in self.per_class.items()
            },
            "source": self.source,
        }


def _n50(lengths: np.ndarray) -> int:
    """Standard sequencing N50: length at which cumulative sorted-descending
    length crosses half of total length."""
    if len(lengths) == 0:
        return 0
    sorted_desc = np.sort(lengths)[::-1]
    cumulative = np.cumsum(sorted_desc)
    half = cumulative[-1] / 2
    idx = np.searchsorted(cumulative, half)
    return int(sorted_desc[min(idx, len(sorted_desc) - 1)])


def _summarize(lengths: np.ndarray) -> LengthStats:
    if len(lengths) == 0:
        # Defensive — caller filters empty groups
        raise AnalysisError("Cannot summarize empty length array")
    return LengthStats(
        n=len(lengths),
        mean=float(np.mean(lengths)),
        median=float(np.median(lengths)),
        std=float(np.std(lengths)),
        min=int(np.min(lengths)),
        max=int(np.max(lengths)),
        p25=float(np.percentile(lengths, 25)),
        p50=float(np.percentile(lengths, 50)),
        p75=float(np.percentile(lengths, 75)),
        p95=float(np.percentile(lengths, 95)),
        p99=float(np.percentile(lengths, 99)),
        n50=_n50(lengths),
    )


def _gather_lengths_by_class(
    records: Iterable[ReadRecord],
) -> dict[str, list[int]]:
    by_class: dict[str, list[int]] = {}
    for r in records:
        if r.length is None or r.length <= 0:
            continue
        by_class.setdefault(r.end_reason, []).append(r.length)
    return by_class


def _from_summary_streaming(path: Path) -> tuple[dict[str, np.ndarray], int]:
    """Read sequencing_summary.txt in chunks and concat lengths per end_reason."""
    by_class: dict[str, list[np.ndarray]] = {}
    n_total = 0
    try:
        for chunk in pd.read_csv(
            path,
            sep="\t",
            usecols=["end_reason", "sequence_length_template"],
            chunksize=200_000,
        ):
            chunk = chunk[chunk["sequence_length_template"] > 0]
            n_total += len(chunk)
            for er, grp in chunk.groupby("end_reason"):
                by_class.setdefault(str(er), []).append(
                    grp["sequence_length_template"].to_numpy(dtype=np.int64)
                )
    except (OSError, ValueError) as exc:
        raise OntIOError(f"Failed to stream {path}: {exc}") from exc

    return ({k: np.concatenate(v) for k, v in by_class.items()}, n_total)


def length(
    source: str | Path | Iterable[ReadRecord],
) -> LengthResult:
    """Per-end_reason length-distribution summary.

    Accepts:
      - A file path (sequencing_summary.txt) — streaming, memory-bounded
      - An iterable of `ReadRecord` — direct mode

    POD5/Fast5 inputs go through ReadRecord iteration because end_reason
    is read at extraction time. For PromethION-scale data, prefer the
    sequencing_summary.txt path.

    ``raw_lengths_by_class`` on the result is the full per-class population,
    never capped or subsampled: it used to be randomly subsampled to 50,000
    reads per class for the ``--plot`` path only, while this function's own
    per-class stats and JSON output were always computed from the full
    population, so the figure and the numbers describe different populations
    (see lib/readdist/CONTRACT.md in ont-ecosystem, which this fix follows).
    """
    if isinstance(source, (str, Path)):
        path = Path(source)
        fmt = detect_format(path) if path.is_file() else "summary"
        if path.is_file() and fmt != "summary":
            # POD5 / Fast5 — go through ReadRecord pipeline
            records = extract_from_summary(path) if fmt == "summary" else None
            if records is None:
                raise AnalysisError(
                    "POD5/Fast5 length analysis: extract end_reason + length "
                    "via extract_from_summary on the run's sequencing_summary.txt"
                )
            return length(records)
        # Stream the summary
        arrays_by_class, n_total = _from_summary_streaming(path)
    else:
        by_class_lists = _gather_lengths_by_class(source)
        arrays_by_class = {k: np.asarray(v, dtype=np.int64) for k, v in by_class_lists.items()}
        n_total = sum(len(v) for v in arrays_by_class.values())

    if n_total == 0:
        raise AnalysisError("No reads with positive length found")

    per_class: dict[str, LengthStats] = {}
    raw: dict[str, list[int]] = {}
    for er, lengths in arrays_by_class.items():
        if len(lengths) == 0:
            continue
        per_class[er] = _summarize(lengths)
        # Full population, no cap and no subsampling: the plot must show
        # exactly what per_class/JSON report, not a random draw from it.
        raw[er] = lengths.tolist()

    return LengthResult(
        total_reads=n_total,
        per_class=per_class,
        raw_lengths_by_class=raw,
        source=str(source) if isinstance(source, (str, Path)) else None,
    )
