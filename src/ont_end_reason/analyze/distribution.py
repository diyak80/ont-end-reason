"""End-reason distribution analysis — the v0.1.0 fully-implemented analysis.

Computes per-end_reason counts and percentages, then applies a quality gate
based on the canonical paper thresholds:

  signal_positive >= 75%  → OK
  signal_positive <  75%  → CHECK
  signal_positive <  50%  → FAIL

This is the workhorse analysis equivalent to `/end-reason analyze` in the lab
ecosystem. It accepts any iterable of `ReadRecord` (typically from
`io.extract_from_*`) or any Manifest produced by `io.discover`.
"""

from __future__ import annotations

import contextlib
import logging
import sys
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..codes import CODES, NAMES
from ..errors import AnalysisError
from ..io.manifest import Manifest, ReadRecord
from ..io.readers import (
    detect_format,
    extract_from_fast5,
    extract_from_pod5,
    extract_from_summary,
)

logger = logging.getLogger(__name__)

DEFAULT_REGISTRY_PATH = Path.home() / ".ont-registry" / "experiments.yaml"
ONT_ECOSYSTEM_PATH = Path.home() / "repos" / "ont-ecosystem"


@dataclass
class DistributionResult:
    """Structured result of an end-reason distribution analysis.

    `counts` is keyed by the canonical lower-case full name
    (`signal_positive`, `unblock_mux_change`, ...). `percentages` mirrors
    counts as ratios of `total_reads`.
    """

    total_reads: int
    counts: dict[str, int] = field(default_factory=dict)
    percentages: dict[str, float] = field(default_factory=dict)
    quality_status: str = "UNKNOWN"  # OK | CHECK | FAIL | UNKNOWN
    signal_positive_pct: float | None = None
    unblock_mux_pct: float | None = None
    data_service_pct: float | None = None
    interpretation: str = ""
    source_format: str | None = None

    @property
    def recognized_reads(self) -> int:
        return sum(
            n for reason, n in self.counts.items() if reason in CODES and reason != "unknown"
        )

    @property
    def unrecognized_reads(self) -> int:
        return self.total_reads - self.recognized_reads

    def to_dict(self) -> dict[str, object]:
        return {
            "total_reads": self.total_reads,
            "quality_status": self.quality_status,
            "recognized_reads": self.recognized_reads,
            "unrecognized_reads": self.unrecognized_reads,
            "signal_positive_pct": round(self.signal_positive_pct, 2)
            if self.signal_positive_pct is not None
            else None,
            "unblock_mux_pct": round(self.unblock_mux_pct, 2)
            if self.unblock_mux_pct is not None
            else None,
            "data_service_pct": round(self.data_service_pct, 2)
            if self.data_service_pct is not None
            else None,
            "counts": self.counts,
            "percentages": {k: round(v, 4) for k, v in self.percentages.items()},
            "interpretation": self.interpretation,
            "source_format": self.source_format,
        }


def _aggregate(records: Iterable[ReadRecord]) -> tuple[int, dict[str, int]]:
    counts: dict[str, int] = {}
    total = 0
    for r in records:
        total += 1
        key = (r.end_reason or "").strip() or "unknown"
        counts[key] = counts.get(key, 0) + 1
    return total, counts


def _interpretation(percentages: Mapping[str, float]) -> str:
    sp_pct = percentages.get("signal_positive", 0.0) * 100
    umc_pct = percentages.get("unblock_mux_change", 0.0) * 100
    parts = []
    if sp_pct >= 95:
        parts.append(f"Excellent: {sp_pct:.1f}% complete reads (top of expected range).")
    elif sp_pct >= 80:
        parts.append(f"Healthy: {sp_pct:.1f}% complete reads.")
    elif sp_pct >= 75:
        parts.append(f"Acceptable: {sp_pct:.1f}% complete reads (lower bound of OK).")
    elif sp_pct >= 50:
        parts.append(f"Concerning: {sp_pct:.1f}% complete reads (below 75% threshold).")
    else:
        parts.append(f"Failed: {sp_pct:.1f}% complete reads (below 50% threshold).")
    if umc_pct >= 10:
        parts.append(f"Adaptive sampling rejection rate is {umc_pct:.1f}%.")
    return " ".join(parts)


def distribution(
    source: str | Path | Manifest | Iterable[ReadRecord],
    *,
    quick: bool = False,
    max_reads: int = 10_000,
) -> DistributionResult:
    """Compute end-reason distribution for any acceptable input.

    Accepts:

    - A file path (POD5/Fast5/sequencing_summary.txt) — auto-detects format
    - A directory path — uses the first viable format found
    - A `Manifest` — concatenates all summaries in it
    - An iterable of `ReadRecord` — direct mode for in-memory analysis

    With `quick=True`, only the first `max_reads` are read. Useful for
    fast preview on multi-million-read POD5 sets.
    """
    src_format: str | None = None
    records: Iterable[ReadRecord]

    if isinstance(source, Manifest):
        if source.summaries:
            src_format = "summary"
            records = list(
                extract_from_summary(source.summaries[0].path, quick=quick, max_reads=max_reads)
            )
        elif source.pod5:
            src_format = "pod5"
            records = extract_from_pod5(source.pod5[0].path, quick=quick, max_reads=max_reads)
        elif source.fast5:
            src_format = "fast5"
            records = extract_from_fast5(source.fast5[0].path, quick=quick, max_reads=max_reads)
        else:
            raise AnalysisError("Manifest has no POD5/Fast5/summary files")
    elif isinstance(source, (str, Path)):
        src_format = detect_format(source)
        if src_format == "pod5":
            records = extract_from_pod5(source, quick=quick, max_reads=max_reads)
        elif src_format == "fast5":
            records = extract_from_fast5(source, quick=quick, max_reads=max_reads)
        elif src_format == "summary":
            records = list(extract_from_summary(source, quick=quick, max_reads=max_reads))
        else:  # pragma: no cover — detect_format raises otherwise
            raise AnalysisError(f"Unsupported format: {src_format}")
    else:
        records = source  # raw iterable

    total, counts = _aggregate(records)
    if total == 0:
        raise AnalysisError("No reads found in source")

    # Normalise keys to canonical lower-case full names. Unknown values pass through.
    norm_counts: dict[str, int] = {}
    for key, count in counts.items():
        canonical = key.lower() if key.lower() in CODES else key
        # Allow short codes too: 'SP' → 'signal_positive'
        if canonical.upper() in NAMES:
            canonical = NAMES[canonical.upper()]
        norm_counts[canonical] = norm_counts.get(canonical, 0) + count

    percentages = {k: v / total for k, v in norm_counts.items()}

    sp_pct: float | None = percentages.get("signal_positive", 0.0) * 100
    umc_pct: float | None = percentages.get("unblock_mux_change", 0.0) * 100
    dumc_pct: float | None = percentages.get("data_service_unblock_mux_change", 0.0) * 100

    unrecognized = sum(
        n for reason, n in norm_counts.items() if reason not in CODES or reason == "unknown"
    )
    interpretation = _interpretation(percentages)
    if unrecognized:
        status = "UNKNOWN"
        sp_pct = umc_pct = dumc_pct = None
        interpretation = (
            f"End-reason quality is unavailable: {total - unrecognized:,}/{total:,} reads "
            "have recognized end-reason metadata. Counts and fractions describe recorded labels; "
            "they do not establish run quality."
        )
    elif sp_pct < 50:
        status = "FAIL"
    elif sp_pct < 75:
        status = "CHECK"
    else:
        status = "OK"

    return DistributionResult(
        total_reads=total,
        counts=norm_counts,
        percentages=percentages,
        quality_status=status,
        signal_positive_pct=sp_pct,
        unblock_mux_pct=umc_pct,
        data_service_pct=dumc_pct,
        interpretation=interpretation,
        source_format=src_format,
    )


# ──────────────────── qc_baseline auto-population (Phase 5) ────────────────────


def _resolve_path(p: str | Path) -> Path | None:
    """Return an absolute, symlink-resolved Path or None if it can't be resolved."""
    try:
        return Path(p).expanduser().resolve()
    except (OSError, RuntimeError):
        return None


def _load_registry_entries(
    registry_path: Path | None = None,
) -> list[dict[str, Any]]:
    """Load the ONT registry experiments list. Returns [] if unavailable."""
    if registry_path is None:
        registry_path = DEFAULT_REGISTRY_PATH
    if not registry_path.exists():
        return []
    try:
        import yaml  # pyyaml is a hard dependency of ont-end-reason
    except ImportError:  # pragma: no cover — pyyaml is in pyproject deps
        return []
    try:
        data = yaml.safe_load(registry_path.read_text()) or {}
    except (yaml.YAMLError, OSError) as exc:
        logger.warning("Could not parse ONT registry at %s: %s", registry_path, exc)
        return []
    entries = data.get("experiments") or []
    if not isinstance(entries, list):
        return []
    return entries


def _match_registry_entry(
    source_path: str | Path,
    entries: Iterable[dict[str, Any]],
) -> dict[str, Any] | None:
    """Return the registry entry whose `location` matches source_path, else None.

    Match is on absolute-resolved paths. A registry entry with a parent path
    that contains the source also counts (the registry stores experiment
    directories, but `analyze distribution` may be invoked on a file within).
    """
    src = _resolve_path(source_path)
    if src is None:
        return None
    for entry in entries:
        loc = entry.get("location")
        if not loc:
            continue
        loc_path = _resolve_path(loc)
        if loc_path is None:
            continue
        if src == loc_path:
            return entry
        # Also accept: source is inside the registered experiment directory
        try:
            src.relative_to(loc_path)
            return entry
        except ValueError:
            continue
    return None


def _parse_timestamp(value: Any) -> float:
    """Parse an ISO-8601 string into a POSIX timestamp; fall back to now()."""
    if isinstance(value, str):
        try:
            dt = datetime.fromisoformat(value)
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
            return dt.timestamp()
        except ValueError:
            pass
    elif isinstance(value, datetime):
        dt = value if value.tzinfo else value.replace(tzinfo=timezone.utc)
        return dt.timestamp()
    return datetime.now(tz=timezone.utc).timestamp()


def maybe_store_baseline(
    result: DistributionResult,
    source_path: str | Path,
    *,
    write: bool = True,
    registry_path: Path | None = None,
) -> str | None:
    """Auto-store a DistributionResult in the qc_baseline store.

    Looks up `source_path` in the ONT registry (`~/.ont-registry/experiments.yaml`).
    If a registry entry's `location` matches (or contains) source_path, builds
    an ExperimentMetadata + QCResult from the entry and stores it via
    `lib.qc_baseline.store_qc_result`.

    Idempotent: result_id is hashed from experiment_id + timestamp; re-running
    with the same registry entry's `discovered` timestamp produces the same id
    and overwrites the same file.

    Behavior:
      * `write=False` → no-op, returns None (the opt-out path).
      * `source_path` doesn't match any registry entry → returns None silently.
      * `lib.qc_baseline` can't be imported → logs a warning, returns None.
      * Successful store → returns the result_id.

    Never raises; on any unexpected error logs a warning and returns None so
    that the distribution analysis itself stays uninterrupted.
    """
    if not write:
        return None
    if (
        result.quality_status == "UNKNOWN"
        or result.unrecognized_reads
        or any(
            value is None
            for value in (
                result.signal_positive_pct,
                result.unblock_mux_pct,
                result.data_service_pct,
            )
        )
    ):
        logger.warning("End-reason quality unavailable; skipping baseline store")
        return None

    entries = _load_registry_entries(registry_path)
    if not entries:
        return None

    entry = _match_registry_entry(source_path, entries)
    if entry is None:
        return None

    # Graceful import of the cross-repo qc_baseline module.
    ont_eco_str = str(ONT_ECOSYSTEM_PATH)
    inserted = False
    if ont_eco_str not in sys.path and ONT_ECOSYSTEM_PATH.exists():
        sys.path.insert(0, ont_eco_str)
        inserted = True
    try:
        from lib.qc_baseline import (  # type: ignore[import-not-found]
            ExperimentMetadata,
            QCResult,
            store_qc_result,
        )
    except ImportError as exc:
        logger.warning(
            "qc_baseline not available (ont-ecosystem missing at %s): %s; skipping baseline store",
            ONT_ECOSYSTEM_PATH,
            exc,
        )
        if inserted:
            with contextlib.suppress(ValueError):
                sys.path.remove(ont_eco_str)
        return None

    try:
        metadata = ExperimentMetadata(
            experiment_id=str(entry.get("id") or entry.get("name") or "unknown"),
            flowcell_type=entry.get("flowcell_type") or entry.get("flowcell"),
            chemistry=entry.get("chemistry"),
            sample_type=entry.get("sample_type"),
            run_duration_hours=entry.get("run_duration_hours"),
            basecaller_model=entry.get("basecaller_model"),
            adaptive_sampling=bool(entry.get("adaptive_sampling", False)),
            notes=entry.get("notes", "") or "",
        )

        metrics: dict[str, float] = {
            "signal_positive_pct": float(result.signal_positive_pct),
            "unblock_mux_pct": float(result.unblock_mux_pct),
            "data_service_pct": float(result.data_service_pct),
            "mux_change_pct": float(result.percentages.get("mux_change", 0.0) * 100.0),
            "signal_negative_pct": float(result.percentages.get("signal_negative", 0.0) * 100.0),
        }

        ts = _parse_timestamp(entry.get("discovered"))

        qc_result = QCResult(
            metadata=metadata,
            metrics=metrics,
            timestamp=ts,
        )
        return store_qc_result(qc_result)
    except Exception as exc:  # pragma: no cover — defensive guard
        logger.warning("Failed to store qc_baseline result: %s", exc)
        return None
