"""Format-specific readers for POD5, Fast5, and sequencing_summary.txt.

Each reader returns a list of `ReadRecord` dataclasses with normalised
field names so downstream analyses don't need format-specific branches.

These functions are direct ports of the equivalent helpers in
ont-ecosystem/skills/end-reason/scripts/end_reason.py (commit ba2f9a4f),
adapted to use the package's `ReadRecord` dataclass and `OntIOError`
instead of plain dicts and ValueError.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import structlog

from ..codes import CODES
from ..errors import IOError as OntIOError
from .manifest import ReadRecord

log = structlog.get_logger(__name__)


def _normalise_end_reason(reason: Any) -> str:
    """Coerce whatever the upstream library returns to the canonical lower-case
    ONT name (signal_positive, unblock_mux_change, ...).

    Handles all observed shapes:
      "signal_positive"
      "EndReason.signal_positive"
      "<EndReason.signal_positive: 4>"
      "(<EndReason.signal_positive: 4>, forced=False)"  ← pod5 NamedTuple repr
      EndReason enum members (any class with .name)
    """
    import re

    if reason is None or (isinstance(reason, float) and reason != reason):
        return "unknown"
    # Enum members have a .name attribute; use it directly when present
    name_attr = getattr(reason, "name", None)
    if name_attr and isinstance(name_attr, str):
        return name_attr.lower()
    s = str(reason).strip()
    # Strip surrounding parens/angle-brackets and pod5's NamedTuple cruft
    # Examples we want to reduce to the bare name:
    #   "(<EndReason.signal_positive: 4>, forced=False)"
    #   "<EndReason.signal_positive: 4>"
    #   "EndReason.signal_positive"
    match = re.search(r"EndReason\.(\w+)", s)
    if match:
        return match.group(1).lower()
    # Fallback: take whatever follows the last dot, strip non-word chars
    if "." in s:
        s = s.rsplit(".", 1)[-1]
    s = re.sub(r"[^\w]", "", s).lower()
    return s or "unknown"


def detect_format(path: str | Path) -> str:
    """Return one of {"pod5", "fast5", "summary"} based on path inspection."""
    p = Path(path)
    if p.is_file():
        suffix = p.suffix.lower()
        if suffix == ".pod5":
            return "pod5"
        if suffix == ".fast5":
            return "fast5"
        if p.name.startswith("sequencing_summary") and p.name.endswith(".txt"):
            return "summary"
        raise OntIOError(f"Cannot detect format of file: {p}")
    if p.is_dir():
        # Directory: prefer POD5 > Fast5 > summary
        if any(p.glob("*.pod5")):
            return "pod5"
        if any(p.glob("*.fast5")):
            return "fast5"
        if any(p.glob("sequencing_summary*.txt")):
            return "summary"
        if any(p.rglob("*.pod5")):
            return "pod5"
        if any(p.rglob("*.fast5")):
            return "fast5"
        raise OntIOError(f"No POD5/Fast5/summary files found under {p}")
    raise OntIOError(f"Path not found: {p}")


def extract_from_pod5(
    path: str | Path,
    *,
    quick: bool = False,
    max_reads: int = 10_000,
) -> list[ReadRecord]:
    """Pull end_reason metadata from a POD5 file or directory of POD5 files.

    If `path` is a directory, every `*.pod5` under it is read. With
    `quick=True`, stops once `max_reads` are extracted across all files.
    """
    try:
        import pod5  # noqa: F401 - imported lazily inside this function
    except ImportError as exc:
        raise OntIOError("POD5 reader requires the `pod5` package: pip install pod5") from exc

    p = Path(path)
    files = sorted(p.rglob("*.pod5")) if p.is_dir() else [p]
    if not files:
        raise OntIOError(f"No POD5 files found under {p}")

    out: list[ReadRecord] = []
    for f in files:
        try:
            from pod5 import Reader as Pod5Reader  # type: ignore[import-not-found]
        except ImportError as exc:
            raise OntIOError("pod5 library missing") from exc

        try:
            with Pod5Reader(str(f)) as reader:
                for read in reader.reads():
                    er = _normalise_end_reason(getattr(read, "end_reason", "unknown"))
                    out.append(
                        ReadRecord(
                            read_id=str(read.read_id),
                            end_reason=er,
                            end_reason_short=CODES.get(er),
                            source_file=str(f),
                            source_format="pod5",
                        )
                    )
                    if quick and len(out) >= max_reads:
                        return out
        except Exception as exc:
            log.warning("pod5 read failed", path=str(f), error=str(exc))
            raise OntIOError(f"Failed to read POD5 {f}: {exc}") from exc

    return out


def extract_from_fast5(
    path: str | Path,
    *,
    quick: bool = False,
    max_reads: int = 10_000,
) -> list[ReadRecord]:
    """Pull end_reason metadata from a Fast5 file or directory of Fast5 files.

    Fast5 is the legacy format; many older runs only have Fast5. This reader
    handles both single- and multi-read Fast5 files.
    """
    try:
        import h5py
    except ImportError as exc:
        raise OntIOError("Fast5 reader requires h5py: pip install h5py") from exc

    import h5py  # type: ignore[no-redef]

    p = Path(path)
    files = sorted(p.rglob("*.fast5")) if p.is_dir() else [p]
    if not files:
        raise OntIOError(f"No Fast5 files found under {p}")

    out: list[ReadRecord] = []
    for f in files:
        try:
            with h5py.File(str(f), "r") as h5:
                for read_group_name in h5:
                    grp = h5[read_group_name]
                    # End reason lives in the Raw/Reads attrs in older versions
                    raw = grp.get("Raw")
                    if raw is None:
                        continue
                    end_reason_value = raw.attrs.get("end_reason", b"unknown")
                    if isinstance(end_reason_value, bytes):
                        end_reason_value = end_reason_value.decode()
                    er = _normalise_end_reason(end_reason_value)
                    rid = grp.attrs.get("read_id", read_group_name)
                    if isinstance(rid, bytes):
                        rid = rid.decode()
                    out.append(
                        ReadRecord(
                            read_id=str(rid),
                            end_reason=er,
                            end_reason_short=CODES.get(er),
                            source_file=str(f),
                            source_format="fast5",
                        )
                    )
                    if quick and len(out) >= max_reads:
                        return out
        except OSError as exc:
            log.warning("fast5 read failed", path=str(f), error=str(exc))
            raise OntIOError(f"Failed to read Fast5 {f}: {exc}") from exc

    return out


def extract_from_summary(
    path: str | Path,
    *,
    quick: bool = False,
    max_reads: int = 10_000,
    chunk_size: int = 100_000,
) -> Iterable[ReadRecord]:
    """Stream `read_id` and `end_reason` columns from sequencing_summary.txt.

    Yields one `ReadRecord` per line. Streaming via pandas chunked TSV reader
    keeps memory bounded even on PromethION-scale (multi-GB) files.

    Yields an iterator rather than a list so callers can stream-aggregate
    without materialising 100M reads in memory.
    """
    import pandas as pd

    p = Path(path)
    if not p.exists():
        raise OntIOError(f"sequencing_summary not found: {p}")

    n = 0
    try:
        for chunk in pd.read_csv(
            p,
            sep="\t",
            usecols=lambda c: (
                c in {"read_id", "end_reason", "sequence_length_template", "mean_qscore_template"}
            ),
            chunksize=chunk_size,
            low_memory=False,
        ):
            for row in chunk.itertuples(index=False):
                er = _normalise_end_reason(getattr(row, "end_reason", "unknown"))
                yield ReadRecord(
                    read_id=str(row.read_id),
                    end_reason=er,
                    end_reason_short=CODES.get(er),
                    length=int(getattr(row, "sequence_length_template", 0) or 0) or None,
                    quality=float(getattr(row, "mean_qscore_template", 0) or 0) or None,
                    source_file=str(p),
                    source_format="summary",
                )
                n += 1
                if quick and n >= max_reads:
                    return
    except (OSError, ValueError) as exc:
        raise OntIOError(f"Failed to read summary {p}: {exc}") from exc
