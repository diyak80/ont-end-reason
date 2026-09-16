"""Table generation for the paper's tables + supplementary tables (TOOL_SPEC 12-13).

Composes results from `analyze.*` into a single table object, then renders
to TSV / CSV / Markdown. The table object is the primary in-memory data
structure; rendering is downstream.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd

from ..errors import AnalysisError
from .distribution import distribution
from .length import length
from .quality import quality


@dataclass
class TableResult:
    name: str = ""
    rows: list[dict[str, Any]] = field(default_factory=list)
    format: str = "tsv"  # tsv | csv | markdown | latex
    path: str | None = None

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)

    def render(self, fmt: str | None = None) -> str:
        fmt = (fmt or self.format).lower()
        df = self.to_dataframe()
        if fmt == "tsv":
            return df.to_csv(sep="\t", index=False)
        if fmt == "csv":
            return df.to_csv(index=False)
        if fmt in ("markdown", "md"):
            return df.to_markdown(index=False) or ""
        if fmt == "latex":
            return df.to_latex(index=False)
        raise AnalysisError(f"Unknown table format: {fmt!r}")


_KNOWN_NAMES = {"summary", "per_class", "quality"}


def generate_tables(source: str | Path, *, name: str, **kwargs: Any) -> TableResult:
    """Generate a named table from a sequencing_summary.txt or POD5 directory.

    Supported names:
      summary    one-row summary of total reads + status + key percentages
      per_class  one row per end_reason with n / pct / median-length / median-q
      quality    one row per end_reason with Q-score summary + GMM k
    """
    fmt = kwargs.pop("fmt", "tsv")
    if name not in _KNOWN_NAMES:
        raise AnalysisError(f"Unknown table name {name!r}. Choose from: {sorted(_KNOWN_NAMES)}")

    if name == "summary":
        dist = distribution(source)
        payload = dist.to_dict()
        rows = [
            {
                "total_reads": dist.total_reads,
                "quality_status": dist.quality_status,
                "signal_positive_pct": payload["signal_positive_pct"],
                "unblock_mux_pct": payload["unblock_mux_pct"],
                "data_service_pct": payload["data_service_pct"],
            }
        ]
    elif name == "per_class":
        dist = distribution(source)
        # Length only meaningful for summary/POD5; skip if missing
        try:
            len_res = length(source)
            len_per_class = dict(len_res.per_class.items())
        except AnalysisError:
            len_per_class = {}
        rows = []
        for er, count in dist.counts.items():
            row: dict[str, Any] = {
                "end_reason": er,
                "n": count,
                "pct": round(dist.percentages.get(er, 0) * 100, 2),
            }
            stats = len_per_class.get(er)
            if stats is not None:
                row.update(
                    {
                        "median_length": round(stats.median, 1),
                        "p95_length": round(stats.p95, 1),
                        "n50": stats.n50,
                    }
                )
            rows.append(row)
    elif name == "quality":
        qual = quality(source)
        rows = []
        for er, s in qual.per_class.items():
            rows.append(
                {
                    "end_reason": er,
                    "n": s.n,
                    "mean_q": round(s.mean, 3),
                    "median_q": round(s.median, 3),
                    "p25_q": round(s.p25, 3),
                    "p75_q": round(s.p75, 3),
                    "gmm_chosen_k": s.gmm_chosen_k,
                }
            )
    else:  # pragma: no cover — guarded above
        raise AnalysisError(name)

    return TableResult(name=name, rows=rows, format=fmt)


def render_table(table: TableResult, output_path: str | Path | None = None) -> str:
    """Render a TableResult to its configured format. If `output_path` is given,
    also writes the rendered string to disk."""
    text = table.render()
    if output_path is not None:
        p = Path(output_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        table.path = str(p)
    return text


# Convenience for the report layer: render every supported table as a dict
def render_all(source: str | Path) -> dict[str, str]:
    """Render summary + per_class + quality tables as markdown strings."""
    out: dict[str, str] = {}
    for name in ("summary", "per_class", "quality"):
        try:
            t = generate_tables(source, name=name, fmt="markdown")
            out[name] = t.render()
        except AnalysisError as exc:
            out[name] = f"_(table not available: {exc})_"
    return out


# Make StringIO not flagged as unused
_ = StringIO
