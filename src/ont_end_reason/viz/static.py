"""matplotlib-backed plotters returning publication-ready `Figure` objects.

Each function accepts the corresponding analysis result dataclass and
returns a `matplotlib.figure.Figure`. The caller saves / shows / composes
as needed — these functions deliberately do NOT call `plt.show()` or
`fig.savefig()`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ..analyze.atlas import END_REASON_METRIC_KEYS, AtlasResult
from ..analyze.distribution import DistributionResult
from ..analyze.length import LengthResult
from ..analyze.temporal import TemporalResult
from ..codes import CODES

if TYPE_CHECKING:
    from matplotlib.figure import Figure


# Short codes + colors for atlas summary's end-reason metric set
_ATLAS_METRIC_LABELS = {
    "signal_positive_pct": "SP",
    "unblock_mux_pct": "UMC",
    "data_service_pct": "DUMC",
    "mux_change_pct": "MC",
    "signal_negative_pct": "SN",
}
_ATLAS_METRIC_COLORS = {
    "signal_positive_pct": "#2ca02c",  # green
    "unblock_mux_pct": "#ff7f0e",  # orange
    "data_service_pct": "#d62728",  # red
    "mux_change_pct": "#1f77b4",  # blue
    "signal_negative_pct": "#7f7f7f",  # gray
}


# Canonical category ordering (most-kept on left → most-rejected on right)
_CATEGORY_ORDER = [
    "signal_positive",
    "unblock_mux_change",
    "data_service_unblock_mux_change",
    "mux_change",
    "partial",
    "unknown",
    "signal_negative",
]

# Canonical colour palette aligned with paper Figure 3.
_PALETTE = {
    "signal_positive": "#2ca02c",
    "unblock_mux_change": "#ff7f0e",
    "data_service_unblock_mux_change": "#d62728",
    "mux_change": "#9467bd",
    "partial": "#8c564b",
    "unknown": "#7f7f7f",
    "signal_negative": "#1f77b4",
}


def plot_distribution(
    result: DistributionResult,
    *,
    title: str | None = None,
    show_pct: bool = True,
) -> Figure:
    """Bar chart of end-reason distribution.

    Parameters
    ----------
    result : DistributionResult
        Output of `analyze.distribution(...)`.
    title : str, optional
        Plot title. Default includes the quality status.
    show_pct : bool
        Annotate bars with percentages. Default True.
    """
    import matplotlib.pyplot as plt  # local import keeps cold-start fast

    ordered_keys = [k for k in _CATEGORY_ORDER if k in result.counts]
    # Append any unexpected categories at the end
    for k in result.counts:
        if k not in ordered_keys:
            ordered_keys.append(k)

    values = [result.counts[k] for k in ordered_keys]
    short_labels = [CODES.get(k, k.upper()[:4]) for k in ordered_keys]
    colours = [_PALETTE.get(k, "#cccccc") for k in ordered_keys]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(short_labels, values, color=colours, edgecolor="black", linewidth=0.5)

    if show_pct:
        for bar, key, _count in zip(bars, ordered_keys, values, strict=True):
            pct = result.percentages.get(key, 0.0) * 100
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{pct:.1f}%",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ax.set_ylabel("Read count")
    ax.set_xlabel("End reason")
    if title is None:
        title = f"End reason distribution — {result.quality_status}  (n = {result.total_reads:,})"
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def plot_length_distribution(
    result: LengthResult,
    *,
    title: str | None = None,
    log_x: bool = True,
) -> Figure:
    """Exact 1 bp read-length curve per end_reason category, no cap, no rebin.

    Every read in ``result.raw_lengths_by_class`` (the full, uncapped
    population) is plotted: :func:`exact_length_counts` builds the exact
    per-base-pair count table and asserts that the plotted counts conserve
    the read total and that no two distinct lengths were merged, the same
    two invariants ont-ecosystem's ``lib.readdist.core`` asserts (see
    ``lib/readdist/CONTRACT.md`` section 3). Previously this histogrammed a
    random 50,000-read subsample per class into 60 log-spaced bins, so the
    figure's largest class disagreed with its own title and JSON output.

    Default x-axis is log-scale because ONT read-length distributions
    typically span 2-3 orders of magnitude; this is an axis-scale choice; it
    does not merge or resample any point. Pass `log_x=False` for linear.
    """
    import matplotlib.pyplot as plt

    from ..exact_counts import exact_length_counts

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ordered = [k for k in _CATEGORY_ORDER if k in result.raw_lengths_by_class]
    for k in result.raw_lengths_by_class:
        if k not in ordered:
            ordered.append(k)

    for key in ordered:
        lengths = result.raw_lengths_by_class[key]
        if not lengths:
            continue
        grid, counts = exact_length_counts(lengths)
        color = _PALETTE.get(key, "#cccccc")
        ax.plot(
            grid,
            counts,
            lw=0.8,
            color=color,
            label=f"{CODES.get(key, key)} (n={len(lengths):,})",
        )
        ax.fill_between(grid, counts, alpha=0.35, color=color, linewidth=0)

    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Read count (exact, 1 bp bins)")
    ax.set_title(title or f"Read length distribution by end_reason (n = {result.total_reads:,})")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def plot_temporal(
    result: TemporalResult,
    *,
    title: str | None = None,
    show_fractions: bool = True,
) -> Figure:
    """Stacked-area plot of end_reason fractions over time.

    If `show_fractions=True` (default), the y-axis is per-bin fraction
    summing to 1.0. If False, raw counts are plotted (useful for spotting
    throughput drops).
    """
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(9, 4.5))
    centers = np.asarray(result.bin_centers, dtype=float)
    ordered = [k for k in _CATEGORY_ORDER if k in result.fractions_by_class]
    for k in result.fractions_by_class:
        if k not in ordered:
            ordered.append(k)

    data_dict = result.fractions_by_class if show_fractions else result.counts_by_class
    stack = np.array([data_dict[k] for k in ordered], dtype=float)
    colours = [_PALETTE.get(k, "#cccccc") for k in ordered]
    labels = [CODES.get(k, k) for k in ordered]
    ax.stackplot(centers, stack, labels=labels, colors=colours, alpha=0.85)

    ax.set_xlabel("Time since run start (hours)")
    ax.set_ylabel("Fraction of reads" if show_fractions else "Read count")
    ax.set_title(title or f"End reason rates over time  (n = {result.total_reads:,})")
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_fractions:
        ax.set_ylim(0, 1.0)
    fig.tight_layout()
    return fig


def plot_quality_violins(result, *, title: str | None = None):
    """Violin plot of Q-score distribution per end_reason category.

    Imports QualityResult lazily to avoid pulling scipy at top-level when
    only distribution/length are used.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from ..analyze.quality import QualityResult

    if not isinstance(result, QualityResult):  # pragma: no cover — type guard
        raise TypeError(f"Expected QualityResult, got {type(result).__name__}")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ordered = [k for k in _CATEGORY_ORDER if k in result.raw_qscores_by_class]
    for k in result.raw_qscores_by_class:
        if k not in ordered:
            ordered.append(k)

    data = [result.raw_qscores_by_class[k] for k in ordered]
    labels = [CODES.get(k, k.upper()[:4]) for k in ordered]
    parts = ax.violinplot(data, showmeans=True, showmedians=False, showextrema=False)
    for body, key in zip(parts["bodies"], ordered, strict=True):
        body.set_facecolor(_PALETTE.get(key, "#cccccc"))
        body.set_alpha(0.7)
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Mean Q-score")
    ax.set_xlabel("End reason")
    ax.set_title(title or f"Q-score distribution by end_reason  (n = {result.total_reads:,})")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def plot_umc_posterior(result, *, title: str | None = None):
    """Two-panel plot of UMC posterior:
    left  — observed length histogram vs posterior expected true length histogram
    right — per-read "bonus" (E[true] − observed) distribution
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from ..analyze.umc_posterior import UMCPosteriorResult

    if not isinstance(result, UMCPosteriorResult):  # pragma: no cover
        raise TypeError(f"Expected UMCPosteriorResult, got {type(result).__name__}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    means = np.asarray(result.posterior_means_per_read)
    bonuses = np.asarray(result.bonus_per_read)
    obs_implied = means - bonuses  # reconstruct observed from per-read means

    bins = np.logspace(np.log10(100), np.log10(50_000), 50)
    ax1.hist(
        obs_implied,
        bins=bins,
        alpha=0.65,
        color="#ff7f0e",
        label=f"Observed (n={result.n_umc_reads:,})",
        edgecolor="none",
    )
    ax1.hist(
        means, bins=bins, alpha=0.65, color="#2ca02c", label="Posterior E[true]", edgecolor="none"
    )
    ax1.set_xscale("log")
    ax1.set_xlabel("Read length (bp)")
    ax1.set_ylabel("UMC read count")
    ax1.set_title("Observed vs posterior")
    ax1.legend(frameon=False, fontsize=9)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ax2.hist(bonuses, bins=40, color="#1f77b4", alpha=0.7, edgecolor="none")
    ax2.axvline(
        result.posterior_bonus_mean,
        color="black",
        linestyle="--",
        label=f"mean = {result.posterior_bonus_mean:,.0f} bp",
    )
    ax2.set_xlabel("Bonus length per read (bp)")
    ax2.set_ylabel("UMC read count")
    ax2.set_title("Per-read posterior bonus")
    ax2.legend(frameon=False, fontsize=9)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.suptitle(title or "UMC posterior length — adaptive-sampling truncation correction")
    fig.tight_layout()
    return fig


def plot_signal_trace(result, *, title: str | None = None):
    """Line plot of raw signal current over time for a single read."""
    import matplotlib.pyplot as plt
    import numpy as np

    from ..analyze.signal_trace import SignalTraceResult

    if not isinstance(result, SignalTraceResult):  # pragma: no cover
        raise TypeError(f"Expected SignalTraceResult, got {type(result).__name__}")

    fig, ax = plt.subplots(figsize=(10, 3.5))
    t = np.arange(result.n_samples) / max(result.samples_per_second, 1)
    colour = _PALETTE.get(result.end_reason, "#444444")
    ax.plot(t, result.signal, color=colour, linewidth=0.4)
    ax.set_xlabel("Time (sec)")
    ax.set_ylabel("Raw current (ADC)")
    er_label = result.end_reason_short or result.end_reason
    ax.set_title(title or f"Signal trace — {result.read_id[:8]}…  end_reason={er_label}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def plot_atlas_summary(
    result: AtlasResult,
    *,
    title: str | None = None,
    min_stratum_size: int = 3,
) -> Figure:
    """Fig8-style cross-run end-reason summary — per-stratum stacked bars.

    Each stratum (flowcell × chemistry × adaptive) becomes a position on the
    x-axis; for each end_reason metric we draw a grouped bar showing the
    mean percentage with a std error bar. Strata with fewer than
    `min_stratum_size` runs render with reduced opacity + hash pattern to
    flag low-confidence aggregates.

    Empty atlas → friendly placeholder figure (1 axis, no data).
    """
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(11, 5))

    if not result.per_stratum:
        ax.text(
            0.5,
            0.5,
            "Atlas is empty.\n\n"
            "Populate the qc_baseline store by running\n"
            "`ont-end-reason analyze distribution` on lab\n"
            "experiments, or `scripts/atlas_backfill.py` to\n"
            "seed from the ONT registry.",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=11,
            color="#444444",
        )
        ax.set_axis_off()
        ax.set_title(title or "Cross-run end-reason atlas (no data)")
        fig.tight_layout()
        return fig

    # Order strata by n_runs descending (most-data first).
    strata = sorted(result.per_stratum, key=lambda s: -s.n_runs)
    stratum_labels = [" / ".join(s.stratum) if s.stratum else "(unstratified)" for s in strata]

    metrics = [m for m in END_REASON_METRIC_KEYS if m in _ATLAS_METRIC_LABELS]
    n_metrics = len(metrics)
    x = np.arange(len(strata))
    bar_width = 0.8 / n_metrics

    for i, metric in enumerate(metrics):
        means = []
        stds = []
        alphas = []
        hatches = []
        for s in strata:
            stat = s.metric_stats.get(metric, {})
            means.append(stat.get("mean", 0.0))
            stds.append(stat.get("std", 0.0))
            low_conf = s.n_runs < min_stratum_size
            alphas.append(0.45 if low_conf else 0.9)
            hatches.append("//" if low_conf else "")
        offsets = x + (i - n_metrics / 2 + 0.5) * bar_width
        # Draw each bar individually so per-bar alpha + hatch apply.
        for xi, mean, std, alpha, hatch in zip(offsets, means, stds, alphas, hatches, strict=True):
            ax.bar(
                xi,
                mean,
                width=bar_width,
                color=_ATLAS_METRIC_COLORS[metric],
                edgecolor="black",
                linewidth=0.4,
                alpha=alpha,
                hatch=hatch,
                yerr=std,
                capsize=2,
                error_kw={"elinewidth": 0.6, "ecolor": "#333333"},
                label=_ATLAS_METRIC_LABELS[metric] if xi == offsets[0] else None,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{lbl}\n(n={s.n_runs})" for lbl, s in zip(stratum_labels, strata, strict=True)],
        rotation=30,
        ha="right",
        fontsize=8,
    )
    ax.set_ylabel("Mean % of reads (± std)")
    ax.set_xlabel("Stratum")
    n_strata = len(strata)
    default_title = (
        f"Cross-run end-reason atlas — "
        f"{result.n_internal} internal + {result.n_external} external runs "
        f"across {n_strata} strata"
    )
    ax.set_title(title or default_title)
    ax.legend(
        loc="upper right",
        frameon=False,
        fontsize=9,
        ncol=n_metrics,
        title="end_reason",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    return fig
