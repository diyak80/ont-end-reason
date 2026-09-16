"""ont-end-reason — command-line interface.

Single click-based dispatcher with one subcommand per capability. Heavy
subpackages (pysam, plotly, matplotlib) are lazy-imported inside the command
bodies so cold start is fast and `--help` works without optional dependencies.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import NoReturn

import click
import structlog

from . import __version__
from .codes import CODES, FAILED, RECOMMENDED_KEEP, TRUNCATED
from .errors import OntEndReasonError


def _configure_logging(*, debug: bool, quiet: bool) -> None:
    import logging

    level = logging.DEBUG if debug else logging.WARNING if quiet else logging.INFO
    logging.basicConfig(format="%(message)s", level=level, stream=sys.stderr)
    structlog.configure(
        processors=[
            structlog.processors.add_log_level,
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.dev.ConsoleRenderer(),
        ],
        wrapper_class=structlog.make_filtering_bound_logger(level),
    )


def _format_percentage(value: float | None) -> str:
    return "unavailable" if value is None else f"{value:.2f}"


def _die(message: str) -> NoReturn:
    click.echo(f"[error] {message}", err=True)
    sys.exit(1)


@click.group(
    context_settings={"help_option_names": ["-h", "--help"]},
    invoke_without_command=True,
)
@click.version_option(__version__, prog_name="ont-end-reason")
@click.option("--debug", is_flag=True, help="Enable DEBUG-level structured logging.")
@click.option("--quiet", is_flag=True, help="Suppress INFO logs; only errors.")
@click.option(
    "--strict",
    is_flag=True,
    help="Escalate warnings to errors (useful in CI).",
)
@click.pass_context
def main(ctx: click.Context, debug: bool, quiet: bool, strict: bool) -> None:
    """Comprehensive CLI for Oxford Nanopore end_reason analysis."""
    if debug and quiet:
        _die("--debug and --quiet are mutually exclusive")
    _configure_logging(debug=debug, quiet=quiet)
    ctx.ensure_object(dict)
    ctx.obj["strict"] = strict
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


# ───────────────────────── discovery ─────────────────────────


@main.command()
@click.argument("path", type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option("--recursive/--no-recursive", default=True, show_default=True)
@click.option(
    "--manifest",
    "manifest_out",
    type=click.Path(file_okay=True, dir_okay=False),
    default=None,
    help="Write the manifest as JSON to this path.",
)
def discover(path: str, recursive: bool, manifest_out: str | None) -> None:
    """Walk PATH and inventory POD5/Fast5/sequencing_summary/BAM/FASTQ files."""
    from .io.discover import discover as do_discover

    try:
        manifest = do_discover(path, recursive=recursive)
    except OntEndReasonError as exc:
        _die(str(exc))
    click.echo(
        f"Found {manifest.total_files()} files "
        f"({manifest.total_size_gb():.2f} GB) under {manifest.root}"
    )
    click.echo(f"  POD5:      {len(manifest.pod5)}")
    click.echo(f"  Fast5:     {len(manifest.fast5)}")
    click.echo(f"  Summaries: {len(manifest.summaries)}")
    click.echo(f"  BAMs:      {len(manifest.bams)}")
    click.echo(f"  FASTQs:    {len(manifest.fastqs)}")
    if manifest_out:
        manifest.to_json(manifest_out)
        click.echo(f"Manifest written: {manifest_out}")


# ───────────────────────── filter operations ─────────────────────────


@main.command()
@click.option("--summary", required=True, type=click.Path(exists=True))
@click.option("--bam", required=True, type=click.Path(exists=True))
@click.option("--out", required=True, type=click.Path())
@click.option("--tag-name", default="ER", show_default=True)
def tag(summary: str, bam: str, out: str, tag_name: str) -> None:
    """Tag BAM reads with end_reason from sequencing_summary.txt."""
    from .filter.tag import tag_bam

    try:
        result = tag_bam(summary, bam, out, tag_name=tag_name)
    except OntEndReasonError as exc:
        _die(str(exc))
    click.echo(
        f"Tagged {result.tagged_reads:,} / {result.input_reads:,} reads "
        f"(missing: {result.missing_reads:,}); wrote {out}"
    )


@main.command()
@click.option("--bam", required=True, type=click.Path(exists=True))
@click.option("--out", required=True, type=click.Path())
@click.option(
    "--keep",
    required=True,
    help="End reason codes to keep, comma-separated (e.g. SP or SP,UMC).",
)
@click.option("--tag-name", default="ER", show_default=True)
@click.option("--threads", default=1, type=int, show_default=True)
@click.option(
    "--shard-size",
    default=100_000,
    type=int,
    show_default=True,
    help="Target reads per shard when threads >= 2 (parallel mode).",
)
def filter(bam: str, out: str, keep: str, tag_name: str, threads: int, shard_size: int) -> None:
    """Filter a tagged BAM by end_reason."""
    from .filter.filter import filter_bam

    try:
        result = filter_bam(
            bam, out, keep, tag_name=tag_name, threads=threads, shard_size=shard_size
        )
    except OntEndReasonError as exc:
        _die(str(exc))
    click.echo(
        f"Kept {result.kept_reads:,} / {result.input_reads:,} reads "
        f"(dropped: {result.dropped_reads:,}); wrote {out}"
    )


@main.command(name="export-fastq")
@click.option("--bam", required=True, type=click.Path(exists=True))
@click.option("--fastq", required=True, type=click.Path())
@click.option("--compress", is_flag=True, help="gzip the output.")
def export_fastq(bam: str, fastq: str, compress: bool) -> None:
    """Export a (filtered) BAM to FASTQ for NanoPack tools."""
    from .filter.export import export_fastq as do_export

    try:
        result = do_export(bam, fastq, compress=compress)
    except OntEndReasonError as exc:
        _die(str(exc))
    click.echo(
        f"Wrote {result.reads_written:,} reads ({result.bytes_written:,} bytes) "
        f"to {result.output_path}"
    )


# ───────────────────────── analyze (group) ─────────────────────────


@main.group()
def analyze() -> None:
    """Run a paper-grade analysis on POD5/Fast5/summary input."""


@analyze.command("distribution")
@click.argument("source", type=click.Path(exists=True))
@click.option("--quick", is_flag=True, help="Sample up to 10k reads.")
@click.option("--max-reads", default=10_000, type=int, show_default=True)
@click.option("--json", "json_out", type=click.Path(), default=None)
@click.option("--plot", "plot_out", type=click.Path(), default=None)
@click.option(
    "--baseline-store/--no-baseline-store",
    "baseline_store",
    default=True,
    show_default=True,
    help=(
        "Auto-populate ~/.ont-qc-baselines when SOURCE matches a registry "
        "experiment. Silent no-op when SOURCE is unknown to the registry."
    ),
)
def analyze_distribution(
    source: str,
    quick: bool,
    max_reads: int,
    json_out: str | None,
    plot_out: str | None,
    baseline_store: bool,
) -> None:
    """End-reason distribution + OK/CHECK/FAIL/UNKNOWN quality gate."""
    from .analyze.distribution import (
        distribution as do_distribution,
    )
    from .analyze.distribution import (
        maybe_store_baseline,
    )

    try:
        result = do_distribution(source, quick=quick, max_reads=max_reads)
    except OntEndReasonError as exc:
        _die(str(exc))

    click.echo(f"Total reads:  {result.total_reads:,}")
    click.echo(f"Status:       {result.quality_status}")
    click.echo(f"Signal+ %:    {_format_percentage(result.signal_positive_pct)}")
    click.echo(f"UMC %:        {_format_percentage(result.unblock_mux_pct)}")
    click.echo(f"DUMC %:       {_format_percentage(result.data_service_pct)}")
    click.echo(f"Interpretation: {result.interpretation}")

    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
        click.echo(f"JSON: {json_out}")

    if plot_out:
        from .viz.static import plot_distribution

        fig = plot_distribution(result)
        Path(plot_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_out, dpi=300, bbox_inches="tight")
        click.echo(f"Plot: {plot_out}")

    # Atlas auto-population: silent when source isn't a registry-known path.
    stored_id = maybe_store_baseline(result, source, write=baseline_store)
    if stored_id:
        click.echo(f"Baseline stored: {stored_id}")


# Scaffolded analyses — each calls its NotImplementedError counterpart and
# the CLI gets a clean error message.


def _scaffold_analyze(cli_name: str, module_name: str) -> click.Command:
    @click.command(cli_name)
    @click.argument("source", type=click.Path(exists=True))
    def cmd(source: str) -> None:
        from importlib import import_module

        mod = import_module(f"ont_end_reason.analyze.{module_name}")
        func = getattr(mod, module_name)
        try:
            func(source)
        except NotImplementedError as exc:
            click.echo(f"[v0.2.0-roadmap] {exc}", err=True)
            sys.exit(2)

    return cmd


@analyze.command("sma-metrics")
@click.argument("source", type=click.Path(exists=True))
@click.option("--json", "json_out", type=click.Path(), default=None)
def analyze_sma_metrics(source: str, json_out: str | None) -> None:
    """SMA per-read metrics (delegates to the smaseq-qc package if installed)."""
    from .analyze.sma_metrics import sma_metrics as do_sma

    result = do_sma(source)
    click.echo(f"Available:    {result.available}")
    click.echo(f"Delegated to: {result.delegated_to}")
    for n in result.notes:
        click.echo(f"  - {n}")
    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2))


@analyze.command("tables")
@click.argument("source", type=click.Path(exists=True))
@click.option("--name", required=True, type=click.Choice(["summary", "per_class", "quality"]))
@click.option(
    "--format",
    "fmt",
    default="markdown",
    type=click.Choice(["tsv", "csv", "markdown", "latex"]),
    show_default=True,
)
@click.option(
    "--out",
    "out_path",
    type=click.Path(),
    default=None,
    help="Write rendered table to this path; omit to print to stdout.",
)
def analyze_tables(source: str, name: str, fmt: str, out_path: str | None) -> None:
    """Generate a paper-table from an analysis result.

    Examples:
      ont-end-reason analyze tables <summary> --name summary --format markdown
      ont-end-reason analyze tables <summary> --name per_class --format tsv --out per_class.tsv
    """
    from .analyze.tables import generate_tables, render_table

    try:
        table = generate_tables(source, name=name, fmt=fmt)
    except OntEndReasonError as exc:
        _die(str(exc))
    text = render_table(table, out_path)
    if out_path:
        click.echo(f"Wrote {out_path}")
    else:
        click.echo(text)


@analyze.command("umc-posterior")
@click.argument("source", type=click.Path(exists=True))
@click.option(
    "--prior-class",
    default="signal_positive",
    show_default=True,
    help="End reason whose length distribution provides the prior.",
)
@click.option("--json", "json_out", type=click.Path(), default=None)
@click.option("--plot", "plot_out", type=click.Path(), default=None)
def analyze_umc_posterior(
    source: str,
    prior_class: str,
    json_out: str | None,
    plot_out: str | None,
) -> None:
    """Posterior over true UMC length given adaptive-sampling truncation.

    The paper's central novel analysis. Returns the expected "bonus" — i.e.
    how much sequence was actually present in the molecule beyond the
    adaptive-sampling cutoff.
    """
    from .analyze.umc_posterior import umc_posterior as do_umc

    try:
        result = do_umc(source, prior_class=prior_class)
    except OntEndReasonError as exc:
        _die(str(exc))

    click.echo(f"UMC reads:              {result.n_umc_reads:,}")
    click.echo(
        f"Prior class:            {result.prior_class}  "
        f"(log μ={result.prior_log_mu:.3f}, log σ={result.prior_log_sigma:.3f})"
    )
    click.echo(f"Observed mean length:   {result.observed_mean:>10.1f} bp")
    click.echo(f"Observed median:        {result.observed_median:>10.1f} bp")
    click.echo(f"Posterior expected mean:{result.posterior_expected_true_mean:>10.1f} bp")
    click.echo(f"Posterior expected med: {result.posterior_expected_true_median:>10.1f} bp")
    click.echo(f"Posterior bonus mean:   {result.posterior_bonus_mean:>10.1f} bp/read")
    click.echo(f"Posterior bonus total:  {result.posterior_bonus_total:>14,.0f} bp")
    ci_lo, ci_hi = result.credible_interval_95_per_read_mean
    click.echo(f"95% CI per-read (mean): [{ci_lo:.1f}, {ci_hi:.1f}]")

    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2))
        click.echo(f"JSON: {json_out}")

    if plot_out:
        from .viz.static import plot_umc_posterior

        fig = plot_umc_posterior(result)
        Path(plot_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_out, dpi=300, bbox_inches="tight")
        click.echo(f"Plot: {plot_out}")


@analyze.command("hypothesis")
@click.argument("source", type=click.Path(exists=True))
@click.option(
    "-a",
    "class_a",
    default="SP",
    show_default=True,
    help="First end_reason class (short code or full name).",
)
@click.option("-b", "class_b", default="UMC", show_default=True, help="Second end_reason class.")
@click.option(
    "--test", default="mann-whitney", type=click.Choice(["mann-whitney", "ks"]), show_default=True
)
@click.option(
    "--column",
    default="sequence_length_template",
    type=click.Choice(["sequence_length_template", "mean_qscore_template"]),
    show_default=True,
)
@click.option("--json", "json_out", type=click.Path(), default=None)
def analyze_hypothesis(
    source: str,
    class_a: str,
    class_b: str,
    test: str,
    column: str,
    json_out: str | None,
) -> None:
    """Mann-Whitney or KS test between two end_reason populations."""
    from .analyze.hypothesis import hypothesis as do_hypothesis

    try:
        result = do_hypothesis(source, a=class_a, b=class_b, test=test, column=column)
    except OntEndReasonError as exc:
        _die(str(exc))

    click.echo(f"Test:           {result.test}")
    click.echo(f"Column:         {result.column}")
    click.echo(f"Comparison:     {result.comparison[0]}  vs  {result.comparison[1]}")
    click.echo(f"n_a / n_b:      {result.n_a:,} / {result.n_b:,}")
    click.echo(f"Statistic:      {result.statistic:.4g}")
    click.echo(f"p-value:        {result.p_value:.4g}")
    click.echo(f"Cliff's Δ:      {result.effect_size:+.4f}")
    click.echo(f"Median a / b:   {result.median_a:.3f} / {result.median_b:.3f}")
    for note in result.notes:
        click.echo(f"Note: {note}")

    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2))
        click.echo(f"JSON: {json_out}")


@analyze.command("quality")
@click.argument("source", type=click.Path(exists=True))
@click.option(
    "--gmm-components",
    default=3,
    type=int,
    show_default=True,
    help="Max GMM components to try (BIC selects the best).",
)
@click.option("--json", "json_out", type=click.Path(), default=None)
@click.option("--plot", "plot_out", type=click.Path(), default=None)
def analyze_quality(
    source: str,
    gmm_components: int,
    json_out: str | None,
    plot_out: str | None,
) -> None:
    """Per-end_reason Q-score distribution + GMM fit (BIC-selected components)."""
    from .analyze.quality import quality as do_quality

    try:
        result = do_quality(source, gmm_components=gmm_components)
    except OntEndReasonError as exc:
        _die(str(exc))

    click.echo(f"Total reads: {result.total_reads:,}")
    click.echo(f"{'End reason':<35}{'n':>8}{'mean Q':>9}{'median':>9}{'k':>4}")
    for er, s in sorted(result.per_class.items(), key=lambda kv: -kv[1].n):
        click.echo(f"  {er:<33}{s.n:>8,d}{s.mean:>9.2f}{s.median:>9.2f}{s.gmm_chosen_k:>4d}")

    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2))
        click.echo(f"JSON: {json_out}")

    if plot_out:
        from .viz.static import plot_quality_violins

        fig = plot_quality_violins(result)
        Path(plot_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_out, dpi=300, bbox_inches="tight")
        click.echo(f"Plot: {plot_out}")


@analyze.command("temporal")
@click.argument("source", type=click.Path(exists=True))
@click.option("--bin-seconds", default=3600.0, type=float, show_default=True)
@click.option("--json", "json_out", type=click.Path(), default=None)
@click.option("--plot", "plot_out", type=click.Path(), default=None)
def analyze_temporal(
    source: str,
    bin_seconds: float,
    json_out: str | None,
    plot_out: str | None,
) -> None:
    """End-reason rates over time. Useful for flowcell-degradation diagnostics."""
    from .analyze.temporal import temporal as do_temporal

    try:
        result = do_temporal(source, bin_seconds=bin_seconds)
    except OntEndReasonError as exc:
        _die(str(exc))

    click.echo(f"Total reads: {result.total_reads:,}  bins: {len(result.bin_centers)}")
    click.echo(f"Bin width:   {result.bin_seconds / 3600:.1f}h")
    # Print SP and UMC fraction time series
    for er in ("signal_positive", "unblock_mux_change"):
        fracs = result.fractions_by_class.get(er, [])
        if not fracs:
            continue
        first, last = fracs[0] * 100, fracs[-1] * 100
        click.echo(f"  {er}: start={first:.1f}% end={last:.1f}% Δ={last - first:+.1f}pp")

    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2))
        click.echo(f"JSON: {json_out}")

    if plot_out:
        from .viz.static import plot_temporal

        fig = plot_temporal(result)
        Path(plot_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_out, dpi=300, bbox_inches="tight")
        click.echo(f"Plot: {plot_out}")


@analyze.command("length")
@click.argument("source", type=click.Path(exists=True))
@click.option("--json", "json_out", type=click.Path(), default=None)
@click.option("--plot", "plot_out", type=click.Path(), default=None)
def analyze_length(source: str, json_out: str | None, plot_out: str | None) -> None:
    """Per-end_reason length distribution: n, mean, median, percentiles, N50."""
    from .analyze.length import length as do_length

    try:
        result = do_length(source)
    except OntEndReasonError as exc:
        _die(str(exc))

    click.echo(f"Total reads:  {result.total_reads:,}")
    click.echo(f"{'End reason':<35}{'n':>8}{'median':>9}{'p95':>9}{'N50':>9}")
    for er, s in sorted(result.per_class.items(), key=lambda kv: -kv[1].n):
        click.echo(f"  {er:<33}{s.n:>8,d}{s.median:>9,.0f}{s.p95:>9,.0f}{s.n50:>9,d}")

    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2))
        click.echo(f"JSON: {json_out}")

    if plot_out:
        from .viz.static import plot_length_distribution

        fig = plot_length_distribution(result)
        Path(plot_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_out, dpi=300, bbox_inches="tight")
        click.echo(f"Plot: {plot_out}")


@analyze.command("signal-trace")
@click.argument("pod5", type=click.Path(exists=True))
@click.option("--read-id", required=True)
@click.option("--json", "json_out", type=click.Path(), default=None)
@click.option("--plot", "plot_out", type=click.Path(), default=None)
def analyze_signal_trace(
    pod5: str, read_id: str, json_out: str | None, plot_out: str | None
) -> None:
    """Extract and (optionally) visualise the raw signal trace for one read."""
    from .analyze.signal_trace import signal_trace

    try:
        result = signal_trace(pod5, read_id=read_id)
    except OntEndReasonError as exc:
        _die(str(exc))

    click.echo(f"read_id:              {result.read_id}")
    click.echo(f"end_reason:           {result.end_reason}  ({result.end_reason_short})")
    click.echo(f"samples:              {result.n_samples:,}")
    click.echo(f"sample_rate:          {result.samples_per_second:,} Hz")
    click.echo(f"duration:             {result.duration_seconds:.3f} sec")
    click.echo(f"source:               {result.source_file}")

    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2))
        click.echo(f"JSON: {json_out}")

    if plot_out:
        from .viz.static import plot_signal_trace

        fig = plot_signal_trace(result)
        Path(plot_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_out, dpi=300, bbox_inches="tight")
        click.echo(f"Plot: {plot_out}")


# ───────────────────────── figure (group) ─────────────────────────


@main.group()
def figure() -> None:
    """Reproduce a named paper figure (fig3, fig5, fig6, supplementary)."""


@figure.command("fig3")
@click.argument("source", type=click.Path(exists=True))
@click.option("--out", required=True, type=click.Path())
@click.option("--quick", is_flag=True)
def figure_fig3(source: str, out: str, quick: bool) -> None:
    """Paper Figure 3 — end-reason distribution bar chart."""
    from .figures.fig3_distribution import fig3_distribution

    try:
        path = fig3_distribution(source, out=out, quick=quick)
    except OntEndReasonError as exc:
        _die(str(exc))
    click.echo(f"Figure 3 written: {path}")


@figure.command("fig5")
@click.argument("source", type=click.Path(exists=True))
@click.option("--out", required=True, type=click.Path())
def figure_fig5(source: str, out: str) -> None:
    """Paper Figure 5 — Q-score violins per end_reason."""
    from .figures.fig5_violin import fig5_violin

    try:
        path = fig5_violin(source, out=out)
    except OntEndReasonError as exc:
        _die(str(exc))
    click.echo(f"Figure 5 written: {path}")


@figure.command("fig6")
@click.argument("source", type=click.Path(exists=True))
@click.option("--out", required=True, type=click.Path())
def figure_fig6(source: str, out: str) -> None:
    """Paper Figure 6 — UMC posterior diagram (observed vs inferred true length)."""
    from .figures.fig6_conceptual import fig6_conceptual

    try:
        path = fig6_conceptual(source, out=out)
    except OntEndReasonError as exc:
        _die(str(exc))
    click.echo(f"Figure 6 written: {path}")


@figure.command("supplementary")
@click.argument("source", type=click.Path(exists=True))
@click.option("--out", required=True, type=click.Path())
def figure_supplementary(source: str, out: str) -> None:
    """[v0.2.0] Paper supplementary figures."""
    from .figures.supplementary import supplementary

    try:
        supplementary(source, out=out)
    except NotImplementedError as exc:
        click.echo(f"[v0.2.0-roadmap] {exc}", err=True)
        sys.exit(2)


@figure.command("atlas")
@click.option("--out", required=True, type=click.Path())
def figure_atlas(out: str) -> None:
    """Paper Figure 8 — cross-run end-reason atlas summary.

    Composes `analyze.atlas` (cross-run aggregation across qc_baseline +
    external peers) and `viz.static.plot_atlas_summary`. Works on an empty
    store — renders a friendly placeholder, never crashes.
    """
    from .figures.atlas import atlas as fig_atlas

    try:
        path = fig_atlas(out=out)
    except OntEndReasonError as exc:
        _die(str(exc))
    click.echo(f"Atlas figure written: {path}")


# ───────────────────────── report (group) ─────────────────────────


@main.group()
def report() -> None:
    """Build composed multi-analysis reports."""


@report.command("interactive")
@click.argument("source", type=click.Path(exists=True))
@click.option("--out", required=True, type=click.Path())
@click.option("--quick", is_flag=True)
def report_interactive(source: str, out: str, quick: bool) -> None:
    """Self-contained interactive HTML report with Plotly figures."""
    from .report.html import build_html_report

    try:
        result = build_html_report(source, output_path=out, quick=quick)
    except OntEndReasonError as exc:
        _die(str(exc))
    click.echo(
        f"HTML report ({result.n_reads:,} reads, "
        f"sections: {', '.join(result.sections)}): {result.output_path}"
    )


@report.command("static")
@click.argument("source", type=click.Path(exists=True))
@click.option("--out", required=True, type=click.Path())
def report_static(source: str, out: str) -> None:
    """[v0.2.0] PDF report composing all analyses."""
    click.echo(
        "[v0.2.0-roadmap] Static PDF report is scheduled for v0.2.0. "
        "Use `report interactive` for the current functionality.",
        err=True,
    )
    sys.exit(2)


# ───────────────────────── utility ─────────────────────────


@main.command()
def codes() -> None:
    """Print the end_reason abbreviation taxonomy."""
    click.echo("End reason codes (canonical lab taxonomy)\n")
    click.echo(f"{'Code':<6}{'Full name':<40}{'Class':<15}")
    click.echo("─" * 61)
    for full, short in CODES.items():
        if short in RECOMMENDED_KEEP:
            cls = "keep"
        elif short in TRUNCATED:
            cls = "truncated"
        elif short in FAILED:
            cls = "failed"
        else:
            cls = "unknown"
        click.echo(f"{short:<6}{full:<40}{cls:<15}")


@main.command()
def schema() -> None:
    """Print the canonical sequencing_summary column set."""
    from .io.schema import RECOMMENDED_COLUMNS, REQUIRED_COLUMNS

    click.echo("sequencing_summary.txt schema (required + recommended columns)\n")
    click.echo("Required:")
    for c in sorted(REQUIRED_COLUMNS):
        click.echo(f"  - {c}")
    click.echo("\nRecommended:")
    for c in sorted(RECOMMENDED_COLUMNS):
        click.echo(f"  - {c}")


@main.command()
@click.option("--summary", required=True, type=click.Path(exists=True))
@click.option("--out", type=click.Path(), default=None)
@click.option("--json", "json_out", type=click.Path(), default=None)
def stats(summary: str, out: str | None, json_out: str | None) -> None:
    """Streaming QC for PromethION-scale sequencing_summary.txt.

    v0.1.0 ships a minimal implementation; the canonical version with
    additional metrics is scheduled for v0.2.0.
    """
    from .analyze.distribution import distribution as do_distribution

    try:
        result = do_distribution(summary)
    except OntEndReasonError as exc:
        _die(str(exc))

    text = (
        f"sequencing_summary stats — {summary}\n"
        f"  total_reads:        {result.total_reads:,}\n"
        f"  quality_status:     {result.quality_status}\n"
        f"  signal_positive_%:  {_format_percentage(result.signal_positive_pct)}\n"
        f"  unblock_mux_%:      {_format_percentage(result.unblock_mux_pct)}\n"
        f"  data_service_%:     {_format_percentage(result.data_service_pct)}\n"
    )
    if out:
        Path(out).write_text(text)
        click.echo(f"Wrote {out}")
    else:
        click.echo(text)

    if json_out:
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
        click.echo(f"Wrote {json_out}")


@analyze.command("atlas")
@click.option("--json", "json_out", type=click.Path(), default=None)
@click.option("--plot", "plot_out", type=click.Path(), default=None)
@click.option(
    "--include-internal/--no-include-internal",
    default=True,
    show_default=True,
    help="Pull internal-lab peers from ~/.ont-qc-baselines.",
)
@click.option(
    "--include-external/--no-include-external",
    default=True,
    show_default=True,
    help="Pull external (public ONT) peers from the fingerprint cache.",
)
@click.option(
    "--strata",
    default="flowcell_type,chemistry,adaptive_sampling",
    show_default=True,
    help="Comma-separated metadata fields to stratify peer cohorts by.",
)
@click.option(
    "--z-threshold",
    default=2.0,
    type=float,
    show_default=True,
    help="Composite anomaly score >= this flags an outlier.",
)
def analyze_atlas(
    json_out: str | None,
    plot_out: str | None,
    include_internal: bool,
    include_external: bool,
    strata: str,
    z_threshold: float,
) -> None:
    """Cross-run end-reason atlas — aggregate peers, flag outliers.

    Reads from the qc_baseline store (internal lab cohort) and the external
    peer fingerprint cache (public ONT datasets), stratifies by metadata
    fields, computes per-stratum baselines, and flags runs whose composite
    anomaly score exceeds --z-threshold. Empty stores degrade gracefully.
    """
    from .analyze.atlas import atlas

    strata_tuple = tuple(s.strip() for s in strata.split(",") if s.strip())
    result = atlas(
        include_internal=include_internal,
        include_external=include_external,
        strata=strata_tuple,
        z_threshold=z_threshold,
    )

    click.echo(
        f"Atlas: {result.n_internal} internal + {result.n_external} external "
        f"runs across {len(result.per_stratum)} strata"
    )
    if result.outliers:
        click.echo("Outliers flagged (top 5 by anomaly_score):")
        for o in result.outliers[:5]:
            click.echo(
                f"  [{o.source:<8}] {o.experiment_id:<40} "
                f"stratum={o.stratum}  score={o.anomaly_score:.2f}"
            )
    else:
        click.echo("No outliers flagged.")

    if result.interpretation:
        click.echo(f"\n  {result.interpretation}")

    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(result.to_dict(), indent=2))
        click.echo(f"\nJSON: {json_out}")

    if plot_out:
        from .viz.static import plot_atlas_summary

        try:
            fig = plot_atlas_summary(result)
        except ImportError:
            _die("matplotlib not available — install ont-end-reason[viz] to plot")
        Path(plot_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_out, dpi=200, bbox_inches="tight")
        import matplotlib.pyplot as plt

        plt.close(fig)
        click.echo(f"Plot: {plot_out}")


if __name__ == "__main__":
    main()
