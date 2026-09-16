<div align="center">

# 🧬 ont-end-reason

**Comprehensive CLI for Oxford Nanopore `end_reason` analysis.**
Discover, tag, filter, analyse, and visualise read-termination patterns.

[![CI](https://github.com/Single-Molecule-Sequencing/ont-end-reason/actions/workflows/ci.yml/badge.svg)](https://github.com/Single-Molecule-Sequencing/ont-end-reason/actions/workflows/ci.yml)
[![Docs](https://github.com/Single-Molecule-Sequencing/ont-end-reason/actions/workflows/docs.yml/badge.svg)](https://github.com/Single-Molecule-Sequencing/ont-end-reason/actions/workflows/docs.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-143%20passing-brightgreen.svg)](#testing)
[![Coverage](https://img.shields.io/badge/coverage-63%25-brightgreen.svg)](#testing)
[![Version](https://img.shields.io/badge/version-0.2.0a1-blue.svg)](CHANGELOG.md)

### 🚀 [**→ Interactive dashboard & tutorials**](https://rdlu0053.tail46dbe4.ts.net:8464/artifacts/mirror/ont-end-reason/index.html)

Companion to the [end-reason paper](https://github.com/Single-Molecule-Sequencing/end_reason_6_5_26).

</div>

---

## Table of contents

- [Why this tool](#why-this-tool)
- [Install](#install)
- [Quickstart](#quickstart)
- [The headline result](#the-headline-result)
- [CLI surface](#cli-surface)
- [Python API](#python-api)
- [End_reason taxonomy](#end_reason-taxonomy)
- [How the UMC posterior works](#how-the-umc-posterior-works)
- [Testing](#testing)
- [Lab infrastructure integration](#lab-infrastructure-integration)
- [Status / roadmap](#status--roadmap)
- [Citing](#citing)
- [License](#license)

---

## Why this tool

Oxford Nanopore sequencers tag every read with an `end_reason` explaining
why sequencing stopped. **A read can have high base quality (Q>20) and
still be truncated or rejected by adaptive sampling** — filtering by
Q-score alone is not enough for accurate downstream analysis.

`ont-end-reason` unifies the eight published analyses from the end-reason
paper into a single PyPI-installable CLI, including the paper's novel
*posterior length model* for adaptive-sampling-truncated reads.

Before this tool, the analyses lived in scattered scripts inside
[End_Reason_Manuscript/pipeline/bin/](https://github.com/Single-Molecule-Sequencing/end_reason_6_5_26)
(now archived). Every script was promoted to this repo with provenance
headers crediting commit `b47166a` of the source. The package is the
canonical implementation going forward.

---

## Install

```bash
# Static figures only (matplotlib)
pip install ont-end-reason

# + Plotly for interactive HTML reports
pip install "ont-end-reason[interactive]"

# Development (from source)
git clone https://github.com/Single-Molecule-Sequencing/ont-end-reason.git
cd ont-end-reason
pip install -e ".[dev,interactive]"
```

Python 3.10+ required. Tested on Linux + macOS, Python 3.10 through 3.13.

---

## Quickstart

Five commands cover the canonical pipeline:

```bash
# 1. Inventory what's in a sequencing-run directory
ont-end-reason discover /path/to/run --manifest run.json

# 2. Tag a BAM with end_reason from sequencing_summary.txt
ont-end-reason tag --summary sequencing_summary.txt \
                   --bam aligned.bam --out tagged.bam

# 3. Filter to complete reads only (signal_positive)
ont-end-reason filter --bam tagged.bam --keep SP --out complete.bam

# 4. Run the paper's central novel analysis
ont-end-reason analyze umc-posterior sequencing_summary.txt --plot umc.pdf

# 5. Build a self-contained 6-section HTML report
ont-end-reason report interactive sequencing_summary.txt --out report.html
```

→ Full walkthrough with live charts on the [**dashboard**](https://rdlu0053.tail46dbe4.ts.net:8464/artifacts/mirror/ont-end-reason/index.html).

---

## The headline result

On the [synthetic 5000-read test fixture](tests/fixtures/sequencing_summary_synthetic.txt):

```
$ ont-end-reason analyze umc-posterior tests/fixtures/sequencing_summary_synthetic.txt
UMC reads:              600
Prior class:            signal_positive  (log μ=8.488, log σ=0.600)
Observed mean length:        926.2 bp
Posterior expected mean:    5868.1 bp
Posterior bonus mean:       4941.9 bp/read
Posterior bonus total:       2,965,111 bp     ← ~3 Mb of unobserved sequence
```

**Adaptive-sampling truncation hides ~5× more sequence than the observed read length suggests.** Scaled to a real PromethION run with millions of UMC reads, the recovered-sequence estimate grows linearly. This is exactly what the paper's central analysis is for — and the tool surfaces it as one command on any `sequencing_summary.txt`.

![UMC posterior](docs/figures/umc_posterior.png)

---

## CLI surface

### Discovery + filter operations

| Command | Purpose |
|---|---|
| `ont-end-reason discover <path>` | Walk a directory, inventory POD5 / Fast5 / summary / BAM / FASTQ files |
| `ont-end-reason tag` | Add end_reason tag to BAM reads from sequencing_summary.txt |
| `ont-end-reason filter` | Keep / drop BAM reads by end_reason short code (parallel sharded, `--threads N`) |
| `ont-end-reason export-fastq` | Convert filtered BAM → FASTQ for NanoPack tools |
| `ont-end-reason stats` | Streaming QC summary from sequencing_summary.txt |

### Analysis (9 subcommands)

| Command | What it does |
|---|---|
| `ont-end-reason analyze distribution` | Per-end_reason counts + OK/CHECK/FAIL quality gate |
| `ont-end-reason analyze length` | Length distributions per end_reason (N50, percentiles) |
| `ont-end-reason analyze quality` | Q-score distributions with Gaussian Mixture Model fit |
| `ont-end-reason analyze temporal` | End_reason rates over sequencing-run time |
| `ont-end-reason analyze hypothesis` | Mann-Whitney U / KS tests with Cliff's Δ effect size |
| **`ont-end-reason analyze umc-posterior`** ⭐ | Bayesian posterior on truncated UMC length (paper's central analysis) |
| `ont-end-reason analyze signal-trace` | Raw POD5 current trace extraction for a single read |
| `ont-end-reason analyze sma-metrics` | Optional bridge to the `smaseq-qc` package |
| `ont-end-reason analyze tables` | Generate summary/per-class/quality tables (TSV/CSV/md/LaTeX) |

### Paper-figure reproducers + reports

| Command | Output |
|---|---|
| `ont-end-reason figure fig3 <source>` | Paper Figure 3 — distribution bar chart |
| `ont-end-reason figure fig5 <source>` | Paper Figure 5 — Q-score violins |
| `ont-end-reason figure fig6 <source>` | Paper Figure 6 — UMC posterior diagram |
| `ont-end-reason report interactive` | 6-section self-contained HTML report with embedded Plotly |
| `ont-end-reason report static` | Paginated PDF report (v0.3.0 roadmap) |

Run `ont-end-reason <cmd> --help` for full flag documentation. Examples and screenshots: [dashboard](https://rdlu0053.tail46dbe4.ts.net:8464/artifacts/mirror/ont-end-reason/index.html).

---

## Python API

Every CLI subcommand has a public Python API equivalent. Functions return typed `dataclass`es so callers can compose, persist, or pipe results without re-parsing CLI output:

```python
from ont_end_reason import discover, classify
from ont_end_reason.analyze.distribution import distribution
from ont_end_reason.analyze.umc_posterior import umc_posterior
from ont_end_reason.viz.static import plot_umc_posterior

# Discovery → Manifest
manifest = discover("/path/to/sequencing_run")
print(f"Found {manifest.total_files()} files")

# Analysis → typed result
result = umc_posterior("sequencing_summary.txt")
print(f"Posterior bonus total: {result.posterior_bonus_total:,.0f} bp")

# Visualisation → matplotlib Figure
fig = plot_umc_posterior(result)
fig.savefig("umc.pdf")
```

Each analysis result has a `.to_dict()` for JSON serialisation and roundtrip.

---

## End_reason taxonomy

The lab's canonical 7-class taxonomy. Print from the CLI any time with `ont-end-reason codes`:

| Code | Full name | Class | Action |
|---|---|---|---|
| `SP` | signal_positive | **keep** | Complete read — always keep |
| `UMC` | unblock_mux_change | truncated | Filter unless studying artifacts |
| `MC` | mux_change | truncated | Filter |
| `DUMC` | data_service_unblock_mux_change | truncated | Filter (software-triggered) |
| `PART` | partial | truncated | Filter |
| `SN` | signal_negative | **failed** | Always filter |
| `UNK` | unknown | unknown | Investigate distribution |

`--keep SP` is the canonical recommendation (Table 1 of end-reason-paper). Use `--keep SP,UMC` to retain truncated reads for artifact studies.

---

## How the UMC posterior works

The paper's novel analytic contribution, in one paragraph:

Given an observed UMC read of length `o`, the molecule's true length `L` is *unknown but at least `o`* (it was truncated, not foreshortened). Fitting a lognormal prior `L ~ Lognormal(μ, σ²)` to `signal_positive` reads gives the prior on what completed reads look like; the posterior on a UMC read's true length is then the prior left-truncated at the observation:

```
P(L | L ≥ o)  ∝  Lognormal(L; μ, σ²) · 𝟙[L ≥ o]
```

The truncated mean has a closed form via the normal CDF's Mills ratio:

```
E[L | L ≥ o]  =  exp(μ + σ²/2) · Φ(σ - z) / (1 - F(o))    where  z = (log o − μ)/σ
```

Implementation: `scipy.stats.lognorm`, vectorised over all UMC reads. O(n).
Aggregated, this is the paper's headline "sequence lost to adaptive sampling" estimate — runnable on any sequencing_summary.txt with one command.

---

## Testing

```bash
pytest                       # 175 tests, ~40s
pytest --cov=ont_end_reason  # with coverage (currently 71%)
ruff check .                 # lint
mypy src/ont_end_reason      # type-check
```

Coverage gate is 69% in CI (1 pp below actual, to absorb cross-platform variance).

Tests run against:
- **Synthetic fixture** ([5000 reads, deterministic distributions](tests/fixtures/sequencing_summary_synthetic.txt)) for every analysis
- **Hypothesis property tests** for the SP/UMC/MC taxonomy (round-trips, classification disjointness)
- **CliRunner integration tests** for every subcommand's `--help` and dispatch
- **Real-data smoke** against the AWG074 MinION run (1,571 tagged reads → 1,451 SP / 89 UMC / 31 SN)

---

## Performance

The `filter` subcommand runs sequential by default; `--threads N` engages a parallel
sharded path for inputs above ~50k reads:

- **Shard boundaries** are placed by virtual offset (`bam.tell()`) during a single
  sequential scan — workers `seek()` directly to their slice, avoiding the original
  port's O(N²/2) linear-skip pattern.
- **Shard merge** uses `pysam.cat`, which splices BGZF blocks without re-decompressing.
- **Auto-fallback** to sequential below `MIN_READS_FOR_PARALLEL` (50k reads) — worker-pool
  setup outweighs the gain on small inputs.
- **Bit-identical output**: enforced by an integration test that compares kept-read
  sets across both paths on every CI run.

Synthetic micro-bench (`bench/bench_parallel_filter.py`, ONT-shaped 2 kb reads,
8-core dev machine):

| n_reads | threads | shard_size | shards | seq_s | par_s | speedup |
|---:|:---:|---:|:---:|---:|---:|:---:|
| 20,000 | 4 | 2,000 | 2 | 0.02 | 0.05 | 0.45× *(below threshold)* |
| 100,000 | 4 | 12,500 | 9 | 0.78 | 0.72 | 1.08× |
| 300,000 | 4 | 37,500 | 9 | 2.14 | 1.91 | 1.12× |

Speedup is intentionally modest at this workload — pysam already pipelines BGZF
decompression internally, so the worker pool only parallelizes tag-lookup + write.
Larger gains expected on multi-GB real ONT BAMs where per-record CPU cost is higher.

---

## Cross-run atlas

The `analyze atlas` subcommand answers a question single-run analysis can't:
**"is THIS run normal relative to all the lab's prior runs (and the public ONT
community)?"**

```bash
# Aggregate across the qc_baseline store + external peer cache
ont-end-reason analyze atlas --json atlas.json --plot atlas.png

# Stratify on different metadata dimensions
ont-end-reason analyze atlas --strata flowcell_type,basecaller_model

# Tighten outlier flagging (default z >= 2.0)
ont-end-reason analyze atlas --z-threshold 3.0
```

**Data sources:**
- **Internal lab peers** — auto-populated into `~/.ont-qc-baselines/` by every
  `ont-end-reason analyze distribution` invocation (see `--baseline-store`).
- **External peers** — public ONT datasets (GIAB, hereditary-cancer ONT Open
  Data) cached as Parquet fingerprints at `~/.ont-qc-baselines/external_peers/`,
  refreshed by the lab's `/ont-public-data` skill.

**Backfill:** one-time `scripts/atlas_backfill.py --dry-run` lists every
registry experiment eligible for ingest; drop the `--dry-run` to run them all.

**Output shape:** AtlasResult JSON with `per_stratum` (mean/median/std/min/max
for all 5 end-reason metrics per stratum), `outliers` (runs with composite
`anomaly_score = max(|z_i|) >= --z-threshold`), and a human-readable
`interpretation`. Designed-for: paper figure regenerators, dashboard panels,
batch QC gates.

Spec: [docs/superpowers/specs/2026-05-12-end-reason-atlas-design.md](docs/superpowers/specs/2026-05-12-end-reason-atlas-design.md)

---

## Lab infrastructure integration

`ont-end-reason` is part of the Single-Molecule-Sequencing org's analytic toolchain:

| Repo | How it integrates |
|---|---|
| [end-reason-paper](https://github.com/Single-Molecule-Sequencing/end_reason_6_5_26) | Companion paper. Claim atoms (`results.alignment_rate_filtered`, `results.snv_f1_filtered`, etc.) pin to this tool for reproducibility. |
| [ont-ecosystem](https://github.com/Single-Molecule-Sequencing/ont-ecosystem) | Lab Claude Code skills `/end-reason` and `/end-reason-filter` will become thin wrappers that `pip install ont-end-reason` (tracked in [issue #6](https://github.com/Single-Molecule-Sequencing/ont-end-reason/issues/6)). |
| [lab-onboarding](https://github.com/Single-Molecule-Sequencing/lab-system) | Bundled in the canonical lab-repo manifest. Cloned automatically by `bash wsl/bootstrap.sh` on every new lab device. |
| [End_Reason_Manuscript](https://github.com/Single-Molecule-Sequencing/end_reason_6_5_26) | **Archived.** Each script in this repo carries a provenance header crediting commit `b47166a` of that source. |
| [smaseq-qc](https://github.com/Single-Molecule-Sequencing/smaseq-qc) | Optional dependency for `analyze sma-metrics`. Tool detects-and-skips when missing. |

---

## Status / roadmap

**Current: v0.2.0a1 (alpha)**

- ✅ 9 analysis subcommands fully implemented
- ✅ Bayesian posterior model for UMC truncation (paper's central novel analysis)
- ✅ Interactive HTML reports with embedded Plotly
- ✅ 143 tests, CI matrix on Python 3.10–3.13 × Ubuntu/macOS
- ✅ [Interactive dashboard](https://rdlu0053.tail46dbe4.ts.net:8464/artifacts/mirror/ont-end-reason/index.html) with live examples
- 🚧 Reproducibility CI against end-reason-paper claim atoms ([#4](https://github.com/Single-Molecule-Sequencing/ont-end-reason/issues/4))
- 🚧 Parallel sharded BAM filtering ([#5](https://github.com/Single-Molecule-Sequencing/ont-end-reason/issues/5))
- 🚧 Lab-skill thin-wrap migration after PyPI release ([#6](https://github.com/Single-Molecule-Sequencing/ont-end-reason/issues/6))
- ⏳ conda-forge feedstock (post-v0.1.0 PyPI)

See [`CHANGELOG.md`](CHANGELOG.md) for per-release detail and [open issues](https://github.com/Single-Molecule-Sequencing/ont-end-reason/issues) for roadmap items.

---

## Citing

If you use `ont-end-reason` in published work, please cite the companion paper:

```
Athey BD et al. (in preparation). End reason filtering for accurate analysis
of Oxford Nanopore sequencing data. Single-Molecule-Sequencing Lab,
University of Michigan.
https://github.com/Single-Molecule-Sequencing/end_reason_6_5_26
```

Machine-readable citation metadata is in [`CITATION.cff`](CITATION.cff).

---

## License

MIT — see [`LICENSE`](LICENSE).

---

<div align="center">

Built by the [Athey Lab](https://github.com/Single-Molecule-Sequencing) at the University of Michigan.

[Dashboard](https://rdlu0053.tail46dbe4.ts.net:8464/artifacts/mirror/ont-end-reason/index.html) ·
[Issues](https://github.com/Single-Molecule-Sequencing/ont-end-reason/issues) ·
[CHANGELOG](CHANGELOG.md) ·
[Design spec](docs/superpowers/specs/2026-05-12-ont-end-reason-design.md)

</div>

See [missing end-reason metadata](docs/missing-end-reason.md) for the UNKNOWN
status, nullable quality metrics, and baseline-store behavior.
