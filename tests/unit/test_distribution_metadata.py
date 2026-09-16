"""Missing or incomplete end-reason metadata is not measured poor quality."""

from __future__ import annotations

import importlib
import json

import pytest
from click.testing import CliRunner

from ont_end_reason.analyze.distribution import distribution, maybe_store_baseline
from ont_end_reason.analyze.tables import generate_tables
from ont_end_reason.cli import main
from ont_end_reason.errors import AnalysisError
from ont_end_reason.io.manifest import ReadRecord


def _records(reasons):
    return [ReadRecord(read_id=str(i), end_reason=value) for i, value in enumerate(reasons)]


@pytest.mark.parametrize(
    "reasons",
    [
        ["unknown"],
        ["UNK"],
        [""],
        [" "],
        ["new_reason"],
        ["signal_positive", "unknown"],
        ["signal_positive", "new_reason"],
    ],
)
def test_unavailable_or_partial_metadata_is_unknown(reasons):
    result = distribution(_records(reasons))
    assert result.quality_status == "UNKNOWN"
    assert result.signal_positive_pct is None
    assert result.unblock_mux_pct is None
    assert result.data_service_pct is None
    assert sum(result.counts.values()) == len(reasons)
    assert sum(result.percentages.values()) == pytest.approx(1)
    assert result.recognized_reads == reasons.count("signal_positive")
    assert result.unrecognized_reads == len(reasons) - result.recognized_reads
    assert "unavailable" in result.interpretation.lower()
    assert result.to_dict()["signal_positive_pct"] is None


@pytest.mark.parametrize(
    "reasons,status,percent",
    [
        (["signal_positive"] * 4 + ["signal_negative"], "OK", 80),
        (["signal_positive"] * 3 + ["signal_negative"] * 2, "CHECK", 60),
        (["signal_negative"] * 5, "FAIL", 0),
        (["partial"], "FAIL", 0),
        (["SP"], "OK", 100),
    ],
)
def test_complete_metadata_preserves_quality_gates(reasons, status, percent):
    result = distribution(_records(reasons))
    assert result.quality_status == status
    assert result.signal_positive_pct == percent
    assert result.recognized_reads == len(reasons)
    assert result.unrecognized_reads == 0


def test_empty_input_remains_error():
    with pytest.raises(AnalysisError, match="No reads"):
        distribution([])


@pytest.fixture(params=[None, "", "unknown", "signal_positive\nr2\tunknown"])
def summary(tmp_path, request):
    path = tmp_path / "sequencing_summary.txt"
    if request.param is None:
        text = "read_id\tsequence_length_template\nr1\t100\n"
    else:
        text = "read_id\tend_reason\nr1\t" + request.param + "\n"
    path.write_text(text, encoding="utf-8")
    return path


def test_summary_missing_blank_unknown_or_mixed(summary):
    result = distribution(summary)
    assert result.quality_status == "UNKNOWN"
    assert result.to_dict()["signal_positive_pct"] is None


def test_baseline_unknown_skips_before_registry_access(summary, monkeypatch):
    module = importlib.import_module("ont_end_reason.analyze.distribution")

    def forbidden(*args, **kwargs):
        pytest.fail("Incomplete metadata reached the baseline store path")

    monkeypatch.setattr(module, "_load_registry_entries", forbidden)
    assert maybe_store_baseline(distribution(summary), summary) is None


def test_summary_table_keeps_unavailable_values(summary):
    table = generate_tables(summary, name="summary")
    assert table.rows[0]["quality_status"] == "UNKNOWN"
    assert table.rows[0]["signal_positive_pct"] is None
    assert "UNKNOWN" in table.render()


@pytest.mark.parametrize("command", ["distribution", "stats"])
def test_cli_unknown_text_and_json(summary, tmp_path, command):
    out = tmp_path / "result.json"
    args = (
        ["analyze", "distribution", str(summary), "--no-baseline-store"]
        if command == "distribution"
        else ["stats", "--summary", str(summary)]
    )
    result = CliRunner().invoke(main, [*args, "--json", str(out)])
    assert result.exit_code == 0, result.output
    assert "UNKNOWN" in result.output
    assert "unavailable" in result.output.lower()
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["quality_status"] == "UNKNOWN"
    assert payload["signal_positive_pct"] is None


def test_cli_missing_file_remains_error(tmp_path):
    result = CliRunner().invoke(
        main, ["analyze", "distribution", str(tmp_path / "missing.txt"), "--no-baseline-store"]
    )
    assert result.exit_code != 0
    assert "Traceback" not in result.output


def test_case_whitespace_short_codes_normalize_without_unknown():
    result = distribution(_records([" SP ", "SIGNAL_POSITIVE", "  uMc ", " Partial "]))
    assert result.counts == {"signal_positive": 2, "unblock_mux_change": 1, "partial": 1}
    assert result.quality_status == "CHECK"
    assert result.recognized_reads == 4


def test_unknown_aliases_normalize_and_preserve_unrecognized_labels():
    result = distribution(_records([" UnK ", " UNKNOWN ", "new_reason"]))
    assert result.counts == {"unknown": 2, "new_reason": 1}
    assert result.quality_status == "UNKNOWN"


def test_blank_summary_values_are_unknown_not_nan(tmp_path):
    path = tmp_path / "sequencing_summary.txt"
    path.write_text("read_id\tend_reason\nr1\t\n", encoding="utf-8")
    assert distribution(path).counts == {"unknown": 1}


def test_distribution_renderers_preserve_unknown(summary):
    import matplotlib.pyplot as plt

    from ont_end_reason.report.html import _distribution_section
    from ont_end_reason.viz.static import plot_distribution

    result = distribution(summary)
    fig = plot_distribution(result)
    assert "UNKNOWN" in fig.axes[0].get_title()
    plt.close(fig)
    _, html = _distribution_section(str(summary), quick=False)
    assert "status-UNKNOWN" in html
    assert "quality is unavailable" in html
    assert "Failed:" not in html
