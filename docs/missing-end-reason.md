# Missing end-reason metadata

Distribution analysis reports `quality_status: UNKNOWN` when any read has a
missing, blank, explicitly unknown, or unrecognized end reason. It cannot infer
run quality from those reads. The three derived metrics `signal_positive_pct`,
`unblock_mux_pct`, and `data_service_pct` are `null` in JSON and unavailable in
CLI text. Incomplete metadata is never written to the QC baseline store.

`recognized_reads` and `unrecognized_reads` report coverage. Counts and per-label
fractions remain available for inspection, with the total read count as their
denominator. They describe recorded labels, including unknown values; they do not
estimate the missing end reasons. A mixture of known and unknown values remains
UNKNOWN, even if every known value is signal-positive.

Canonical full names and short codes are accepted regardless of case or surrounding
whitespace. The canonical `partial`/`PART` category is a recorded truncation reason,
not missing metadata. Complete metadata continues to use the existing OK, CHECK,
and FAIL thresholds. A fully classified run with no signal-positive reads still
has a measured signal-positive percentage of zero and a FAIL status.

A readable input with unavailable metadata exits successfully because the counts
were computed; consumers must inspect `quality_status`. Missing files and inputs
with no reads remain errors. Summary tables preserve null metrics, and reports
show UNKNOWN with an explanation instead of a failed-quality interpretation.
