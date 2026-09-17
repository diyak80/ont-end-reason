"""Minimal exact 1 bp length-count invariants for plotting.

Reimplemented locally rather than importing ``lib.readdist`` from the
ont-ecosystem repository: ont-end-reason ships and installs as its own
standalone package (see ``pyproject.toml``), and making every install depend
on a private sibling git checkout for one plotting helper would be more
fragile than the ~20 lines this needs. The two invariants enforced here are
the same ones ``lib.readdist.core.LengthCounts.validate`` asserts in that
repository (see ``ont-ecosystem/lib/readdist/CONTRACT.md`` section 3):

* Conservation — the plotted counts sum to exactly the number of input
  reads. Nothing is dropped, capped or subsampled.
* Resolution — the number of distinct plotted x values equals the number of
  distinct integer lengths in the input. Nothing is rebinned.

Bin width is implicitly 1 (the support grid steps by 1), so there is no
separate bin-width invariant to check.
"""

from __future__ import annotations

import numpy as np


class ExactCountInvariantError(RuntimeError):
    """A conservation or resolution invariant failed. The plot would lie."""


def exact_length_counts(lengths) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(lengths, counts)`` at exact 1 bp resolution, no cap, no bin.

    ``lengths`` is the full per-read length population for one class (no
    subsampling upstream). The returned arrays span every integer from the
    observed minimum to the observed maximum inclusive, with explicit zeros
    where no read was observed, so a line or step plot drawn from them shows
    every base pair rather than a merged bin.

    Raises :class:`ExactCountInvariantError` if the counting failed to
    conserve every read or collapsed two distinct lengths together, which
    would otherwise be a silent data-reduction bug.
    """
    arr = np.asarray(lengths, dtype=np.int64)
    n = int(arr.size)
    if n == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)

    lo, hi = int(arr.min()), int(arr.max())
    counts = np.bincount(arr - lo, minlength=hi - lo + 1)

    plotted_total = int(counts.sum())
    if plotted_total != n:
        raise ExactCountInvariantError(
            f"conservation invariant failed: plotted counts sum to {plotted_total} "
            f"but the input has {n} reads (difference {plotted_total - n})"
        )

    n_distinct_source = int(np.unique(arr).size)
    n_distinct_plotted = int(np.count_nonzero(counts))
    if n_distinct_plotted != n_distinct_source:
        raise ExactCountInvariantError(
            f"resolution invariant failed: {n_distinct_plotted} distinct plotted "
            f"lengths but the input has {n_distinct_source} distinct lengths; "
            "two or more lengths were merged"
        )

    grid = np.arange(lo, hi + 1, dtype=np.int64)
    return grid, counts
