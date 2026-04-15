#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lcs_benchmark.py
================
LCS Algorithm Benchmark Suite — IEEE Manuscript Figure Generator

Produces two publication-quality figures (SVG + PNG + EPS) comparing four LCS
algorithms across five sequence pairs of increasing lengths.

Default plot style  : side-by-side sub-panels  (time | memory)  ← IEEE-recommended
Optional plot style : dual y-axis                                (--dual-axis flag)

Usage
-----
  # Fresh benchmark run (overwrites cache):
  python lcs_benchmark.py

  # Re-plot from previously saved results, skip benchmarking:
  python lcs_benchmark.py --reuse

  # Run only one chart:
  python lcs_benchmark.py --chart 1

  # Use the original dual-axis style instead of the default:
  python lcs_benchmark.py --dual-axis

  # Quick smoke-test; no .fa files required:
  python lcs_benchmark.py --dry-run
"""
from __future__ import annotations

import argparse
import json
import random
import time
import tracemalloc
import sys
import os
from pathlib import Path
from collections import deque
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

# VS Code's relative path problem, learn why: https://stackoverflow.com/questions/58709973/
sys.path[0] = str(Path(sys.path[0]).parent.parent)
from lcs_dp import (BUILD_LCS, GET_LCS)
from lcs_hs import get_lcs_hs
from lcs_lsbs import lcs_lsbs
from lis_funcs import (lis_dp_2d, lis_ps_2d)

# Embed fonts as live text in SVG (not outlines) → crisp zoom in Word / LaTeX
matplotlib.rcParams["svg.fonttype"] = "none"
# Embed TrueType fonts (Type 42) in EPS — required for clean Word import
matplotlib.rcParams["ps.fonttype"]  = 42


# =============================================================================
# §0  CONFIGURATION  ── edit this block to suit your environment
# =============================================================================

# VS Code
FA_DIR      = Path("benchmark/python/data")   # directory containing your .fa files
RESULTS_DIR = os.getcwd() / Path("benchmark/python/results")          # cached benchmark JSON
FIGURES_DIR = os.getcwd() / Path("benchmark/python/figures")          # output SVG / PNG figures

# PyCharm
# FA_DIR      = Path(f"{os.getcwd()}/data/sequences")   # directory containing your .fa files
# RESULTS_DIR = os.getcwd() / Path(f"{os.getcwd()}/results")          # cached benchmark JSON
# FIGURES_DIR = os.getcwd() / Path(f"{os.getcwd()}/figures")          # output SVG / PNG figures

# Exactly ten .fa filenames.
# Pairs: (FA_FILES[0], FA_FILES[1]), (FA_FILES[2], FA_FILES[3]), …
FA_FILES: list[str] = [
    "chr1_GL383518v1_alt.fa", "chr2_GL383521v1_alt.fa",
    "chr3_GL383526v1_alt.fa", "chr4_GL000257v2_alt.fa",
    "chr5_GL339449v2_alt.fa", "chr6_GL000250v2_alt.fa",
    "chr7_GL383534v2_alt.fa", "chr8_KI270810v1_alt.fa",
    "chr9_GL383539v1_alt.fa", "chr10_GL383545v1_alt.fa",
]

# Sequence length (characters) per pair
CHART1_LENGTHS: list[int] = [50, 100, 150, 200, 250]
CHART2_LENGTHS: list[int] = [50, 100, 150, 200, 250]   # mirrors Chart 1 — same .fa source data

# Characters used for the random half of Chart-2 strings
GENOME_ALPHABET: str = "ACGT"

# N-heavy alphabet for the dense random fill in Chart 2
# Weighted toward 'N' (unknown sequence indicator) to maximise match density
CHART2_DENSE_ALPHABET: str = "N" * 6 + "ACG"   # ≈ 67 % N, remainder A/C/G

RANDOM_SEED: int = 42

# ── Visual style ─────────────────────────────────────────────────────────────
#   False → IEEE-recommended default: two side-by-side sub-panels
#   True  → user's original vision:  one plot with dual y-axes
#   (override at runtime with --dual-axis)
DUAL_AXIS_MODE: bool = False

# Algorithm display names — order MUST match ALGORITHMS list in §1
ALGO_LABELS: list[str] = [
    "Std-DP",
    "Hunt-Szymanski",
    "LSBS-DP (ALGO 2)",
    "LSBS-PS (ALGO 3)",
]

_MARKERS:     list[str] = ["o", "s", "^", "D"]
_LINESTYLES:  list[str] = ["-", "--", "-.", ":"]


# =============================================================================
# §1  ALGORITHM PLACEHOLDERS  ── wire your implementations here
# =============================================================================
#
# Replace the `raise NotImplementedError` bodies with imports / calls to your
# actual modules.  The lambdas in ALGORITHMS are the sole invocation point:
# edit only those lambdas (and the LIS callback) when hooking things up.
#
# =============================================================================


def algo_standard_dp(s1: str, s2: str) -> int:
    """Algorithm 1 — Standard O(|A|·|B|) Dynamic Programming LCS."""
    # raise NotImplementedError("Supply algo_standard_dp")
    return GET_LCS(s1, s2, BUILD_LCS(s1, s2, {}))


def algo_hunt_szymanski(s1: str, s2: str) -> int:
    """Algorithm 2 — Hunt-Szymanski LCS."""
    # raise NotImplementedError("Supply algo_hunt_szymanski")
    return get_lcs_hs(s1, s2)


def algo_lsbs_a2(s1: str, s2: str) -> int:
    """Algorithm 3 — LSBS with Algorithm-2 variant as the LIS callback."""
    # raise NotImplementedError("Supply algo_lsbs_a2")
    return lcs_lsbs(s1, s2, lis_dp_2d)[1]


def algo_lsbs_a3(s1: str, s2: str) -> int:
    """Algorithm 4 — LSBS with Algorithm-3 variant as the LIS callback."""
    # raise NotImplementedError("Supply algo_lsbs_a3")
    return lcs_lsbs(s1, s2, lis_ps_2d)[1]


# ── Single invocation point ── edit lambdas / callbacks here only ─────────────
ALGORITHMS: list[tuple[str, Callable[[str, str], int]]] = [
    (ALGO_LABELS[0], lambda a, b: algo_standard_dp(a, b)),
    (ALGO_LABELS[1], lambda a, b: algo_hunt_szymanski(a, b)),
    (ALGO_LABELS[2], lambda a, b: algo_lsbs_a2(a, b)),
    (ALGO_LABELS[3], lambda a, b: algo_lsbs_a3(a, b)),
]


# =============================================================================
# §2  FASTA UTILITIES  (Python 3.10-style; two-pass, memory-efficient)
# =============================================================================

_HEADER_CHARS: frozenset[str] = frozenset((">", ";"))


def _is_fasta_header(line: str) -> bool:
    return bool(line) and line[0] in _HEADER_CHARS


def _count_fasta_length(path: Path) -> int:
    """
    Pass-1: count total sequence characters without storing anything.
    1 MB read buffer for throughput on large chromosomal files.
    """
    total = 0
    with path.open("r", buffering=1 << 20) as fh:
        for raw in fh:
            if not _is_fasta_header(raw):
                total += len(raw.rstrip())
    return total


def sample_fasta_evenly(path: Path, target: int) -> str:
    """
    Return exactly `target` characters sampled at arithmetically-spaced
    positions across the full FASTA sequence.

    WHY EVENLY SPACED?
    Chromosomal sequences have position-dependent nucleotide composition
    (GC-rich centromeres, AT-rich telomeres, N-padded gaps).  Sampling only
    the first N characters would bias all pairs toward the same structural
    region, distorting match-density and invalidating cross-algorithm
    comparisons.  Evenly-spaced sampling gives each pair a representative
    cross-section of the chromosome regardless of its length target.

    Complexity: O(2 · total_len) time,  O(target) peak memory (two-pass).
    """
    total = _count_fasta_length(path)
    if total == 0:
        raise ValueError(f"No sequence data found in {path}")

    if total <= target:
        # Sequence shorter than requested — return everything
        parts: list[str] = []
        with path.open("r", buffering=1 << 20) as fh:
            for raw in fh:
                if not _is_fasta_header(raw):
                    parts.append(raw.rstrip())
        return "".join(parts)

    # Sorted, de-duplicated global character indices to collect
    sample_positions: list[int] = sorted(
        dict.fromkeys(
            np.linspace(0, total - 1, target, dtype=np.intp).tolist()
        )
    )
    n_targets = len(sample_positions)

    result: list[str] = []
    pos    = 0
    ptr    = 0      # pointer into sample_positions

    with path.open("r", buffering=1 << 20) as fh:
        for raw in fh:
            if _is_fasta_header(raw):
                continue
            line     = raw.rstrip()
            line_end = pos + len(line)

            # Collect all sample positions that fall within this line
            while ptr < n_targets and sample_positions[ptr] < line_end:
                result.append(line[sample_positions[ptr] - pos])
                ptr += 1

            pos = line_end
            if ptr >= n_targets:
                break   # no need to continue reading

    return "".join(result)


def sample_fasta_tail_mixed(path: Path, target: int, rng: random.Random) -> str:
    """
    Return a `target`-length string for the Dense-CET / high-overlap scenario:

        [ random ACGT (target//2) ] ++ [ N-heavy tail (target - target//2) ]

    RATIONALE
    The terminal tail of GRCh38 chromosomes carries extended N-runs (unknown
    assembly bases).  When two N-heavy strings are compared, almost every
    character in A matches every 'N' in B, producing a Dense Compressed
    Equivalence Table.  Algorithms that enumerate all matches before invoking
    their LIS subroutine (Hunt-Szymanski, LSBS variants) degenerate toward
    O(|A|·|B|) in this regime, exposing the contrast with Std-DP.

    IMPLEMENTATION
    Rolling deque — O(tail_len) peak memory regardless of chromosome size.
    No full-file load required.
    """
    tail_len   = target // 2
    random_len = target - tail_len

    buf: deque[str] = deque(maxlen=tail_len)
    with path.open("r", buffering=1 << 20) as fh:
        for raw in fh:
            if not _is_fasta_header(raw):
                buf.extend(raw.rstrip())

    tail_str   = "".join(buf)
    random_str = "".join(rng.choice(GENOME_ALPHABET) for _ in range(random_len))

    # Random prefix + N-heavy suffix — overlap is densest at the aligned end
    return random_str + tail_str


# =============================================================================
# §3  SEQUENCE PREPARATION
# =============================================================================

def _fa_path(name: str) -> Path:
    p = FA_DIR / name
    if not p.exists():
        raise FileNotFoundError(
            f"FASTA file not found: {p}\n"
            f"  → Verify FA_DIR ({FA_DIR.resolve()}) and FA_FILES in §0."
        )
    return p


def prepare_chart1_pairs(rng: random.Random) -> list[tuple[str, str]]:
    """
    Five pairs; each pair uses evenly-sampled portions of increasing length.
    Files are consumed two at a time (pair_index determines the file pair).
    """
    pairs: list[tuple[str, str]] = []
    for i, length in enumerate(CHART1_LENGTHS):
        fa_a = _fa_path(FA_FILES[i * 2])
        fa_b = _fa_path(FA_FILES[i * 2 + 1])
        print(f"    Pair {i+1} | L={length:>6,} | {fa_a.name}  ×  {fa_b.name}")
        pairs.append((
            sample_fasta_evenly(fa_a, length),
            sample_fasta_evenly(fa_b, length),
        ))
    return pairs


def prepare_chart2_pairs(rng: random.Random) -> list[tuple[str, str]]:
    """
    Five pairs built from the same .fa files as Chart 1.

    Each string of target length L is constructed as:
        [ evenly-sampled FASTA chars  (L // 2) ]
     ++ [ N-dense random fill         (L - L // 2) ]

    The N-dense suffix (drawn from CHART2_DENSE_ALPHABET, ~67 % 'N') creates
    a high-overlap tail that degenerates match-enumeration algorithms by
    producing a Dense Compressed Equivalence Table while the first half retains
    realistic nucleotide composition for a fair cross-algorithm comparison.
    """
    pairs: list[tuple[str, str]] = []
    for i, length in enumerate(CHART2_LENGTHS):
        fa_a = _fa_path(FA_FILES[i * 2])
        fa_b = _fa_path(FA_FILES[i * 2 + 1])
        half = length // 2
        fill = length - half
        print(f"    Pair {i+1} | L={length:>6,} | {fa_a.name}  ×  {fa_b.name}  [50 % FASTA + 50 % N-dense]")
        fasta_a = sample_fasta_evenly(fa_a, half)
        fasta_b = sample_fasta_evenly(fa_b, half)
        fill_a  = "".join(rng.choice(CHART2_DENSE_ALPHABET) for _ in range(fill))
        fill_b  = "".join(rng.choice(CHART2_DENSE_ALPHABET) for _ in range(fill))
        pairs.append((fasta_a + fill_a, fasta_b + fill_b))
    return pairs


# ── Dry-run: synthetic pairs (no .fa files needed) ───────────────────────────

def _make_dry_run_pairs(lengths: list[int],
                         rng: random.Random,
                         n_heavy: bool = False) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for n in lengths:
        if n_heavy:
            half = n // 2
            s1 = "".join(rng.choice(GENOME_ALPHABET) for _ in range(half)) + "N" * (n - half)
            s2 = "".join(rng.choice(GENOME_ALPHABET) for _ in range(half)) + "N" * (n - half)
        else:
            s1 = "".join(rng.choice(GENOME_ALPHABET) for _ in range(n))
            s2 = "".join(rng.choice(GENOME_ALPHABET) for _ in range(n))
        pairs.append((s1, s2))
    return pairs


# =============================================================================
# §4  BENCHMARK RUNNER
# =============================================================================

@dataclass
class RunRecord:
    algo:        str
    pair_idx:    int
    length:      int
    elapsed_ms:  float   # wall-clock time in milliseconds
    peak_mem_kb: float   # peak tracemalloc delta in kilobytes
    lcs_len:     int     # sanity-check: LCS length returned by algorithm

    def to_dict(self) -> dict:
        return asdict(self)


def _bench_once(label: str,
                fn: Callable[[str, str], int],
                s1: str, s2: str,
                pair_idx: int,
                length: int) -> RunRecord:
    """Single timed + memory-traced invocation of one algorithm on one pair."""
    was_tracing = tracemalloc.is_tracing()
    if not was_tracing:
        tracemalloc.start()

    t0 = time.perf_counter()
    result = fn(s1, s2)
    elapsed_ms = (time.perf_counter() - t0) * 1_000.0

    _, peak_bytes = tracemalloc.get_traced_memory()
    if not was_tracing:
        tracemalloc.stop()

    return RunRecord(
        algo=label,
        pair_idx=pair_idx,
        length=length,
        elapsed_ms=elapsed_ms,
        peak_mem_kb=int(peak_bytes) / 1_024.0,
        lcs_len=result  #len(result)  # int(result),
    )


def run_benchmark(chart: int,
                  pairs: list[tuple[str, str]],
                  lengths: list[int]) -> list[RunRecord]:
    """
    Iterate over all pairs × all algorithms; collect RunRecords.
    Algorithms that raise NotImplementedError are silently skipped.
    """
    records: list[RunRecord] = []
    for pair_idx, (s1, s2) in enumerate(pairs):
        length = lengths[pair_idx]
        print(f"\n  Pair {pair_idx+1}/{len(pairs)}  L={length:,}", flush=True)

        for label, fn in ALGORITHMS:
            try:
                rec = _bench_once(label, fn, s1, s2, pair_idx, length)
                records.append(rec)
                print(
                    f"    {label:<22}  "
                    f"{rec.elapsed_ms:>10.2f} ms   "
                    f"{rec.peak_mem_kb:>10.1f} KB",
                    flush=True,
                )
            except NotImplementedError:
                print(f"    {label:<22}  [NOT IMPLEMENTED — skipped]")
    return records


# =============================================================================
# §5  CACHING
# =============================================================================

_CACHE_FILES: dict[int, str] = {
    1: "chart1_benchmark.json",
    2: "chart2_benchmark.json",
}


def save_cache(chart: int, records: list[RunRecord]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / _CACHE_FILES[chart]
    with out.open("w") as fh:
        json.dump([r.to_dict() for r in records], fh, indent=2)
    print(f"\n  Benchmark data cached → {out}")


def load_cache(chart: int) -> list[RunRecord] | None:
    path = RESULTS_DIR / _CACHE_FILES[chart]
    if not path.exists():
        return None
    with path.open() as fh:
        data = json.load(fh)
    return [RunRecord(**d) for d in data]


# =============================================================================
# §6  PLOTTING
# =============================================================================
#
# DEFAULT  (DUAL_AXIS_MODE = False):
#   Two side-by-side sub-panels per figure.
#   Panel (a) — elapsed time  (ms)
#   Panel (b) — peak memory   (KB)
#   One line per algorithm; 4-entry legend; clean IEEE ticks + despined axes.
#
# OPTIONAL (DUAL_AXIS_MODE = True  /  --dual-axis flag):
#   Single axes, primary y = elapsed time, secondary y = peak memory.
#   8 lines total (4 algos × 2 metrics).
#
# WHY DUAL-AXIS IS NOT THE DEFAULT — honest assessment
# ─────────────────────────────────────────────────────
# 1.  8 simultaneous lines at IEEE two-column width (≈ 3.5 in) is genuinely
#     unreadable — especially when printed in greyscale.
# 2.  The two independent y-scales create a spurious visual correlation:
#     a time curve "crossing" a memory curve implies no real relationship.
# 3.  IEEE editorial guidance (§9.5) explicitly discourages twin-axis figures
#     when the axes do not share a physical relationship.
# 4.  A legend with 8 entries at 6.5 pt competes with the data for ink budget.
# 5.  Algorithm reviewers routinely flag this layout; it delays acceptance.
#
# The dual-axis version is retained as a fully functional configurable option
# for posters, supplementary material, or repository README figures where
# compactness outweighs legibility constraints.
# =============================================================================

_IEEE_RC: dict = {
    "font.family":    "serif",
    "font.serif":     ["Lato", "Times New Roman", "CMU Serif", "DejaVu Serif", "Georgia"],
    "font.size":       9,
    "axes.titlesize":  8,
    "axes.labelsize":  7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "figure.dpi":      150,
    "svg.fonttype":   "none",
}

_PALETTE = sns.color_palette("colorblind", n_colors=4)


def _records_to_df(records: list[RunRecord]) -> pd.DataFrame:
    return pd.DataFrame([r.to_dict() for r in records])


def _xtick_labels(lengths: list[int]) -> list[str]:
    """
    Magnitude-aware tick formatter.
    Values ≥ 1 000 → '1.5k' style; values < 1 000 → plain integer.
    Always suffixed with the pair label '(P{n})' on a second line.
    Adapts automatically as CHART1_LENGTHS / CHART2_LENGTHS are changed.
    """
    def _fmt(n: int) -> str:
        return f"{n / 1_000:g}k" if n >= 1_000 else str(n)
    return [f"{_fmt(n)}\n(P{i+1})" for i, n in enumerate(lengths)]


def _fmt_axis(ax: plt.Axes, xticks: list[int], xlabels: list[str],
               ylabel: str) -> None:
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_xlabel("Sequence Length", labelpad=4)
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f"{v:,.0f}")
    )
    sns.despine(ax=ax)


# ── §6-A  Default: side-by-side sub-panels ────────────────────────────────────

def _add_stats_table(
    ax: plt.Axes,
    df: pd.DataFrame,
    metric: str,
    fmt: str = "{:.2f}",
    bbox: list | None = None,
) -> None:
    """
    Overlay a compact MIN / MAX / AVG summary table in the top-left area of `ax`.

    Intended for the free "verandah" space produced by reverse-Ogive curves
    (low values at short lengths, steep rise toward the right) — the top-left
    quadrant is naturally unoccupied.

    Parameters
    ----------
    ax     : target axes
    df     : full records DataFrame (must contain 'algo' and `metric` columns)
    metric : column name to summarise ('elapsed_ms' or 'peak_mem_kb')
    fmt    : Python format string applied to each cell value
    bbox   : [left, bottom, width, height] in axes-fraction coordinates;
             defaults to a compact top-left position
    """
    if bbox is None:
        bbox = [0.28, 0.63, 0.5, 0.35]

    rows, cells = [], []
    for label, _ in ALGORITHMS:
        vals = df[df["algo"] == label][metric]
        if vals.empty:
            continue
        rows.append(label)
        cells.append([
            fmt.format(vals.min()),
            fmt.format(vals.max()),
            fmt.format(vals.mean()),
        ])

    if not rows:
        return

    tbl = ax.table(
        cellText=cells,
        rowLabels=rows,
        colLabels=["MIN", "MAX", "AVG"],
        loc="upper left",
        cellLoc="right",
        bbox=bbox,
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(5)

    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.3)
        cell.set_edgecolor("#cccccc")
        if r == 0:                              # column-header row
            cell.set_text_props(fontweight="bold", color="#333333")
            cell.set_facecolor("#e8e8e8")
        else:
            cell.set_facecolor("white")
        cell.set_alpha(0.88)


def _plot_side_by_side(records: list[RunRecord], lengths: list[int],
                        chart_title: str, subtitle: str) -> plt.Figure:
    """
    IEEE-recommended layout.  Each metric gets its own panel — no scale
    ambiguity, no cognitive overload.  Shared colour + marker coding makes
    cross-panel algorithm tracking trivial.  Figure width = 7.16 in (IEEE
    double-column).
    """
    with plt.rc_context(_IEEE_RC):
        sns.set_theme(
            style="ticks", context="paper", font_scale=1.0,
            rc={"axes.spines.right": False, "axes.spines.top": False},
        )

        fig, (ax_t, ax_m) = plt.subplots(
            1, 2,
            figsize=(7.16, 3.2),
            constrained_layout=True,
        )

        df      = _records_to_df(records)
        xticks  = list(range(len(lengths)))
        xlabels = _xtick_labels(lengths)

        for i, (label, _) in enumerate(ALGORITHMS):
            sub = df[df["algo"] == label].sort_values("pair_idx")
            if sub.empty:
                continue
            kw = dict(
                color=_PALETTE[i],
                marker=_MARKERS[i],
                linestyle=_LINESTYLES[i],
                linewidth=1.6,
                markersize=5,
                label=label,
            )
            ax_t.plot(sub["pair_idx"].tolist(), sub["elapsed_ms"].tolist(),  **kw)
            ax_m.plot(sub["pair_idx"].tolist(), sub["peak_mem_kb"].tolist(), **kw)

        _fmt_axis(ax_t, xticks, xlabels, "Elapsed Time (ms)\n")
        _fmt_axis(ax_m, xticks, xlabels, "Peak Memory (KB)\n")
        ax_t.set_title("(a)  Runtime",       pad=6)
        ax_m.set_title("(b)  Memory Usage",  pad=6)

        _add_stats_table(ax_t, df, "elapsed_ms",  fmt="{:.3f}")
        _add_stats_table(ax_m, df, "peak_mem_kb", fmt="{:.1f}")

        handles, labels = ax_t.get_legend_handles_labels()
        fig.legend(
            handles, labels,
            loc="lower center",
            ncol=len(ALGORITHMS),
            frameon=False,
            bbox_to_anchor=(0.5, -0.08),
        )
        fig.suptitle(chart_title, fontsize=9, fontweight="bold", y=1.03)
        fig.text(
            0.5, 1.05, subtitle,
            ha="center", va="bottom",
            fontsize=6.5, style="italic", color="#555555",
            transform=fig.transFigure,
        )
        return fig


# ── §6-B  Optional: dual y-axis ───────────────────────────────────────────────

def _plot_dual_axis(records: list[RunRecord], lengths: list[int],
                     chart_title: str, subtitle: str) -> plt.Figure:
    """
    User's original vision: single axes, twin y-scales.
    Solid lines → elapsed time (left axis).
    Dashed lines → peak memory (right axis, grey ticks).
    """
    with plt.rc_context(_IEEE_RC):
        sns.set_theme(style="ticks", context="paper", font_scale=1.0)

        fig, ax1 = plt.subplots(figsize=(7.16, 3.6), constrained_layout=True)
        ax2 = ax1.twinx()

        df      = _records_to_df(records)
        xticks  = list(range(len(lengths)))
        xlabels = _xtick_labels(lengths)

        all_handles: list       = []
        all_labels:  list[str]  = []

        for i, (label, _) in enumerate(ALGORITHMS):
            sub = df[df["algo"] == label].sort_values("pair_idx")
            if sub.empty:
                continue
            xs = sub["pair_idx"].tolist()
            l_t, = ax1.plot(
                xs, sub["elapsed_ms"].tolist(),
                color=_PALETTE[i], marker=_MARKERS[i],
                linestyle="-", linewidth=1.6, markersize=5,
                label=f"{label} [time]",
            )
            l_m, = ax2.plot(
                xs, sub["peak_mem_kb"].tolist(),
                color=_PALETTE[i], marker=_MARKERS[i],
                linestyle="--", linewidth=1.0, markersize=4, alpha=0.65,
                label=f"{label} [mem]",
            )
            all_handles += [l_t, l_m]
            all_labels  += [f"{label} – time", f"{label} – mem"]

        ax1.set_xticks(xticks)
        ax1.set_xticklabels(xlabels)
        ax1.set_xlabel("Sequence Length", labelpad=4)
        ax1.set_ylabel("Elapsed Time (ms)")
        ax1.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))

        ax2.set_ylabel("Peak Memory (KB)", color="#777777")
        ax2.tick_params(axis="y", labelcolor="#777777")
        ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        sns.despine(ax=ax1, right=False)

        # Two stacked tables in the top-left verandah; time above, memory below
        _add_stats_table(ax1, df, "elapsed_ms",  fmt="{:.3f}", bbox=[0.01, 0.65, 0.37, 0.29])
        _add_stats_table(ax1, df, "peak_mem_kb", fmt="{:.1f}", bbox=[0.01, 0.34, 0.37, 0.29])

        fig.legend(
            all_handles, all_labels,
            loc="lower center",
            ncol=4, frameon=False,
            bbox_to_anchor=(0.5, -0.12),
            fontsize=6.5,
        )
        fig.suptitle(chart_title, fontsize=9, fontweight="bold")
        fig.text(
            0.5, 0.97, subtitle,
            ha="center", va="top",
            fontsize=6.5, style="italic", color="#555555",
            transform=fig.transFigure,
        )
        return fig


# ── Dispatcher ────────────────────────────────────────────────────────────────

def build_figure(records: list[RunRecord], lengths: list[int],
                  chart_title: str, subtitle: str) -> plt.Figure:
    if DUAL_AXIS_MODE:
        return _plot_dual_axis(records, lengths, chart_title, subtitle)
    return _plot_side_by_side(records, lengths, chart_title, subtitle)


# =============================================================================
# §7  EXPORT  (SVG — Word/LaTeX lossless;  PNG 300 dpi — raster fallback;
#              EPS — PostScript vector, preferred for older Word / journal CMS)
# =============================================================================

def export_figure(fig: plt.Figure, stem: str) -> None:
    """
    SVG output notes
    ----------------
    • 'svg.fonttype = none' (set at module level) keeps font glyphs as <text>
      elements rather than converting them to paths.  Word 2016+ renders this
      correctly and the text remains selectable/searchable in PDF viewers.
    • Insert into Word via  Insert → Pictures → (select .svg).
      The figure scales losslessly at any zoom level.
    • For LaTeX: \\includegraphics{figures/chartN.svg}  with the svg package,
      or convert via Inkscape CLI: inkscape chart1.svg --export-pdf=chart1.pdf

    EPS output notes
    ----------------
    • 'ps.fonttype = 42' (set at module level) embeds TrueType fonts as Type 42
      PostScript fonts — the format Word's EPS importer expects.
    • No additional Python packages are required; Matplotlib ships its own
      PostScript/EPS backend (matplotlib.backends.backend_ps).
    • RECOMMENDED system tool: Ghostscript (gs / gswin64c).
      Without it the output is still valid EPS, but Ghostscript post-processes
      font embedding for maximum compatibility with Word, Acrobat, and journal
      submission portals.  Install via your OS package manager:
          Ubuntu/Debian : sudo apt install ghostscript
          macOS (Brew)  : brew install ghostscript
          Windows       : https://www.ghostscript.com/releases/
    • Transparency caveat: EPS does not support alpha compositing natively.
      In --dual-axis mode the dashed memory lines (alpha=0.65) are rasterised
      at 300 dpi inside the EPS envelope — quality remains high but those
      elements are no longer pure vectors.
    • Insert into Word via  Insert → Object → Create from File → (select .eps).
      Word 2019+ and Microsoft 365 render EPS natively at full resolution.
    """
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    svg_path = FIGURES_DIR / f"{stem}.svg"
    png_path = FIGURES_DIR / f"{stem}.png"
    eps_path = FIGURES_DIR / f"{stem}.eps"

    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
    fig.savefig(eps_path, format="eps", bbox_inches="tight")

    print(f"  → {svg_path}   (SVG — lossless, Word 2016+ / LaTeX ready)")
    print(f"  → {png_path}   (PNG 300 dpi — raster fallback)")
    print(f"  → {eps_path}   (EPS — PostScript vector, Word / journal CMS)")


# =============================================================================
# §8  CHART METADATA
# =============================================================================

_CHART_META: dict[int, tuple[str, str]] = {
    1: (
        "LCS Benchmark — Uniform Sequence Distribution",
        "Five chromosome pairs at increasing lengths (50–250 chars); evenly sampled from GRCh38",
    ),
    2: (
        "LCS Benchmark — High-Overlap (N-Dense) Sequences",
        (
            "50% evenly-sampled GRCh38 sequence + 50% N-dense random fill (≈67 % N, remainder A/C/G); "
            "exposes Dense-CET degeneration"
        ),
    ),
}


# =============================================================================
# §9  MAIN / CLI
# =============================================================================

def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="LCS Algorithm Benchmark — IEEE Manuscript Figure Generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "--reuse", action="store_true",
        help="Skip benchmarking; re-plot from cached JSON in results/",
    )
    ap.add_argument(
        "--chart", type=int, choices=[1, 2], default=None,
        help="Process only one chart (default: both)",
    )
    ap.add_argument(
        "--dual-axis", action="store_true",
        help="Use dual y-axis style (not recommended for IEEE manuscripts)",
    )
    ap.add_argument(
        "--dry-run", action="store_true",
        help="Use tiny synthetic sequences; tests pipeline without .fa files",
    )
    return ap.parse_args()


def main() -> None:
    args = _parse_args()

    global DUAL_AXIS_MODE
    if args.dual_axis:
        DUAL_AXIS_MODE = True
        print(
            "\n[NOTE] Dual-axis mode enabled.  "
            "See §6 in the script and BENCHMARK_ADDENDUM.md §3 for a detailed "
            "explanation of why the side-by-side layout is preferred for "
            "IEEE submissions.\n"
        )

    charts_to_run = [1, 2] if args.chart is None else [args.chart]
    rng           = random.Random(RANDOM_SEED)

    for chart in charts_to_run:
        title, subtitle = _CHART_META[chart]
        lengths         = CHART1_LENGTHS if chart == 1 else CHART2_LENGTHS

        print(f"\n{'━' * 62}")
        print(f"  CHART {chart} — {title}")
        print(f"{'━' * 62}")

        # Decide: fresh benchmark or load from cache
        run_fresh = not args.reuse
        records: list[RunRecord] = []

        if args.reuse:
            cached = load_cache(chart)
            if cached is None:
                print(f"  [!] No cache found for chart {chart} — running fresh.")
                run_fresh = True
            else:
                records = cached
                print(f"  Loaded {len(records)} records from cache.")

        if run_fresh:
            if args.dry_run:
                dry_lengths = [100, 200, 300, 400, 500]
                print(f"  [dry-run] Synthetic sequences, lengths: {dry_lengths}")
                pairs   = _make_dry_run_pairs(dry_lengths, rng, n_heavy=(chart == 2))
                lengths = dry_lengths
            else:
                print("  Preparing sequences from .fa files …")
                pairs = prepare_chart1_pairs(rng) if chart == 1 \
                        else prepare_chart2_pairs(rng)

            print("\n  Benchmarking …")
            records = run_benchmark(chart, pairs, lengths)

            if not args.dry_run:
                save_cache(chart, records)

        if not records:
            print(
                "\n  [!] No benchmark records (all algorithms unimplemented?).\n"
                "      Wire your implementations in §1 and re-run."
            )
            continue

        print("\n  Building figure …")
        fig = build_figure(records, lengths, title, subtitle)

        stem = f"chart{chart}{'_dualaxis' if DUAL_AXIS_MODE else ''}"
        export_figure(fig, stem)
        plt.close(fig)

    print(f"\n✓  Done.  Figures saved in: {FIGURES_DIR.resolve()}")


if __name__ == "__main__":
    main()