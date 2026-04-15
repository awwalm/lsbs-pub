# Benchmarks of LCS & LSBS

## § Purpose

* §0 All configurations — FA_DIR, file list, lengths, DUAL_AXIS_MODE

* §1 The four algorithm placeholders + _lis_placeholder + ALGORITHMS (single invocation point)

* §2 FASTA utilities — sample_fasta_evenly and sample_fasta_tail_mixed

* §3 Sequence preparation for Chart 1 (even) and Chart 2 (N-heavy + random)

* §4 Benchmark runner using tracemalloc + time.perf_counter()

* §5 JSON cache — save_cache / load_cache

* §6 Two plot styles (default side-by-side; `--dual-axis`; optional)

* §7 SVG + PNG export

* §8–9 Chart metadata + CLI (`--reuse, --chart, --dual-axis, --dry-run`)

Dataset directly downloadable from https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/ 
(match the names from any resulting error messages), and position them in the required folders