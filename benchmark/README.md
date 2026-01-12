# Benchmark

This folder contains scripts and results for benchmarking the performance of
**Metagenomic-DeepFRI** against **EggNOG-mapper** (v2.1.13).

## Contents

- `time_benchmark_cpu.py`: Python script that:
    1. Generates random protein sequences (10, 100, 1000, 10000) sampled from
       SwissProt.
    2. Introduces random mutations (5-80% of sequence length) to simulate
       random metagenomic proteins.
    3. Runs EggNOG-mapper (fast mode) and mDeepFRI on these sequences.
    4. Measures execution time and logs results.

- `benchmaprk_results.ipynb`: Jupyter Notebook for analyzing the benchmark
  logs and plotting the results. It separates mDeepFRI runtime into "Search"
  and "Inference" phases.

- `benchmark_results.tsv`: Tab-separated values file containing the raw
  execution times.

- `benchmark_time_cpu_comparison.png`: A plot comparing the execution times
  (log-log scale).

## Results

![Benchmark Result](benchmark_time_cpu_comparison.png)

### Comparison of Execution Time (CPU)

| Sequences | EggNOG (fast) [s] | mDeepFRI [s] | Speedup |
|-----------|-------------------|--------------|---------|
| 10        | 507.66            | 22.55        | 22.5x   |
| 100       | 591.51            | 68.68        | 8.6x    |
| 1000      | 1466.45           | 595.07       | 2.5x    |
| 10000     | 7504.77           | 6181.24      | 1.2x    |

mDeepFRI demonstrates superior performance compared to EggNOG-mapper (fast
mode) on CPU, particularly for smaller query sets. The initial overhead of
mDeepFRI is significantly lower. As the number of sequences increases, the
alignment/search phase dominates the runtime, bringing the total time closer to
that of EggNOG-mapper, yet mDeepFRI remains faster up to the tested 10,000
sequences. The stacked area in the plot for mDeepFRI distinguishes between the
database search time (dark blue) and the GCN inference time (light blue).
