# Changelog

All notable changes to MEIsensor are documented here.

## 2.1.0 — annotation revision (2026-11-27)

优化了 MEI 的注释方法（Improved the MEI annotation method）：TSD、poly(A/T) 与 EN motif
改为由 MEIsensor 内部按 TSD 边界直接计算。

### Changed
- TSD and poly(A/T) are computed natively by `src/annotate.py`; annotation is on by
  default when `--reference` is supplied and can be disabled with `--annotate false`.
- EN motif is anchored at the physical TSD boundary instead of being searched for in a
  ±25 bp window: a fixed 12-bp window (positions −6…−1 and +1…+6) is read on the nicked
  strand; the 6-bp core (−4…−1, +1, +2) is reported in `EN_MOTIF`, the nicked strand in
  `EN_STRAND`, and `EN_OFFSET` is 0.
- The `--en-motif-window` option and the candidate-ranking motif search were removed.
- TSD: longest exact (0-mismatch) match between the inserted sequence and either flank,
  ≥4 bp, capped at 60 bp, never summed.
- poly(A/T): terminal A-rich (3′) or T-rich (5′) tract after TSD removal, ≥10 bp,
  ≥80% purity in every 10-bp window, no gap bridging.
- `--allow-overwrite` documented for reruns into existing output paths.

### License
- Added the MIT License (`LICENSE`), copyright (c) 2026 Yudong Wang, Kai Ye, Xiaofei Yang, and Xi'an Jiaotong University.
- `README.md` now declares the MIT License and links to `LICENSE`.


