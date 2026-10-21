# Changelog

All notable changes to `versor-java` are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-10-21

First stable release.

### Added

- Arbitrary `Cl(p,q,r)` with up to 14 generators; precomputed sign tables for
  `n <= 10` and on-demand signs above that.
- Dense `Multivector`, sparse `SparseMultivector`, `MultivectorBatch`, and
  `MultivectorMatrix`.
- Geometric, outer, inner, regressive, scalar and sandwich products; reverse,
  conjugate, involution, dual/undual, norm, versor inverse and checked general
  inverse.
- `E3`, `P3` (PGA3), `Conformal` (Cl(4,1)), `Motion`, `Field`, `Calculus`, and
  `Manifold` layers.
- Optional Blaze3D/JOML adapter, built separately.
- Self-contained test harness: 454 checks, no external test dependency.
- Benchmarks, fuzz checks, and a reproducible performance report.

### Changed

- `Signature.metric` is now private. Read entries with `metric(int)` or get a
  defensive copy with `metricCopy()`.
