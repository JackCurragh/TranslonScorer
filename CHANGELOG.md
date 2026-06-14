# Changelog

All notable changes to TranslonScorer are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Scalable matrix query engine (`pipeline/matrix_rollup.py`): offset-independent
  alignment index, map-reduce frame rollup, multimapping strategies
  (`unique` / guilty-by-association), per-sample offset calibration, oxbow
  reader, and region-scoped positional profiles for clustering. See
  `docs/matrix_query_engine.md`.
- Per-read-length frame-dominance emission in `pipeline/matrix_qc.py`
  (`frame_dominance_matrix`).
- Release tooling: bump2version config, CI + PyPI Trusted-Publishing workflows,
  this changelog, and `CITATION.cff`.

## [0.1.1]

- Initial alpha: CLI for scoring translational events from Ribo-seq data.
