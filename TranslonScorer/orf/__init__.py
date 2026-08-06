"""ORF table handling: import, mapping, schema and benchmark panels.

- orfs_import:    BED12 → transcript-mapped canonical ORF tables
- map_orfs:       ORFs → per-ORF feature chains and slices
- score_schema:   score-table schema helpers (stable ORF keys, column sets)
- panel_manifest: frozen benchmark panel + sha256 provenance

Scoring itself lives on the spine (events → scoring/ → report), not here.
"""
