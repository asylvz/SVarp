# Changelog

## 1.3.0 (September 2026)

Signals and clustering
- Breakpoint signals from clipped single alignments (`--min-clip`, default 500 bp).
- `generate_sv_node` rewritten as a path walk with strand-consistent positions; node names containing `:` handled.
- Neighbouring clusters merged by link orientation and cluster extent.
- Genome-wide depth for clusters on alt nodes; depth filters counted once per cluster.
- Reads without a GAF record are counted and, with `--write-unmapped`, listed; the unmapped-read assembly is gone.

Assembly
- Clusters assembled in parallel (`--asm-jobs`), each in its own work directory, with per-step time limits.
- Clusters of untagged reads are not assembled unless `--keep-untagged` is given.

Svtig output (Step 5)
- Svtigs are remapped with GraphAligner and kept when anchored in the graph (`--min-graph-cov`, records with identity >= `--min-identity`) and at least `--min-svtig-len` long.
- Noisy contig ends, 1 kb windows aligning below `--trim-identity`, are trimmed before these gates (`trim=` in the header).
- `graph_explained=yes|no` records whether a graph path reproduces the svtig within 50 bp (`max_gap`, `max_indel`); it does not filter.
- Svtigs whose path is the plain reference are dropped (`--keep-reference-svtigs` keeps them).
- Headers carry `path`, `graph_cov`, `graph_identity`, `max_gap`, `max_indel`, `graph_explained`, `alt_nodes` and `trim`.
- `--no-remap` writes every assembled contig.

Logs and robustness
- `<sample>_assembly.log` and `<sample>_remap.log` are always written, with a reason for every dropped cluster or svtig; `--keep-remap` keeps the remap GAF; tool versions are logged.
- Argument range checks; hardened tool lookup, sequence fetches and phase-file parsing.

Tests cover the command line, step time limits, reverse-strand `generate_sv_node`, cluster merging, trimming, reference paths, `--no-remap`, assembly, reference, variant and phasing code.

## 1.2.0

Baseline of the Genome Biology submission (`v1.2.0-paper-baseline`).
