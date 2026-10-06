# CLAUDE.md proposal — APPLIED 2026-10-06

The proposed edits were refreshed and applied to `/home/daniel/CLAUDE.md` on 2026-10-06, after
`exp/repro-audit` was merged (`474387f`). They live there now; this file is kept only as a pointer.
Where to find them in CLAUDE.md:

- `# STABLE — reproducibility harness (repro/)`: how to run and read the harness, the mandatory
  8-core CPU cap, the claims registry and the generated pointer chunks.
- `# STABLE — code layout`: script 6 is the single producer of appendix §5.1-5.2.
- `2026-09-13` entry: `hill_size_helper.R` is now sourced only by `fig4-build`.
- `2026-09-29` HMSC entry: the probit stays the uncapped 2026-07-21 fit (capping tested, no help).
- `2026-10-05` DECISIONS entry: audit outcome, cold-refit results, probit decision, open items.
- `# CURRENT FOCUS`: script 5's object is now 55 samples; §3.1 figures stale (TODO #7).
