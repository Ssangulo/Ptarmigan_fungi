# repro/ — is every appendix number reproducible?

One command regenerates every table and figure the Supplementary Appendix shows,
and checks them against what is committed and against the numbers quoted in the
prose.

```bash
cd Scripts_server
./repro/run_repro.sh cached          # a few hours; fitted models reused (as copies)
./repro/run_repro.sh cold            # ~15-20 h; every model refitted from scratch
./repro/run_repro.sh checks <RUN_ID> # re-run only the checks of an earlier run
```

Each run writes `repro/logs/<YYYY-MM-DD_HHMM>_<mode>/`. Start with `summary.md`.

**CPU budget (shared server).** Cap a run with `REPRO_CORES`: it re-launches under
`taskset` (only those cores) at `nice 19`. Match the thread counts to it, e.g. for 8 cores:

```bash
REPRO_CORES=88-95 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 HMSC_BLAS_THREADS=2 HMSC_NPARALLEL=4 \
  setsid nohup ./repro/run_repro.sh cold > ~/Ptarmigan/repro_sandbox/cold_console.log 2>&1 &
```

`setsid nohup` keeps a long run alive after the terminal or Claude session closes; the
run's PID is in `repro/logs/<RUN_ID>/run.pid` (`kill -- -<PID>` stops it and its children).

## What a run does

1. Records the environment (`env.txt`) and checksums the upstream inputs
   (`inputs.sha256`). Scripts 1-3 (demultiplexing, DADA2/LULU, taxonomy) are
   **not** re-run; their outputs are pinned by these checksums instead.
2. Runs `4_data_prep.R` into a sandbox and compares the rebuilt workspace with
   the canonical `eco_analysis.RData`, object by object. **Gate:** if `alldat`,
   `alldat.rfy`, `alldat_full` or `alldat_rg2_sensitivity` differ, the run stops,
   because every later script reads the canonical workspace.
3. Runs scripts 5, 8, 9, 6, 7 and 10. Each one writes its models/plots/tables to
   the sandbox (`$REPRO_SANDBOX/<RUN_ID>/`, default `~/Ptarmigan/repro_sandbox/`)
   and stages its appendix files into **this checkout's** `Supplementary/`.
   Scripts run from frozen copies, because `Rscript` reads a file incrementally.
4. Runs the appendix build chunks (`sec2/4/5-build`, `fig1`-`fig5-build`) with
   `Supplementary/build_main_figure.R`, and `model_development/hmsc_convergence/04_v1_disagreement_table.R`.
5. Renders both appendix profiles, greps for stray `:::` fences, and checks that
   no pointer chunk sits inside a docs-only block.
6. Runs the checks below.

The canonical `~/Ptarmigan/models|plots|tables` and `eco_analysis.RData` are only
ever **read**. Because everything is staged into the checkout, the run refuses to
start unless `Supplementary/`, `Scripts/` and `model_development/` are committed, and
`git diff` afterwards shows exactly what the run changed.

## The checks

| Check | File in the log | Fails when |
|---|---|---|
| Outputs | `diff_outputs.csv` | a committed table/figure comes back different: CSVs by value (tolerance 1e-8 cached, 0.05 cold), PNGs by pixel, PDFs by rendered page |
| Manifest | `manifest_check.csv` | a file the appendix shows has no producer, its producer never names it, it was not rewritten during the run, or a pointer chunk lists a file its script does not write (or misses one) |
| Quoted numbers | `claims_result.csv` | a number in `claims.csv` no longer matches its table **or** its sentence is no longer in the qmd |
| Number sweep | `unverified_numbers.csv` | never fails; lists every other number in the prose that no table backs, for triage |

## Files you edit

- **`manifest.csv`**: one row per file under `Supplementary/tables|figures`. It names
  the producer (`Scripts/N_*.R`, `qmd:<chunk>` or `model_development/hmsc_convergence/04_*.R`), the script
  section, the appendix section and the role (`shown`, `input` = read only by a build
  chunk, `unused`, `deliverable`). A non-empty `exception` marks a documented
  reproducibility exception (reported as EXCEPTION rather than FAIL when the file
  is not regenerated). Section 6's rows carry one; see below.
- **`claims.csv`**: one row per quoted number. `expr` is R evaluated with the table
  bound to `d` (`tab("tables/x.csv")` reads a second one); `digits` and `quoted` are
  the rounding and the number as printed; `snippet` is a verbatim piece of the
  sentence. When you **rewrite prose** that contains a registered number, update its
  snippet. When a **number legitimately changes**, update `quoted` (and the prose).
- **`cached_fits.txt`**: the fitted models copied into the sandbox in cached mode.

## Pointer chunks in the appendix

Every script-backed section ends with a collapsed `<id>-source` chunk: source
script and sections, inputs, outputs, run command, runtime. They are generated:

```bash
python3 repro/make_stubs.py --run-id <RUN_ID>   # idempotent; Output lines come from manifest.csv
python3 repro/stub_visibility.py                # every pointer must render in both builds
```

Edit `STUBS` in `make_stubs.py`, not the chunks themselves.

## Known, documented exceptions

- **Section 6 (H1 GLLVM)** reproduces from the saved 2026-07-02 fits and the
  pre-depth-fix backup object, not from the canonical data. The check confirms that
  script 7 regenerates every Section 6 table from those fits exactly; it does not
  refit them.
- **`hmsc_v1_chain_disagreement.csv`** is read from the saved diagnostics of the
  superseded 2026-07-20 HMSC fit (historical record).
- **Scripts 1-3** are pinned by checksum, not re-run (~2 days from raw FASTQ).
