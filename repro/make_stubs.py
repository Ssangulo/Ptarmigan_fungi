#!/usr/bin/env python3
"""repro/make_stubs.py -- (re)write the pointer-stub chunks in the appendix.

Each script-backed appendix section carries one collapsed, non-evaluated chunk
labelled `<id>-source` that says which script produces the section's figures
and tables, what it reads, what it writes and how to re-run it -- instead of
duplicating the script's code. The `# Output :` lines are generated from
repro/manifest.csv, so a stub cannot list a file its script does not write
(repro/check_manifest.R verifies this on every run).

Usage (from the repo root):  python3 repro/make_stubs.py [--run-id RUN_ID]
  Idempotent: every existing `*-source` chunk is removed, then each stub is
  inserted immediately before its anchor line (or before the docs-only wrapper
  that opens just above it). Check placement with repro/stub_visibility.py.
"""
import csv, re, sys, argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
QMD = ROOT / "Supplementary" / "Supplementary_Appendix.qmd"
MANIFEST = ROOT / "repro" / "manifest.csv"

# id, source, sections described, inputs, run command, runtime note, anchor
# (the stub goes just before the first line that starts with `anchor`), and
# an optional filter on the manifest `sections` column so a script that serves
# several sections lists only the files of this one.
STUBS = [
    dict(id="sec1", source="Scripts/4_data_prep.R",
         summary="Source of the analysis objects: Scripts/4_data_prep.R",
         what="Sections 1-9 -- builds the analysis objects every section below reads",
         inputs=["nopoolps.dada2.rds / poolps / pspoolps (2_DADA2_lulu.R + 3_taxonomic_assignment.R)",
                 "fullps_nopool.rds (PCR-replicate level), rg2.*.tax.rds (sensitivity)",
                 "matches_out_taxonomy_{nopool,pool,pspool}.csv (PlutoF non-fungal screen)"],
         outputs_free=["eco_analysis.RData (alldat, alldat.rfy, alldat_full, alldat_rg2_sensitivity;",
                       "  kept outside the repository with the other large data objects)",
                       "phyloseq_comparison_summary.csv (read by the Section 2.2 build chunk)"],
         run="conda run -n r_env Rscript Scripts/4_data_prep.R", runtime="< 1 min",
         anchor="```{r load-objects}", section_filter=None),
    dict(id="sec3-1", source="Scripts/5_community_composition.R",
         summary="Source: Scripts/5_community_composition.R -- inputs, outputs, how to re-run",
         what="Section 8 (taxonomic bar plots, Figures S2-S5); Section 9 (assignment summary)",
         inputs=["eco_analysis.RData: alldat_full$nopool, PCR replicates summed in the",
                 "  script's Section 1 (its own object: no depth filter, see Section 3.1)"],
         run="conda run -n r_env Rscript Scripts/5_community_composition.R", runtime="< 1 min",
         anchor="### 3.2 ", section_filter=None),
    dict(id="sec3-4", source="Scripts/9_dark_taxa_SH_matching.R",
         summary="Source: Scripts/9_dark_taxa_SH_matching.R -- inputs, outputs, how to re-run",
         what="Sections 1-7 (UNITE Species Hypothesis matching of every OTU)",
         inputs=["eco_analysis.RData: alldat$nopool reference sequences",
                 "sh_matching container output outdata/source_1.zip (Scripts/9_setup_sh_matching.sh",
                 "  builds the container; re-running it takes ~90 min, the script reuses the zip)",
                 "PlutoF live-database run models/source_67536.zip (external, not regenerable)"],
         run="conda run -n r_env Rscript Scripts/9_dark_taxa_SH_matching.R", runtime="< 1 min (container output reused)",
         anchor="#### Documentation notes — Section 3", section_filter=None),
    dict(id="sec5", source="Scripts/6_diversity_analyses.R",
         summary="Source: Scripts/6_diversity_analyses.R (Sections 1-2) -- inputs, outputs, how to re-run",
         what="Sections 1-2 (pooled iNEXT3D, per-sample and per-replicate Hill numbers, Tables S8-S9)",
         inputs=["eco_analysis.RData: alldat$nopool (52 samples) and alldat_full$nopool",
                 "size-standardised Hill numbers: verbatim copy of Supplementary/hill_size_helper.R"],
         run="S6_STOP_AFTER=H2 conda run -n r_env Rscript Scripts/6_diversity_analyses.R",
         runtime="~5 min (S6_STOP_AFTER=H2 stops after Section 2; the full script incl. Part B takes ~45 min)",
         anchor="```{r sec5-build}", section_filter=r"^5"),
    dict(id="sec6", source="Scripts/7_gllvm.R",
         summary="Source: Scripts/7_gllvm.R (Sections 1-7) -- inputs, outputs, how to re-run",
         what="Sections 1-7 (H1 GLLVM: model comparison, per-OTU Season coefficients, ordination, genus summaries)",
         inputs=["models/fit_nb_2.rds, fit_nb_2_2.rds, fit_poisson.rds (saved fits, 2026-07-02)",
                 "eco_analysis.RData.bak_20260701_partial (the pre-depth-fix object they were fitted on;",
                 "  see the caveat above -- these are not refitted on the canonical object)"],
         run="conda run -n r_env Rscript Scripts/7_gllvm.R", runtime="~1 min with the saved fits",
         anchor="#### Documentation notes — Section 6", section_filter=r"^6"),
    dict(id="sec7-1", source="Scripts/6_diversity_analyses.R",
         summary="Source: Scripts/6_diversity_analyses.R (Sections 3, 3b) -- inputs, outputs, how to re-run",
         what="Sections 3 and 3b (H2 mechanism: Gamma brms models, family choice, sensitivity, robustness)",
         inputs=["eco_analysis.RData: alldat$nopool; per-sample Hill numbers from the script's Section 2",
                 "plant_ITS/phyloseq_plant_ITS.rds (paired plant ITS2 diet data)",
                 "brms fits cached as models/H2_brms_diet_richness_gamma{,_noyear}.rds"],
         run="S6_STOP_AFTER=H2MECH conda run -n r_env Rscript Scripts/6_diversity_analyses.R",
         runtime="~10 min with the cached brms fits",
         anchor="### 7.2 ", section_filter=r"^7\.1"),
    dict(id="sec7-2", source="Scripts/6_diversity_analyses.R",
         summary="Source: Scripts/6_diversity_analyses.R (Section 3c) -- inputs, outputs, how to re-run",
         what="Section 3c (H3: hurdle negative-binomial brms model, specificity contrast)",
         inputs=["eco_analysis.RData: alldat_full$nopool (PCR-replicate level, matched subset)",
                 "plant_ITS/phyloseq_plant_ITS.rds",
                 "brms fits cached as models/H3_brms_hurdle_{species,genus}.rds (~2.4 h cold)"],
         run="S6_STOP_AFTER=H3 conda run -n r_env Rscript Scripts/6_diversity_analyses.R",
         runtime="~15 min with the cached brms fits (Sections 1-3c)",
         anchor="### 7.3 ", section_filter=r"^7\.2"),
    dict(id="sec7-3", source="Scripts/7_gllvm.R",
         summary="Source: Scripts/7_gllvm.R (Section 8) -- inputs, outputs, how to re-run",
         what="Section 8 (H3 GLLVM cross-check on the same predictors)",
         inputs=["eco_analysis.RData (canonical), plant_ITS/phyloseq_plant_ITS.rds",
                 "models/H3_perOTU_plant_coef.csv from script 6 Section 3c (agreement tables)",
                 "GLLVM fits cached as models/gllvm_H3sp_*.rds"],
         run="conda run -n r_env Rscript Scripts/7_gllvm.R",
         runtime="~1 min with the cached GLLVM fits (runs Sections 1-8)",
         anchor="#### Documentation notes — Section 7", section_filter=r"^7\.3"),
    dict(id="sec8", source="Scripts/8_functional_guilds.R",
         summary="Source: Scripts/8_functional_guilds.R -- inputs, outputs, how to re-run",
         what="Sections 1-7 (FUNGuild lookup, guild rule, seasonal read shares, Figure S-GUILDa)",
         inputs=["eco_analysis.RData: alldat$nopool (52 samples, all 1129 OTUs)",
                 "funguild_input.guilds.txt (cached FUNGuild CLI run, used only as an OTU -> guild lookup)"],
         run="conda run -n r_env Rscript Scripts/8_functional_guilds.R", runtime="< 1 min",
         anchor="#### Documentation notes — Section 8", section_filter=None),
    dict(id="sec9", source="Scripts/10_hmsc.R",
         summary="Source: Scripts/10_hmsc.R -- inputs, outputs, how to re-run",
         what="Sections 1-7 (three HMSC models, phylogeny variant, PCR-replicate variance)",
         inputs=["eco_analysis.RData: alldat_full$nopool (226 OTUs present in >= 5 replicates)",
                 "tree.rds (ML ITS2 tree, Scripts/6_diversity_analyses.R Section 4)",
                 "models/gllvm_perOTU_season_coef.csv (Section 6, for the cross-check)",
                 "fits cached as models/hmsc_{clr,rep,pa,clr_phylo}_fit.rds"],
         run="OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 HMSC_BLAS_THREADS=4 HMSC_NPARALLEL=4 \\\n#          HMSC_RUN_PHYLO=1 HMSC_RUN_CV=0 conda run -n r_env Rscript Scripts/10_hmsc.R",
         runtime="~6 min with the cached fits; cold: ~5.5 h fits + ~2.5 h CV + ~5.7 h phylogeny",
         anchor="#### Documentation notes — Section 9", section_filter=None),
    dict(id="sec9-v1", source="hmsc_converge/04_v1_disagreement_table.R",
         summary="Source: hmsc_converge/04_v1_disagreement_table.R -- the superseded model's chain table",
         what="whole script (Table S-HMSCb2: per-chain values of the superseded 2026-07-20 fit)",
         inputs=["models/hmsc_v2/hmsc_clr_fit_diag.rds (saved diagnostics of the superseded fit;",
                 "  historical, not regenerable from the current scripts)"],
         run="conda run -n r_env Rscript hmsc_converge/04_v1_disagreement_table.R", runtime="seconds",
         anchor="#### Documentation notes — Section 9", section_filter=None),
]


def manifest_outputs(source, section_filter):
    out = []
    for r in csv.DictReader(open(MANIFEST)):
        if r["producer"] != source:
            continue
        if section_filter and not re.search(section_filter, r["sections"]):
            continue
        out.append(r["file"])
    return sorted(out, key=lambda f: (not f.startswith("figures/"), f))


def stub_text(s, run_id):
    files = s.get("outputs_free") or manifest_outputs(s["source"], s["section_filter"])
    L = [f"```{{r {s['id']}-source}}", "#| eval: false",
         "#| echo: !expr knitr::is_html_output()", "#| code-fold: true",
         f'#| code-summary: "{s["summary"]}"', "",
         f"# Source : {s['source']}, {s['what']}"]
    for i, x in enumerate(s["inputs"]):
        L.append(("# Input  : " if i == 0 else "#          ") + x)
    for i, x in enumerate(files):
        L.append(("# Output : " if i == 0 else "#          ") + x)
    L.append(f"# Run    : (from Scripts_server/) {s['run']}")
    L.append(f"#          {s['runtime']}")
    if run_id:
        L.append(f"# Checked: repro/logs/{run_id}/summary.md")
    L.append("```")
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", default="")
    a = ap.parse_args()
    lines, i = QMD.read_text().split("\n"), 0
    while i < len(lines):                                   # drop existing stubs
        if re.match(r"^```\{r [a-z0-9-]+-source\}$", lines[i]):
            j = next(k for k in range(i + 1, len(lines)) if lines[k] == "```")
            del lines[i:j + 2 if j + 1 < len(lines) and lines[j + 1] == "" else j + 1]
            continue
        i += 1
    for s in STUBS:
        new = stub_text(s, a.run_id)
        if True:
            idx = [k for k, l in enumerate(lines) if l.startswith(s["anchor"])]
            if not idx:
                sys.exit(f"anchor not found for {s['id']}: {s['anchor']!r}")
            k = idx[0]
            # "#### Documentation notes" bins sit inside a docs-only wrapper that
            # opens just above them; the stub must go before the wrapper or it
            # would vanish from the submission build.
            while k > 0 and (lines[k - 1].strip() == "" or lines[k - 1].startswith(":::")):
                k -= 1
                if lines[k].startswith(":::"):
                    break
            if lines[k].startswith(":::"):
                pass
            else:
                k = idx[0]
            lines[k:k] = new + [""]
    QMD.write_text("\n".join(lines))
    print(f"wrote {len(STUBS)} pointer stubs into {QMD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
