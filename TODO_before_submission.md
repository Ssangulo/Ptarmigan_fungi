# TODO before submission

Deferred decisions and open questions. Each entry records what was checked, why it was
deferred, and what has to be decided before the manuscript goes out.

---

## 1. Functional-guild assignment rule — decide how much to disclose

**Status:** deferred 2026-08-04. Currently named in one clause of the Table S-HMSCa caption;
no other text changes made. **The rule itself is deliberately unchanged.**

**Where the rule lives:** `Scripts/10_hmsc.R:193-213`. Script 8 reuses those lines *verbatim*
(commit d8b94f3) so appendix §9 and §10 classify identically. Any change here must be made in
both places or the two sections desynchronise.

**What it does:** priority `dung_saprotroph > plant_associated > pathotroph > dark_unassigned >
other`, applied by `grepl()` to the whole FUNGuild `Guild` string. FUNGuild strings are compound
(`|Dung Saprotroph|-Endophyte-Plant Saprotroph-Undefined Saprotroph`), with `|…|` marking
FUNGuild's own highest-confidence call. The regexes ignore the pipes, so an OTU qualifies for a
category if that guild appears **anywhere** in its string.

**Verified consequences** (all reproduced from `tables/hmsc_otu_guild.csv`, 226 modelled OTUs):

- All **14** dung saprotrophs matched via the FUNGuild string. The `COPRO_GENERA` genus clause
  added **zero** OTUs beyond it (7 overlap, union = 14) — the genus list is inert on this data.
- Only **2 of 14** (OTU1, OTU213, both *Thelebolus*) have "Dung Saprotroph" as FUNGuild's own
  highlighted call. For *Sporormiella* (OTU430/636/1206/1230) and *Coniochaeta* (OTU1225) the
  highlighted call is *Plant Saprotroph*; they qualify because "Dung Saprotroph" appears among
  4–7 listed alternatives.
- **OTU457 (*Penicillium*)** is classified `dung_saprotroph` off a 7-guild string whose
  highlighted call is "Undefined Saprotroph". 6 of the 14 have no genus assignment at all.
- **51 OTUs** match more than one category regex, so their class depends purely on priority order
  (46 plant∩patho, 14 dung∩plant, 9 dung∩patho).
- Because `is_plant` is tested before `is_patho`, **`pathotroph` (n = 16) is a residual class**:
  animal/fungal/lichen parasites only (*Malassezia*, Fungal Parasite ×7, Lichen Parasite ×4).
  It is not "plant pathogens", and reading a seasonal ecology off it is loose.
- `is_dark` requires `is.na(genus) AND is.na(guild)`, so the 4 OTUs with a guild but no genus and
  the 10 with a genus but no guild land in **`other`** (n = 67): 41 "Undefined Saprotroph"
  basidiomycete yeasts (*Piskurozyma* ×20, *Phenoliferia* ×5, *Vishniacozyma*, *Mrakia*,
  *Bannozyma*, *Pseudohyphozyma*) plus 8 lichenized taxa. Script 8 already splits this for display
  ("Other assigned" vs "Unassigned / dark", reconciled by `guild_scheme_crosswalk.csv`); §10 does
  not, because there it is a model trait rather than a read-share category.

**Decide before submission:**

1. Whether the methods text should state that guild membership is "guild appears anywhere in the
   compound FUNGuild string", rather than FUNGuild's highest-confidence call — a reviewer who
   knows FUNGuild will ask, and *Penicillium* as a dung saprotroph is the visible edge case.
2. Whether `pathotroph` should be relabelled (e.g. "animal/fungal/lichen parasites") so it is not
   read as plant pathogens.
3. Whether to run the highest-confidence-call-only variant as a one-line sensitivity check.
   Cheap: it is a regex change on `hmsc_otu_guild.csv`, no refit. Note script 8 already reports a
   related robustness check (COPRO_GENERA-genus-only: summer/winter dung share 45.7 %/8.5 %,
   p = 0.00094), which is reassuring for §9 but does not cover §10's trait vector.

**Not urgent because:** the guild vector enters Section 10 only as a *trait*, and the corrected
Gamma analysis finds guild does **not** credibly modify the Season response — so no Section 10
conclusion currently rests on the boundary cases above. It would matter if a later draft revived
a guild → season claim.

---

## 2. Rename script 8's figure so it obeys the `S`-prefix rule

**Status:** deferred 2026-08-05, at the same time as appendix Section 11 (main-text figures) was
created. Purely a filename/reference change — no analysis is affected.

**The rule (now in CLAUDE.md):** every *supplementary* figure filename starts with `S`; *main-text*
figures are `Fig<N>_<slug>` and live in `figures/main/`.

**The violation:** `8_functional_guilds.R` writes `plots/Fig1_guild_season.{png,pdf}` and stages it
to `Supplementary/figures/`. Despite the `Fig1_` prefix it is the **supplementary** Figure S-GUILDa
of Section 9.2, not a main-text figure. Main-text Figure 1 is now a different figure
(`figures/main/Fig1_composition_function.png`, appendix §11.1), so the two currently collide in
name only.

**What the rename requires:**

1. Edit the `fig_png` / `fig_pdf` / `figs` names in `8_functional_guilds.R` Sections 5 and 7
   (e.g. to `S_GUILDa_guild_season.{png,pdf}`).
2. Re-run `conda run -n r_env Rscript 8_functional_guilds.R` (it re-stages into `Supplementary/`;
   its Section 2d guard means it will not clobber `funguild_otu` / `fungaltraits_otu`).
3. Update the `![](figures/Fig1_guild_season.png)` reference in appendix §9.2.
4. Delete the old-named files from `plots/` and `Supplementary/figures/`.
5. Drop the parenthetical "legacy exception" note from appendix §11's file-naming paragraph.

**Not urgent because:** it is cosmetic, and both figures are already unambiguously labelled in the
prose (Figure S-GUILDa vs Figure 1) and separated by directory.
