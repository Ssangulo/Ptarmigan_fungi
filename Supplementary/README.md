# Supplementary/: the Supplementary Appendix

Published version: <https://ssangulo.github.io/Ptarmigan_fungi/>

| File / folder | What it is |
|---|---|
| `Supplementary_Appendix.qmd` | The appendix source. Every table and figure in it reads a committed file from `figures/` or `tables/`, using relative paths. |
| `figures/main/` | **Main-text Figures 1–5**, saved as PNG and JPG at 600 dpi plus a vector PDF. Their build code is the `fig<N>-build` chunk in appendix Section 10. There is no separate script for them. |
| `figures/` | Supplementary figures (PNG). Each is written ("staged") here by the numbered script that produces it. |
| `tables/` | Result tables (CSV) shown or cited in the appendix, staged the same way. [`repro/manifest.csv`](../repro/manifest.csv) names the producer of every file in `figures/` and `tables/`. |
| `build_main_figure.R` | Runs one build chunk of the `.qmd` by its label, e.g. `fig1-build`. |
| `hill_size_helper.R` | Size-standardised Hill numbers (iNEXT.3D), sourced by `fig4-build`. `Scripts/6_diversity_analyses.R` carries a verbatim copy. |
| `publish_pages.sh` | Renders the submission version and publishes it to GitHub Pages (the `gh-pages` branch). |
| `appendix_qr.{png,svg,pdf}` | QR code linking to the published appendix. |

## Two versions from one source

The `.qmd` renders two ways, using Quarto profiles:

```bash
cd Supplementary
quarto render Supplementary_Appendix.qmd --to html                 # submission -> Appendix_HTML.html
quarto render Supplementary_Appendix.qmd --to html --profile docs \
  --output Appendix_documentation.html                             # documentation (full record)
```

Content wrapped in `:::: {.content-visible when-profile="docs"}` appears only in the
documentation version. The rendered HTML files are not committed.

## Rebuilding a main-text figure

Run from this folder, so the chunk's relative paths resolve the same way they do at render
time:

```bash
Rscript build_main_figure.R fig1-build      # fig1-build ... fig5-build
```

## Naming

Main-text figures are `figures/main/Fig<N>_<slug>.*`. Supplementary figures and tables keep the
file name their producing script gives them; the appendix caption carries the S-number.
