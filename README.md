# Ptarmigan_fungi

ITS2 metabarcoding of the willow ptarmigan (*Lagopus lagopus*) dung mycobiome.
Lifjellet, Lierne, Norway; winter and summer, 2022–2024.

**Supplementary appendix:** <https://ssangulo.github.io/Ptarmigan_fungi/>

## Layout

- `Scripts/` — the pipeline, numbered 1–10: demultiplexing, DADA2 + LULU, UNITE taxonomy,
  data preparation, community composition, diversity, GLLVM, functional guilds, UNITE
  species-hypothesis matching, HMSC.
- `Supplementary/` — `Supplementary_Appendix.qmd` (appendix source, and the source of truth
  for the main-text figures) with the `figures/` and `tables/` it reads.

Sequence data and intermediate workspaces are not in this repository.

## Rendering the appendix

    cd Supplementary
    quarto render Supplementary_Appendix.qmd --to html
    quarto render Supplementary_Appendix.qmd --to html --profile docs \
      --output Appendix_documentation.html

R 4.3.3; see script headers for packages.
