#!/usr/bin/env bash
#
# Render the submission appendix and publish it to GitHub Pages.
#
#   cd Scripts_server/Supplementary && ./publish_pages.sh
#
# Result: https://ssangulo.github.io/Ptarmigan_fungi/
#
# The site is the DEFAULT (submission) profile only -- the documentation build
# stays local. It is published to an orphan `gh-pages` branch that is rebuilt as
# a single commit every time, so re-publishing never accumulates layers in the
# repository and `main` is never touched.
#
# Why it renders from a throwaway copy of the .qmd:
#   The committed .qmd sets `embed-resources: true`, which inlines all 31 figures
#   and yields a 32 MB single file -- correct for emailing a reviewer, poor as a
#   web page. `--output-dir` alone does not change that, and `-M embed-resources:
#   false` does NOT work: quarto puts `-M` values in top-level metadata, where
#   `format.html.embed-resources` still overrides them (verified 2026-09-09,
#   quarto 1.9.38). So the script copies the .qmd to `index.qmd` in this same
#   directory (same directory, so the relative figures/ and tables/ paths in the
#   chunks still resolve), flips that one line, renders, and deletes the copy.
#   The tracked .qmd is never modified. Naming the copy `index.qmd` also makes
#   quarto name the support directory `index_files/` rather than something
#   starting with a dot, which GitHub Pages would not serve.
#   Result: 655 KB of HTML plus figures fetched as the reader scrolls.
#
# Quarto copies exactly the figures the document references into the output
# directory, so there is deliberately no `cp -a figures` step here; the
# verification below is what guards against a missed asset.

set -euo pipefail

cd "$(dirname "$0")"
SUPP="$PWD"
REPO="$(git rev-parse --show-toplevel)"
REMOTE="$(git remote get-url origin)"
SRC="Supplementary_Appendix.qmd"
TMPQMD="$SUPP/index.qmd"
BUILD="$(mktemp -d)"
CLONE="$(mktemp -d)"

cleanup() { rm -f "$TMPQMD"; rm -rf "$BUILD" "$CLONE"; }
trap cleanup EXIT

# ---- 1. render the submission profile, non-embedded --------------------------
cp "$SRC" "$TMPQMD"
sed -i '0,/embed-resources: true/s//embed-resources: false/' "$TMPQMD"
grep -q 'embed-resources: false' "$TMPQMD" \
  || { echo "ERROR: could not flip embed-resources in the copy -- has the YAML changed?"; exit 1; }

conda run -n r_env quarto render "$TMPQMD" --to html --output index.html --output-dir "$BUILD"
rm -f "$TMPQMD"

touch "$BUILD/.nojekyll"

# ---- 2. verify before anything is pushed -------------------------------------
[ -s "$BUILD/index.html" ] || { echo "ERROR: no index.html produced"; exit 1; }

size=$(stat -c %s "$BUILD/index.html")
if [ "$size" -gt 5000000 ]; then
  echo "ERROR: index.html is $((size/1048576)) MB -- embed-resources was not disabled."; exit 1
fi

# Every local src=/href= must resolve to a real file in the build.
missing=0
while read -r f; do
  [ -z "$f" ] && continue
  case "$f" in data:*|http*|"#"*|mailto:*) continue;; esac
  [ -e "$BUILD/$f" ] || { echo "MISSING ASSET: $f"; missing=1; }
done < <(
  { grep -o 'src="[^"]*"'  "$BUILD/index.html" | sed 's/src="//;s/"$//'
    grep -o 'href="[^"]*"' "$BUILD/index.html" | sed 's/href="//;s/"$//'; } | sort -u
)
[ "$missing" -eq 0 ] || { echo "ERROR: unresolved assets, not publishing"; exit 1; }

# Wrong profile: the submission build must contain no documentation bins.
if [ "$(grep -c 'Documentation notes' "$BUILD/index.html")" -ne 0 ]; then
  echo "ERROR: documentation-only content leaked into the submission build"; exit 1
fi

# A mis-counted ::: fence renders as a literal one rather than failing the build.
# vegan:::.calc_rclr is the only legitimate occurrence.
strays=$(grep -o ':::[^:]' "$BUILD/index.html" | wc -l)
legit=$(grep -o 'vegan:::.calc_rclr' "$BUILD/index.html" | wc -l)
[ "$strays" -le "$legit" ] || echo "WARNING: $((strays-legit)) literal ':::' in the output -- check the profile fences"

echo "Built $(du -sh "$BUILD" | cut -f1) site; index.html $((size/1024)) KB."

# ---- 3. publish to an orphan gh-pages branch ---------------------------------
# Done in a throwaway clone so the real working tree and .git are never touched.
git clone --quiet --depth 1 --branch main "file://$REPO" "$CLONE"
cd "$CLONE"
git checkout --quiet --orphan gh-pages
git rm -rq --cached .
rm -rf ./*
cp -a "$BUILD"/. .
git add -A
git -c user.name="$(git -C "$REPO" config user.name)" \
    -c user.email="$(git -C "$REPO" config user.email)" \
    commit -qm "pages: publish submission appendix ($(date +%Y-%m-%d))"
git push -q --force "$REMOTE" gh-pages

echo "Published to gh-pages. Site: https://ssangulo.github.io/Ptarmigan_fungi/"
