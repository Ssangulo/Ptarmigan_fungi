#!/usr/bin/env python3
"""Report, for every *-source pointer chunk, whether it sits inside a
docs-only block (`.content-visible when-profile="docs"`) -- it must not, or
the pointer would be missing from the submitted appendix."""
import re, sys
from pathlib import Path
L = (Path(__file__).resolve().parent.parent / "Supplementary/Supplementary_Appendix.qmd").read_text().split("\n")
stack, bad, in_code = [], 0, False
for i, l in enumerate(L, 1):
    if l.startswith("```"):
        m = re.match(r"^```\{r ([a-z0-9-]+-source)\}$", l)
        if m:
            docs = any(stack)
            print(f"{m.group(1):16s} line {i:5d}  {'DOCS-ONLY  <-- wrong' if docs else 'both builds'}")
            bad += docs
        in_code = not in_code if (l == "```" or l.startswith("```{")) else in_code
        continue
    if in_code:
        continue
    m = re.match(r"^(:{3,})\s*(.*)$", l)
    if m:
        if m.group(2):                       # opening fence with attributes
            stack.append('when-profile="docs"' in m.group(2))
        elif stack:                          # bare fence closes innermost div
            stack.pop()
sys.exit(1 if bad else 0)
