#!/usr/bin/env python3
"""Remove Figure S8 (NTI_aw diagnostic) from Analysis_pipeline.qmd
and renumber the remaining Figure S9 (soil heatmap) to S8.
"""
from pathlib import Path

src = Path("Analysis_pipeline.qmd")
text = src.read_text()
n = 0

# 1) Drop the sentence pointing to Fig. S8
sentence = " Model diagnostics for the NTI_aw model are available in Fig. S8."
if sentence in text:
    text = text.replace(sentence, "")
    n += 1
    print("[1/3] Removed Fig. S8 sentence")
else:
    print("[1/3] Sentence not found (may already be removed)")

# 2) Drop the image + caption block for Figure S8
block = (
    "![](figures/Fig_S9_NTI_aw_diagnostics.png)\n"
    "\n"
    "**Figure S8**: Diagnostic plots for the NTI_aw linear model.\n"
    "\n"
)
if block in text:
    text = text.replace(block, "")
    n += 1
    print("[2/3] Removed Figure S8 image + caption")
else:
    print("[2/3] Figure S8 block not found")

# 3) Renumber the soil heatmap caption: S9 -> S8
old_caption = "**Figure S9**: Correlation heatmap"
new_caption = "**Figure S8**: Correlation heatmap"
if old_caption in text:
    text = text.replace(old_caption, new_caption)
    n += 1
    print("[3/3] Renamed Figure S9 -> Figure S8 (soil heatmap)")
else:
    print("[3/3] Figure S9 caption not found")

src.write_text(text)
print(f"\nDone. Applied {n}/3 changes to Analysis_pipeline.qmd")