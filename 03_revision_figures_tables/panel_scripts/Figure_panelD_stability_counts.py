# -*- coding: utf-8 -*-
"""Preview: Figure 6 panel D as a pathway-stability count instead of the constant 100% block.

Narrative it has to carry (Results, pathway paragraph):
  full model 19 sets -> retained in >=5 of 7 omission models -> retained in all 7 -> 0 after
  exact donor-label calibration.  That is the sentence "pathway-level findings are
  therefore reported as exploratory", drawn.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"D:\Deep harness\revision_support\fig4d"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.pad_inches": 0.03,
})

LEVELS = ["Full model", "\u22655 of 7\nomission models", "All 7\nomission models",
          "Exact donor-label\ncalibration"]
CD4 = [19, 17, 15, 0]
CD8 = [19, 13, 8, 0]
C4, C8 = "#4E79A7", "#DD8452"          # categorical pair, outside the NES red-blue ramp
GREY, DARK, GRID = "#8e98a6", "#3f4650", "#e9ecf1"

FS_TICK, FS_VAL, FS_TITLE, FS_SUB, FS_NOTE = 9.0, 9.0, 11.0, 8.4, 7.6
W_IN, H_IN = 5.40, 4.30
DPI = 400

fig = plt.figure(figsize=(W_IN, H_IN))
ax = fig.add_axes([0.125, 0.255, 0.855, 0.575])

x = range(len(LEVELS))
bw = 0.36
b1 = ax.bar([i - bw / 2 for i in x], CD4, width=bw, color=C4, edgecolor=DARK,
            linewidth=0.5, label="CD4 NC", zorder=3)
b2 = ax.bar([i + bw / 2 for i in x], CD8, width=bw, color=C8, edgecolor=DARK,
            linewidth=0.5, label="CD8 NC", zorder=3)

for bars, vals in ((b1, CD4), (b2, CD8)):
    for rect, v in zip(bars, vals):
        ax.text(rect.get_x() + rect.get_width() / 2, v + 0.45, str(v),
                ha="center", va="bottom", fontsize=FS_VAL, color=DARK, zorder=4)

# reference: the starting count
ax.axhline(19, color=GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)

ax.annotate("no set survives\nexact calibration",
            xy=(2.99, 0.3), xytext=(2.87, 4.8), fontsize=7.6, color="#B03A2E",
            ha="center", va="bottom", linespacing=1.3,
            arrowprops=dict(arrowstyle="-", color="#B03A2E", lw=0.8,
                            connectionstyle="arc3,rad=0.2"))

ax.set_xticks(list(x))
ax.set_xticklabels(LEVELS, fontsize=FS_TICK, color=DARK, linespacing=1.35)
ax.set_ylim(0, 22.5)
ax.set_yticks([0, 5, 10, 15, 20])
ax.set_ylabel("Hallmark sets with FDR < 0.05", fontsize=FS_TICK, color=DARK, labelpad=4)
ax.tick_params(axis="both", length=2.5, width=0.7, colors=DARK, labelsize=FS_TICK)
ax.yaxis.grid(True, color=GRID, lw=0.8, zorder=0)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color("#b9c0c9")
    ax.spines[s].set_linewidth(0.8)

leg = ax.legend(frameon=False, fontsize=FS_TICK, loc="upper right",
                bbox_to_anchor=(1.005, 0.995), handlelength=1.1, handleheight=0.85,
                borderpad=0.0, labelspacing=0.35)
for t in leg.get_texts():
    t.set_color(DARK)

fig.text(0.5, 0.965, "Pathway stability across donor resampling",
         ha="center", va="top", fontsize=FS_TITLE, fontweight="bold", color=DARK)
fig.text(0.5, 0.902, "of 50 Hallmark sets tested in each population",
         ha="center", va="top", fontsize=FS_SUB, color=GREY)
fig.text(0.5, 0.035,
         "Counts of Hallmark sets with FDR < 0.05; exact calibration enumerates all 35 donor-label allocations.",
         ha="center", va="bottom", fontsize=FS_NOTE, color=GREY)

for ext, kw in (("png", dict(dpi=DPI)), ("pdf", {})):
    p = os.path.join(OUT, "Figure4D_stability_count_preview.%s" % ext)
    fig.savefig(p, **kw)
    print("saved", p, os.path.getsize(p))
plt.close(fig)
