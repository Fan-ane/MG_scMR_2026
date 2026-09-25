import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

# =============================== DATA INPUT ==================================
# One row per pathway, in top-to-bottom figure order:
#     (pathway label, NES CD4 NC, FDR CD4 NC, NES CD8 NC, FDR CD8 NC)
# FDR = FDR_across_all_evaluable_tests; < 0.05 is drawn as a filled marker.
# Edit / replace these numbers to plot your own data.
DATA = [
    ("Tnfa Signaling via Nfkb",             2.70425037551519,  5.82443797026295e-21,  2.32062612615513,  1.97438767606764e-09),
    ("Oxidative Phosphorylation",           2.34970843259259,  7.29954052663204e-14,  2.42371679490278,  3.09342060534999e-15),
    ("Mitotic Spindle",                    -2.17397887197996,  4.78114573929654e-09, -2.00864610399126,  7.65086466766829e-06),
    ("Wnt Beta Catenin Signaling",         -1.94265426332669,  0.00142950517885907,  -1.77983209495236,  0.0182740428142702),
    ("Hypoxia",                             1.6623426414694,   0.002622362760098,     1.84871594175405,  0.00047231267413896),
    ("Inflammatory Response",               1.83834262803245,  0.000176474168082345,  1.48946342183239,  0.0388811963333092),
    ("Reactive Oxygen Species Pathway",     1.7789247722909,   0.00692490562180769,   1.82688166779151,  0.00532804852899175),
    ("P53 Pathway",                         1.79610668326817,  0.00010178544152102,   1.38378591514935,  0.0812981951694294),
    ("Myc Targets V1",                      1.40977861201043,  0.0222003683333706,    1.7210386446839,   0.000511114866058124),
    ("Epithelial Mesenchymal Transition",   1.63007531083558,  0.0132970649176946,    1.66810707463098,  0.0227938583311535),
    ("Apoptosis",                           1.50331817091373,  0.0170361284061738,    1.63303633466629,  0.00692490562180769),
    ("Glycolysis",                          1.32076705206371,  0.0834107409789942,    1.62679783618294,  0.00948239414735671),
    ("Apical Junction",                    -1.15793804905164,  0.284487385936661,    -1.61060107313798,  0.0182740428142702),
    ("Complement",                          1.56806625079612,  0.00490888266235209,   1.54993542261843,  0.0170361284061738),
    ("Dna Repair",                          1.55511598260884,  0.00706326342555587,   1.50091892904213,  0.0210331548391597),
    ("Kras Signaling Up",                   1.52159339070823,  0.0121640716089936,    1.53854262709496,  0.0285956652778816),
    ("Uv Response Up",                      1.44201081590792,  0.0326774719206393,    1.52453838742761,  0.0210331548391597),
    ("Uv Response Dn",                     -1.48551465857868,  0.0222151685373483,   -1.51631956407459,  0.0388811963333092),
    ("G2m Checkpoint",                     -1.50770683093112,  0.00939613868092737,  -1.49727369622999,  0.0182740428142702),
    ("Estrogen Response Late",              1.44327970388187,  0.0246592073980635,    1.43174816575428,  0.0586439180457674),
    ("Allograft Rejection",                 1.436753,          0.0227939,             1.416863,          0.0395517),
]

# Column headers shown in the figure (left to right, matches the tuple order).
GROUP_LABELS = ["CD4 NC", "CD8 NC"]

# Optional: read the raw numbers from a csv instead of DATA.
# Columns required: pathway_label, cell_type, NES, FDR_across_all_evaluable_tests
# (cell_type values CD4_NC / CD8_NC). DATA still supplies the row order.
CSV_PATH = None          # e.g. "Figure4C_Hallmark_all50.csv"
# =============================================================================

# --------------------------------------------------------------- settings ---
ASPECT = 16.0 / 9.0     # content width / height.  Use 1.0 for a square canvas.
FS_LAB = 8.0            # pathway label size (pt)
FS_HDR = 8.4            # group header size (pt)
K_MARK = 5.4            # marker size: s = (|NES| * K_MARK) ** 2
ROWS_PER_BLOCK = 11     # 21 pathways -> blocks of 11 + 10, side by side
GAP_BLOCK = 0.26        # inches between a block and the next block's labels
GAP_LEG = 0.30          # inches between the last matrix and the legend column
LEG_W = 1.06            # width of the right-hand legend column (inches)
CBAR_W = 0.16           # colour bar width (inches)
BOT = 0.16              # inches below the matrices
TITLE_TOP = 0.13        # inches above the title
V_MIN, V_MAX = -2.5, 2.5
FDR_CUT = 0.05
DPI = 400

TITLE = "Hallmark pathway activity"
SUBTITLE = ("all 21 Hallmark sets reaching FDR < 0.05 in at least one population; "
            "pathways 1\u201311 (left) and 12\u201321 (right), "
            "in the original order")

CMAP = plt.get_cmap("RdBu_r")
NORM = Normalize(V_MIN, V_MAX)
BAND = "#f1f3f6"        # alternating row background
GREY = "#8e98a6"
DARK = "#3f4650"

BASE = os.path.dirname(os.path.abspath(__file__))
OUTD = r"D:\Deep harness\revision_support\fig4c"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42,        # embed TrueType so text stays editable in Illustrator
    "ps.fonttype": 42,
    "savefig.pad_inches": 0.03,
})


# ----------------------------------------------------------------- helpers ---
def load_data():
    """Return (row_order, nes, fdr) with one column per entry of GROUP_LABELS."""
    order = [row[0] for row in DATA]

    if CSV_PATH:
        df = pd.read_csv(os.path.join(BASE, CSV_PATH))
        keep = df[df["pathway_label"].isin(order)]
        cells = ["CD4_NC", "CD8_NC"]
        nes = keep.pivot(index="pathway_label", columns="cell_type", values="NES")
        fdr = keep.pivot(index="pathway_label",
                         columns="cell_type",
                         values="FDR_across_all_evaluable_tests")
        nes = nes.reindex(index=order, columns=cells)
        fdr = fdr.reindex(index=order, columns=cells)
        nes.columns = GROUP_LABELS
        fdr.columns = GROUP_LABELS
    else:
        nes = pd.DataFrame({g: [row[1 + 2 * i] for row in DATA]
                            for i, g in enumerate(GROUP_LABELS)}, index=order)
        fdr = pd.DataFrame({g: [row[2 + 2 * i] for row in DATA]
                            for i, g in enumerate(GROUP_LABELS)}, index=order)

    bad = [r for r in order if nes.loc[r].isna().any() or fdr.loc[r].isna().any()]
    if bad:
        raise SystemExit("missing NES/FDR for: %s" % bad)
    if not ((fdr.values > 0) & (fdr.values <= 1)).all():
        raise SystemExit("FDR values must lie in (0, 1]")
    return order, nes, fdr


def text_w_in(s, fs, weight="normal"):
    """Rendered width of a string, in inches. Used to size the columns from the
    real font metrics instead of guessing."""
    f = plt.figure(figsize=(2, 2), dpi=200)
    t = f.text(0.1, 0.5, s, fontsize=fs, fontweight=weight)
    f.canvas.draw()
    w = t.get_window_extent(renderer=f.canvas.get_renderer()).width / f.dpi
    plt.close(f)
    return w


def draw_bubbles(ax, rows, nes, fdr, k):
    """One block: pathway labels on the left, one dot column per group."""
    n = len(rows)
    for i in range(n):
        if i % 2 == 0:
            ax.add_patch(Rectangle((-0.5, i - 0.5), 2.0, 1.0, facecolor=BAND,
                                   edgecolor="none", zorder=0))
    ax.axvline(0.5, color="#c3c9d2", lw=0.8, zorder=1)
    for i, r in enumerate(rows):
        for j, g in enumerate(GROUP_LABELS):
            v = float(nes.loc[r, g])
            f = float(fdr.loc[r, g])
            s = (abs(v) * k) ** 2
            if f < FDR_CUT:
                ax.scatter(j, i, s=s, marker="o", facecolor=CMAP(NORM(v)),
                           edgecolor=DARK, linewidths=0.55, zorder=3)
            else:
                ax.scatter(j, i, s=s, marker="o", facecolor="white",
                           edgecolor=CMAP(NORM(v)), linewidths=1.25, zorder=3)
    ax.set_xlim(-0.5, len(GROUP_LABELS) - 0.5)
    ax.set_ylim(n - 0.5, -0.5)              # row 0 on top
    ax.set_yticks(range(n))
    ax.set_yticklabels(rows, fontsize=FS_LAB)
    ax.set_xticks(range(len(GROUP_LABELS)))
    ax.set_xticklabels(GROUP_LABELS, fontsize=FS_HDR, fontweight="bold")
    ax.xaxis.set_ticks_position("top")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=0, pad=4)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.margins(0)


def add_colorbar(fig, x, y, w, h, fw, fh):
    cax = fig.add_axes([x / fw, y / fh, w / fw, h / fh])
    sm = ScalarMappable(norm=NORM, cmap=CMAP)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("NES", fontsize=8.6, labelpad=3)
    cb.ax.tick_params(labelsize=7.6, length=2, width=0.6)
    cb.outline.set_linewidth(0.5)
    cb.set_ticks([-2, -1, 0, 1, 2])


def add_legend_block(fig, x, y, w, h, fw, fh, vals, k):
    """|NES| size legend + filled/hollow legend. Axis units are inches."""
    lax = fig.add_axes([x / fw, y / fh, w / fw, h / fh])
    lax.set_xlim(0, w)
    lax.set_ylim(0, h)
    lax.axis("off")
    cx, tx = 0.17, 0.33
    lax.text(0.0, h, "|NES|", fontsize=8.6, fontweight="bold", va="top",
             ha="left", color=DARK)
    for v, yy in zip(vals, [h * 0.66, h * 0.47, h * 0.28]):
        lax.scatter(cx, yy, s=(v * k) ** 2, c=GREY, edgecolor=DARK,
                    linewidths=0.5, clip_on=False, zorder=3)
        lax.text(tx, yy, f"{v:.1f}", fontsize=7.8, va="center", color=DARK)
    lax.scatter(cx, h * 0.115, s=46, c=GREY, edgecolor=DARK, linewidths=0.5,
                clip_on=False, zorder=3)
    lax.text(tx, h * 0.115, f"FDR < {FDR_CUT:g}", fontsize=7.4, va="center", color=DARK)
    lax.scatter(cx, h * 0.025, s=46, c="white", edgecolor=GREY, linewidths=1.25,
                clip_on=False, zorder=3)
    lax.text(tx, h * 0.025, f"FDR $\\geq$ {FDR_CUT:g}", fontsize=7.4, va="center",
             color=DARK)


# ------------------------------------------------------------------- build ---
def main():
    order, nes, fdr = load_data()
    print("input: %d pathways x %d groups" % (len(order), len(GROUP_LABELS)))

    blocks = [order[i:i + ROWS_PER_BLOCK]
              for i in range(0, len(order), ROWS_PER_BLOCK)]

    # --- geometry: derive the canvas from the real text metrics ---------------
    lab_w = max(text_w_in(r, FS_LAB) for r in order) + 0.09
    hdr_w = max(text_w_in(g, FS_HDR, "bold") for g in GROUP_LABELS)
    ax_w = 2.0 * (hdr_w + 0.06)
    content_w = (len(blocks) * (lab_w + ax_w)
                 + (len(blocks) - 1) * GAP_BLOCK + GAP_LEG + LEG_W)
    content_h = content_w / ASPECT
    fh = content_h + TITLE_TOP + BOT
    fw = content_w + 0.10
    head = 0.16 + (0.16 if SUBTITLE else 0.0) + 0.06 + 0.14
    ax_h = fh - TITLE_TOP - head - BOT

    fig = plt.figure(figsize=(fw, fh))

    def rect(x, y, w, h):
        return [x / fw, y / fh, w / fw, h / fh]

    x = lab_w
    for rows in blocks:
        draw_bubbles(fig.add_axes(rect(x, BOT, ax_w, ax_h)), rows, nes, fdr, K_MARK)
        x += ax_w + GAP_BLOCK + lab_w

    # right-hand column: colour bar stacked above the size / significance legend
    leg_x = content_w - LEG_W
    cb_h, gap_l, lg_h = 1.55, 0.28, 2.05
    y0 = BOT + max(0.0, (ax_h - (cb_h + gap_l + lg_h)) / 2.0)
    add_colorbar(fig, leg_x + 0.15, y0 + lg_h + gap_l, CBAR_W, cb_h, fw, fh)
    add_legend_block(fig, leg_x, y0, LEG_W, lg_h, fw, fh, [1.5, 2.0, 2.5], K_MARK)

    fig.text(0.5, (fh - TITLE_TOP) / fh, TITLE, ha="center", va="top",
             fontsize=11, fontweight="bold")
    if SUBTITLE:
        fig.text(0.5, (fh - TITLE_TOP - 0.20) / fh, SUBTITLE, ha="center",
                 va="top", fontsize=7.4, color=DARK)

    # ------------------------------------------------------------------ save --
    os.makedirs(OUTD, exist_ok=True)
    for ext, kw in (("png", dict(dpi=DPI)), ("pdf", {})):
        p = os.path.join(OUTD, "Figure4C_A_bubble21.%s" % ext)
        fig.savefig(p, bbox_inches="tight", **kw)
        print("saved", p)
    plt.close(fig)

    from PIL import Image
    w, h = Image.open(os.path.join(OUTD, "Figure4C_A_bubble20.png")).size
    print("layout: lab_w=%.2f ax_w=%.2f ax_h=%.2f  row=%.1fpt  content=%.2fx%.2f"
          % (lab_w, ax_w, ax_h, ax_h * 72 / ROWS_PER_BLOCK, content_w, content_h))
    print("output: %dx%d px  ratio=%.3f (target %.3f)" % (w, h, w / h, ASPECT))


if __name__ == "__main__":
    main()
