# -*- coding: utf-8 -*-
"""Figure 1 rebuilt in the original graphical-abstract style, v4 (layout defects fixed).

Design at final print size: 1 unit = 1 pt, canvas 510 x 690 pt = 180 x 243 mm at 600 dpi.
All text is >= 7 pt at print size; graphics are 300-dpi assets cropped from the authors' original panels.
"""
import os, math, random
from PIL import Image, ImageDraw, ImageFont
import fitz

ASSETS = r"D:\Deep harness\dump\fig1_assets"
OUT = r"D:\Deep harness\revision_support"
DPI = 600
K = DPI / 72.0
W_PT, H_PT = 510.0, 690.0
W, H = int(W_PT * K), int(H_PT * K)

REG = r"C:\Windows\Fonts\arial.ttf"
BOLD = r"C:\Windows\Fonts\arialbd.ttf"
ITAL = r"C:\Windows\Fonts\ariali.ttf"
_cache = {}


def font(pt, bold=False, italic=False):
    key = (round(pt, 1), bold, italic)
    if key not in _cache:
        _cache[key] = ImageFont.truetype(BOLD if bold else (ITAL if italic else REG),
                                         int(round(pt * K)))
    return _cache[key]


img = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(img)
P = lambda pt: int(round(pt * K))


def text(x, y, s, pt=7.5, bold=False, fill="#1b1b1b", anchor="la", italic=False):
    d.text((P(x), P(y)), s, font=font(pt, bold, italic), fill=fill, anchor=anchor)


def rounded(x0, y0, x1, y1, r=7, fill="white", outline="#C9CED6", width=0.8):
    d.rounded_rectangle([P(x0), P(y0), P(x1), P(y1)], radius=P(r), fill=fill,
                        outline=outline, width=max(1, P(width)))


def arrow(x0, y0, x1, y1, colour="#2B2B2B", w=1.6, head=4.0):
    d.line([P(x0), P(y0), P(x1), P(y1)], fill=colour, width=max(1, P(w)))
    ang = math.atan2(y1 - y0, x1 - x0)
    for s in (+1, -1):
        a = ang + s * math.radians(150)
        d.line([P(x1), P(y1), P(x1 + head * math.cos(a)), P(y1 + head * math.sin(a))],
               fill=colour, width=max(1, P(w)))


def paste(name, x, y, w_pt, max_h=None, centre_x=None):
    p = os.path.join(ASSETS, name + ".png")
    if not os.path.exists(p):
        print("  !! missing asset", name); return 0
    im = Image.open(p).convert("RGB")
    w_px = P(w_pt); h_px = int(im.size[1] * w_px / im.size[0])
    if max_h and h_px > P(max_h):
        h_px = P(max_h); w_px = int(im.size[0] * h_px / im.size[1])
    im = im.resize((w_px, h_px), Image.LANCZOS)
    if centre_x is not None:
        x = centre_x - (w_px / K) / 2
    img.paste(im, (P(x), P(y)))
    return h_px / K


COLS = [("A", 8.0, "Main data source", "#1F4E99"),
        ("B", 176.0, "Primary analyses", "#7030A0"),
        ("C", 344.0, "Follow-up analyses", "#375623")]
CW = 158.0
ROWS = [(44.0, 248.0), (256.0, 460.0), (468.0, 672.0)]

for tag, x, title, colour in COLS:
    d.rounded_rectangle([P(x), P(8), P(x + CW), P(34)], radius=P(6), fill=colour)
    text(x + CW / 2, 14.6, "%s. %s" % (tag, title), pt=10, bold=True, fill="white", anchor="ma")
arrow(8 + CW + 1, 21, 175, 21, w=2.0, head=4.5)
arrow(176 + CW + 1, 21, 343, 21, w=2.0, head=4.5)

TITLES = {(0, 0): "Single-cell eQTL datasets", (0, 1): "GWAS summary statistics",
          (0, 2): "Study scope and abbreviations",
          (1, 0): "Two-sample Mendelian randomization", (1, 1): "Bayesian colocalization",
          (1, 2): "Single-cell transcriptomics",
          (2, 0): "Phenome-wide association analysis", (2, 1): "Drug annotation",
          (2, 2): "Prioritized genes and programs"}
for ci, (_, x, _, _) in enumerate(COLS):
    for ri, (y0, y1) in enumerate(ROWS):
        rounded(x, y0, x + CW, y1)
        text(x + 9, y0 + 7, TITLES[(ci, ri)], pt=8.5, bold=True)
        if ri < 2:
            arrow(x + CW / 2, y1 + 1.5, x + CW / 2, y1 + 7.5, colour="#6A3FA0", w=1.8, head=4.2)

# ---------------- A1 ----------------
x, y0 = 8.0, 44.0
paste("a1_motif", x + 12, y0 + 20, 46, max_h=30)
text(x + 64, y0 + 22, "cis-eQTLs from OneK1K", pt=7.5)
text(x + 64, y0 + 32, "PBMC of healthy donors", pt=7.5)
text(x + 12, y0 + 56, "Instrument selection", pt=7.5, bold=True)
text(x + 12, y0 + 66, "cis-eQTL P < 5 \u00d7 10\u207b\u2078;  F > 10", pt=7.5)
text(x + 12, y0 + 76, "conditionally independent records", pt=7.5)
paste("a1_cells", x + 12, y0 + 90, 134)
text(x + 12, y0 + 118, "14 peripheral blood immune cell types", pt=7.5)
text(x + 12, y0 + 130, "7,328 post-MHC gene-cell type tests", pt=7.5)
text(x + 12, y0 + 140, "39 global-BH significant pairs", pt=7.5)
text(x + 12, y0 + 150, "31 strong-colocalization pairs", pt=7.5)
text(x + 12, y0 + 166, "Extended MHC region excluded; effect", pt=7, italic=True, fill="#444444")
text(x + 12, y0 + 175, "scale converted to the per-allele scale", pt=7, italic=True, fill="#444444")

# ---------------- A2 ----------------
y0 = 256.0
heads = ["Phenotype", "Source", "Cases", "Controls"]
vals = ["MG (AChR+)", "European", "5,708", "432,028"]
xs = [x + 12, x + 62, x + 100, x + 126]
for i, (hh, vv) in enumerate(zip(heads, vals)):
    text(xs[i], y0 + 22, hh, pt=7.5, bold=True)
    text(xs[i], y0 + 33, vv, pt=7.5)
d.line([P(x + 11), P(y0 + 45), P(x + CW - 11), P(y0 + 45)], fill="#C9CED6", width=max(1, P(0.7)))
paste("a2_manh", x + 14, y0 + 62, 130, max_h=120, centre_x=x + CW / 2)
text(x + 12, y0 + 164, "European-ancestry meta-analysis", pt=7.5)
text(x + 12, y0 + 174, "5,708 cases / 432,028 controls", pt=7.5)

# ---------------- A3 ----------------
y0 = 468.0
abbr = [("cis-eQTL", "cis-expression QTL"), ("GWAS", "genome-wide association study"),
        ("LD", "linkage disequilibrium"), ("MG", "myasthenia gravis"),
        ("MR", "Mendelian randomization"), ("PheWAS", "phenome-wide assoc. study"),
        ("sc-eQTL", "single-cell eQTL"), ("UMAP", "manifold approximation")]
for i, (a, b) in enumerate(abbr):
    text(x + 12, y0 + 22 + i * 10.2, a, pt=7, bold=True)
    text(x + 52, y0 + 22 + i * 10.2, b, pt=7)
text(x + 12, y0 + 148, "Physical 1-Mb groups are descriptive,", pt=7, italic=True, fill="#444444")
text(x + 12, y0 + 157, "not LD-independent loci.", pt=7, italic=True, fill="#444444")
text(x + 12, y0 + 170, "BH, Benjamini-Hochberg; F, instrument", pt=7, italic=True, fill="#444444")
text(x + 12, y0 + 179, "strength statistic; PP.H4, posterior", pt=7, italic=True, fill="#444444")
text(x + 12, y0 + 188, "probability of a shared causal variant.", pt=7, italic=True, fill="#444444")

# ---------------- B1 ----------------
bx, y0 = 176.0, 44.0
rounded(bx + 14, y0 + 20, bx + CW - 14, y0 + 70, r=5, fill="#F7F7FA", outline="#B9B9C6")
text(bx + CW / 2, y0 + 25, "Confounders", pt=7.5, anchor="ma", fill="#333333")
for xx in (bx + 24, bx + CW - 33):
    d.line([P(xx), P(y0 + 38), P(xx + 8), P(y0 + 48)], fill="#C0392B", width=max(1, P(1.7)))
    d.line([P(xx + 8), P(y0 + 38), P(xx), P(y0 + 48)], fill="#C0392B", width=max(1, P(1.7)))
paste("b1_mr", bx + 14, y0 + 80, 130, max_h=46, centre_x=bx + CW / 2)
text(bx + 12, y0 + 138, "Single-instrument Wald ratios", pt=7.5)
text(bx + 12, y0 + 150, "Direction and rescaled magnitude only;", pt=7.5)
text(bx + 12, y0 + 160, "colocalization carries the shared-variant", pt=7.5)
text(bx + 12, y0 + 170, "interpretation", pt=7.5)

# ---------------- B2 ----------------
y0 = 256.0
half = (CW - 24 - 6) / 2
for i, (lab, col) in enumerate([("GWAS", "#2E6DA4"), ("eQTL", "#D9822B")]):
    x0c = bx + 12 + i * (half + 6)
    rounded(x0c, y0 + 20, x0c + half, y0 + 72, r=4, fill="#FBFBFD", outline="#DDDDE6")
    text(x0c + half / 2, y0 + 23, lab, pt=7.5, anchor="ma")
    for j in range(8):
        px = x0c + 7 + j * (half - 14) / 7
        ph = 30 * math.exp(-((j - 3.5) ** 2) / 3.2) + random.uniform(-2, 2)
        py = y0 + 66 - max(ph, 3) * 0.85
        d.ellipse([P(px - 1.4), P(py - 1.4), P(px + 1.4), P(py + 1.4)], fill=col)
d.line([P(bx + 14), P(y0 + 76), P(bx + CW - 14), P(y0 + 76)], fill="#9A9A9A", width=max(1, P(0.7)))
text(bx + CW / 2, y0 + 78, "Shared genomic region", pt=7, anchor="ma", fill="#444444")
hs = [("H0", "none"), ("H1", "GWAS"), ("H2", "eQTL"), ("H3", "shared"), ("H4", "distinct")]
bw = (CW - 24 - 4 * 3) / 5
for i, (h0, h1) in enumerate(hs):
    x0c = bx + 12 + i * (bw + 3)
    rounded(x0c, y0 + 90, x0c + bw, y0 + 134, r=3, fill="white", outline="#DDDDE6")
    text(x0c + bw / 2, y0 + 92, h0, pt=7, anchor="ma", bold=True)
    text(x0c + bw / 2, y0 + 101, h1, pt=7, anchor="ma", fill="#444444")
    for row, col in enumerate(["#2E6DA4", "#D9822B"]):
        yy = y0 + 112 + row * 12
        flat = (i == 1 and row == 1) or (i == 2 and row == 0) or (i == 0)
        if flat:
            d.line([P(x0c + 3), P(yy), P(x0c + bw - 3), P(yy)], fill="#C9CED6", width=max(1, P(0.7)))
        else:
            n = 15
            pts = [(P(x0c + 3 + t * (bw - 6) / (n - 1)),
                    P(yy - 5.0 * math.exp(-(((t - (n - 1) / 2.0) / 3.0) ** 2)) + 0.5)) for t in range(n)]
            d.line(pts, fill=col, width=max(1, P(1.2)))
text(bx + 12, y0 + 140, "31 strong pairs (PP.H4 > 0.80), 16 genes", pt=7.5)
text(bx + 12, y0 + 152, "Conditional shared-signal probability", pt=7.5)
text(bx + 12, y0 + 162, "PP.H4/(PP.H3+PP.H4) > 0.80 in all 31", pt=7.5)

# ---------------- B3 ----------------
y0 = 468.0
paste("b3_umap", bx + 12, y0 + 22, 96, max_h=112, centre_x=bx + CW / 2)
text(bx + 12, y0 + 140, "111,898 singlets from 7 donors", pt=7.5)
text(bx + 12, y0 + 150, "(3 HC, 4 MG); donor-level pseudobulk", pt=7.5)
text(bx + 12, y0 + 160, "MAIT-like cluster excluded; SESN3", pt=7.5)
text(bx + 12, y0 + 170, "lower in patients in both T-cell", pt=7.5)
text(bx + 12, y0 + 180, "populations", pt=7.5)

# ---------------- C1 ----------------
cx, y0 = 344.0, 44.0
paste("c1_phewas2", cx + 10, y0 + 20, 138, max_h=126, centre_x=cx + CW / 2)
text(cx + 12, y0 + 148, "12 of 28 lead variants with records;", pt=7.5)
text(cx + 12, y0 + 158, "44 distinct SNP-trait pairs", pt=7.5)
text(cx + 12, y0 + 170, "Descriptive coverage only; no", pt=7, italic=True, fill="#444444")
text(cx + 12, y0 + 179, "enrichment test was performed", pt=7, italic=True, fill="#444444")

# ---------------- C2 ----------------
y0 = 256.0
paste("c2_drug", cx + 12, y0 + 22, 116, max_h=110, centre_x=cx + CW / 2)
text(cx + 12, y0 + 140, "Annotations from DSigDB, Enrichr", pt=7.5)
text(cx + 12, y0 + 150, "and Open Targets", pt=7.5)
text(cx + 12, y0 + 162, "Hypothesis-generating only; not", pt=7, italic=True, fill="#444444")
text(cx + 12, y0 + 171, "evidence of efficacy in MG", pt=7, italic=True, fill="#444444")

# ---------------- C3 ----------------
y0 = 468.0
progs = [("c3_prog1", "T-cell regulatory and", "activation programs"),
         ("c3_prog2", "Myeloid, memory B and", "cytokine programs"),
         ("c3_prog3", "B-cell and plasma-cell", "regulatory programs"),
         ("c3_prog4", "Complement-related", "programs")]
yy = y0 + 22
for name, l1, l2 in progs:
    paste(name, cx + 12, yy, 15, max_h=15)
    text(cx + 32, yy + 0.5, l1, pt=7.5)
    text(cx + 32, yy + 10, l2, pt=7.5)
    yy += 23
text(cx + 12, yy + 2, "Prioritized genes", pt=7.5, bold=True)
gx = cx + 12
for g, col in [("SESN3", "#1F4E99"), ("CFD", "#375623")]:
    w = 32 if g == "CFD" else 44
    d.rounded_rectangle([P(gx), P(yy + 12), P(gx + w), P(yy + 24)], radius=P(6),
                        fill="white", outline=col, width=max(1, P(0.9)))
    text(gx + w / 2, yy + 14.5, g, pt=7.5, anchor="ma", fill=col, bold=True)
    gx += w + 6
text(cx + 12, yy + 28, "SESN1/3 axis in naive/central memory", pt=7, italic=True, fill="#444444")
text(cx + 12, yy + 37, "T cells; CFD at 19p13.3", pt=7, italic=True, fill="#444444")

png = os.path.join(OUT, "Figure1_workflow_v5.png")
img.save(png, optimize=True)
doc = fitz.open()
page = doc.new_page(width=W_PT, height=H_PT)
page.insert_image(fitz.Rect(0, 0, W_PT, H_PT), filename=png)
pdf = os.path.join(OUT, "Figure1_workflow_v5.pdf")
doc.save(pdf, deflate=True, garbage=4)
prev = r"D:\Deep harness\dump\fig1\v5_preview.png"
img.resize((int(W / 4.6), int(H / 4.6)), Image.LANCZOS).save(prev)
print("PNG %.0f KB  %dx%d px  == %.0f x %.0f mm at 600 dpi" % (os.path.getsize(png) / 1024, W, H, W_PT / 72 * 25.4, H_PT / 72 * 25.4))
print("PDF %.0f KB  page %.0f x %.0f pt" % (os.path.getsize(pdf) / 1024, W_PT, H_PT))
print("preview", prev)
