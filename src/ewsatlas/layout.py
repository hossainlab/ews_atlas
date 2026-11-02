"""Figure layout: a 12-column grid at fixed width, panels rendered one by one
at their exact slot size, then tiled into a PNG and a vector PDF with bold
uppercase panel letters.

Geometry and type sizes come from the ``figures`` block of ``src/config.json``;
call :func:`setup` once before rendering.

An outline is ``{"tag", "ncol", "row_heights_mm", "panels": [{"letter", "row",
"col", "colspan"}]}``. Panels are drawn without titles (captions carry the
message) onto transparent canvases, so the composer can read each panel's
drawn extent from its alpha channel when placing letters.
"""

from __future__ import annotations

import os
import shutil
import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

from . import config

WIDTH_MM = DPI = GUTTER_MM = LETTER_PT = LETTER_GAP_MM = LETTER_BAND_MM = MM = None


def setup(cfg) -> None:
    """Load grid geometry from config and apply the shared matplotlib style."""
    global WIDTH_MM, DPI, GUTTER_MM, LETTER_PT, LETTER_GAP_MM, LETTER_BAND_MM, MM
    WIDTH_MM = float(config.get(cfg, "figures.width_mm"))
    DPI = int(config.get(cfg, "figures.dpi"))
    GUTTER_MM = float(config.get(cfg, "figures.gutter_mm"))
    LETTER_PT = float(config.get(cfg, "figures.letter_pt"))
    LETTER_GAP_MM = float(config.get(cfg, "figures.letter_gap_mm"))
    LETTER_BAND_MM = float(config.get(cfg, "figures.letter_band_mm"))
    MM = DPI / 25.4
    style()


def outline(tag, rows, panels, ncol=12):
    """Grid outline from (letter, row, col, colspan) tuples."""
    return {"tag": tag, "ncol": ncol, "row_heights_mm": list(rows),
            "panels": [{"letter": L, "row": r, "col": c, "colspan": s} for L, r, c, s in panels]}


def height_mm(outline) -> float:
    return grid(outline)[-1] / MM




def style():
    mpl.rcParams.update({
        "figure.dpi": DPI, "savefig.dpi": DPI, "font.size": 7,
        "axes.titlesize": 7, "axes.labelsize": 6.6, "legend.fontsize": 5.6,
        "xtick.labelsize": 5.6, "ytick.labelsize": 5.6,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.2, "ytick.major.size": 2.2,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "none", "savefig.transparent": True,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def _panel_px(outline, letter):
    """Pixel size of one panel's grid slot."""
    ncol = outline["ncol"]
    total_w = WIDTH_MM * MM
    gut = GUTTER_MM * MM
    colw = (total_w - gut * (ncol - 1)) / ncol
    p = next(x for x in outline["panels"] if x["letter"] == letter)
    w = colw * p["colspan"] + gut * (p["colspan"] - 1)
    h = outline["row_heights_mm"][p["row"]] * MM
    return int(round(w)), int(round(h))


def render(outline, letter, draw, outdir: Path, **kw) -> Path:
    """Render one panel at exactly its grid-slot size, transparent, no bbox trim."""
    w, h = _panel_px(outline, letter)
    fig = plt.figure(figsize=(w / DPI, h / DPI), dpi=DPI)
    ax = fig.add_subplot(111)
    draw(ax, **kw)
    fig.tight_layout(pad=0.45)
    p = outdir / f"{outline['tag']}_{letter}.png"
    fig.savefig(p, dpi=DPI, transparent=True)
    # same figure, same size, as vector: the PDF composition places these pages
    fig.savefig(p.with_suffix(".pdf"), dpi=DPI, transparent=True)
    plt.close(fig)
    return p


def grid(outline):
    """Pixel geometry shared by the PNG and PDF composers."""
    ncol = outline["ncol"]
    gut = int(round(GUTTER_MM * MM))
    total_w = int(round(WIDTH_MM * MM))
    colw = (total_w - gut * (ncol - 1)) / ncol
    # Reserve a band above each row for the panel letter. Stamping the letter
    # over the panel collides with whatever is at its top-left corner (a tall
    # rotated y-label, for instance); adding the band keeps every letter in
    # clean space without shrinking any panel.
    band = int(round(LETTER_BAND_MM * MM)) if len(outline["panels"]) > 1 else 0
    row_y, y = [], 0
    for hmm in outline["row_heights_mm"]:
        y += band
        row_y.append(y)
        y += int(round(hmm * MM)) + gut
    return gut, total_w, colw, band, row_y, y - gut


def _letter_baselines(outline, paths, row_y, band):
    """Baseline (px from page top) of each row's panel letters.

    Letters in one row share a baseline so they line up. It sits just above the
    highest drawn pixel in that row (read from the transparent panel PNGs), so a
    letter hugs its panel instead of floating at the top of an empty margin,
    and it never rises above the reserved band.
    """
    from PIL import Image

    cap = LETTER_PT / 72.0 * DPI * 0.73          # cap height of a bold sans
    gap = LETTER_GAP_MM * MM
    tops = {}
    for p in outline["panels"]:
        bb = Image.open(paths[p["letter"]]).getchannel("A").getbbox()
        t = row_y[p["row"]] + (bb[1] if bb else 0)
        tops[p["row"]] = min(tops.get(p["row"], t), t)
    return {r: max(row_y[r] - band + cap, tops[r] - gap) for r in tops}


def _letter_font(px):
    """Bold sans at px pixels. PIL cannot find fonts by family name on every
    platform (macOS has no DejaVuSans-Bold.ttf on its search path and PIL then
    falls back to a regular-weight bitmap face), so ask matplotlib, which ships
    DejaVu Sans Bold and matches the panels' own font."""
    from matplotlib import font_manager as fm
    from PIL import ImageFont

    path = fm.findfont(fm.FontProperties(family=mpl.rcParams["font.family"], weight="bold"),
                       fallback_to_default=True)
    return ImageFont.truetype(path, px)


def compose_pdf(outline, paths, out: Path) -> Path:
    """Vector twin of :func:`compose`: place each panel's PDF page on one page."""
    from pypdf import PageObject, PdfReader, PdfWriter, Transformation

    gut, total_w, colw, band, row_y, total_h = grid(outline)
    pt = 72.0 / DPI                                   # pixels -> PDF points
    W, H = total_w * pt, total_h * pt
    page = PageObject.create_blank_page(width=W, height=H)
    for p in outline["panels"]:
        src = PdfReader(str(paths[p["letter"]].with_suffix(".pdf"))).pages[0]
        x = round(p["col"] * (colw + gut)) * pt
        top = row_y[p["row"]] * pt
        y = H - top - float(src.mediabox.height)
        page.merge_transformed_page(src, Transformation().translate(x, y))
    if band:
        # letters as real text, drawn by matplotlib on a transparent full page
        base = _letter_baselines(outline, paths, row_y, band)
        lf = plt.figure(figsize=(W / 72.0, H / 72.0), dpi=72)
        for p in outline["panels"]:
            x = round(p["col"] * (colw + gut)) * pt
            yb = H - base[p["row"]] * pt
            lf.text(x / W, yb / H, p["letter"].upper(), ha="left", va="baseline",
                    fontsize=LETTER_PT, fontweight="bold", color="black")
        lp = paths[outline["panels"][0]["letter"]].with_name(f"{outline['tag']}_letters.pdf")
        lf.savefig(lp, transparent=True)
        plt.close(lf)
        page.merge_page(PdfReader(str(lp)).pages[0])
    w = PdfWriter()
    w.add_page(page)
    w.compress_identical_objects()
    tmp = out.with_name(out.stem + ".tmp.pdf")
    with open(tmp, "wb") as fh:
        w.write(fh)
    os.replace(tmp, out)
    return out


def compose(outline, paths, out: Path) -> Path:
    """Tile panels onto the grid and stamp bold panel letters."""
    from PIL import Image, ImageDraw

    gut, total_w, colw, band, row_y, total_h = grid(outline)
    canvas = Image.new("RGBA", (total_w, total_h), (255, 255, 255, 255))
    # ImageFont sizes are PIXELS, so the point size has to be scaled by the
    # output resolution: 14 px at 300 dpi is 3.4 pt, far below the ~8-9 pt
    # bold that journals ask for. Deriving it keeps the letters correct if
    # DPI or WIDTH_MM change.
    px = int(round(LETTER_PT / 72.0 * DPI))
    font = _letter_font(px)
    base = _letter_baselines(outline, paths, row_y, band) if band else {}
    drw = ImageDraw.Draw(canvas)
    for p in outline["panels"]:
        img = Image.open(paths[p["letter"]]).convert("RGBA")
        x = int(round(p["col"] * (colw + gut)))
        canvas.alpha_composite(img, (x, row_y[p["row"]]))
        if band:
            # Cell Press house style sets panel letters in bold UPPERCASE.
            drw.text((x, base[p["row"]]), p["letter"].upper(),
                     fill=(0, 0, 0, 255), font=font, anchor="ls")
    # Write via a temporary file: a transient lock on the target (a viewer, an
    # indexer) otherwise aborts the whole stage part way through the set.
    tmp = out.with_name(out.stem + ".tmp.png")
    canvas.convert("RGB").save(tmp, dpi=(DPI, DPI))
    # A viewer or indexer holding the previous PNG makes os.replace fail with
    # WinError 5; retry briefly rather than abandoning a whole figure set.
    for attempt in range(4):
        try:
            os.replace(tmp, out)
            break
        except PermissionError:
            if attempt < 3:
                time.sleep(0.5)
                continue
            # A file open in a viewer can often still be written in place even
            # though it cannot be renamed over (Windows share modes allow
            # read/write but deny delete). Fall back to an in-place copy so one
            # open figure does not abort the whole set.
            shutil.copyfile(tmp, out)
            tmp.unlink(missing_ok=True)
    return out
