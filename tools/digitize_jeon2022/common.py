"""Shared helpers for digitizing Jeon et al. 2022 (Cell Reports 38:110439) figure panels.

The figures are embedded in the PDF as 300 ppi JPEGs; they are extracted at native resolution with
`pdfimages -j -f <page> -l <page> <pdf> fig_p<page>` (no re-rendering). The PDF is not part of the repo: pass its path
or the extracted image directory with --src. Outputs go to results/emergent_peak/digitized/.
"""
import os

import numpy as np
from PIL import Image

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
OUT = os.path.join(ROOT, 'results', 'emergent_peak', 'digitized')


def load(src, page):
    p = os.path.join(src, f'fig_p{page}-000.jpg')
    if not os.path.exists(p):
        raise SystemExit(f'{p} missing: run  pdfimages -j -f {page} -l {page} <Jeon2022.pdf> {src}/fig_p{page}')
    return np.asarray(Image.open(p).convert('RGB')).astype(int)


def dark(im, thr=200):
    return im.sum(2) < thr


def runs(mask_1d):
    """(start, stop) of True runs."""
    m = np.concatenate([[False], mask_1d, [False]])
    d = np.diff(m.astype(int))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def line_centre(im, x, y0, y1, colour_test):
    """Mean row of pixels in column x (y0..y1) passing colour_test(r, g, b); nan if none."""
    col = im[y0:y1, x]
    ok = colour_test(col[:, 0], col[:, 1], col[:, 2])
    if not ok.any():
        return np.nan
    # if several separate runs, take the longest (a crossing line occludes the other)
    rr = runs(ok)
    a, b = max(rr, key=lambda r: r[1] - r[0])
    return y0 + (a + b - 1) / 2


def fill_gaps(y, max_gap=15):
    """Linear interpolation over gaps (occlusion by another line) up to max_gap columns; longer gaps stay nan."""
    y = np.array(y, float)
    good = np.flatnonzero(np.isfinite(y))
    out = y.copy()
    for a, b in zip(good[:-1], good[1:]):
        if 1 < b - a <= max_gap + 1:
            out[a:b + 1] = np.interp(np.arange(a, b + 1), [a, b], [y[a], y[b]])
    return out


def overlay(im, crop, points, path, scale=3):
    """Save a zoomed crop of the original panel with digitized points drawn on it (points: [(x, y, rgb)])."""
    from PIL import ImageDraw
    x0, y0, x1, y1 = crop
    img = Image.fromarray(im[y0:y1, x0:x1].astype(np.uint8)).resize(((x1 - x0) * scale, (y1 - y0) * scale), Image.NEAREST)
    d = ImageDraw.Draw(img)
    for x, y, c in points:
        if np.isfinite(y):
            X, Y = (x - x0 + 0.5) * scale, (y - y0 + 0.5) * scale
            d.ellipse((X - 2, Y - 2, X + 2, Y + 2), outline=c, fill=c)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    img.save(path)
