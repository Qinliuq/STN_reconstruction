"""Fig. 7I center (expression level along STN axis 1, PV- dotted / PV+ solid, Gria1 red / Gabra1 blue) and
Fig. 7I right (per-cell log2(Gria1:Gabra1) box plots, PV- vs PV+).

Center: no y ticks (A.U.); values are pixel heights above the x-axis line, then divided by the overall mean of the
four curves (shared arbitrary unit). Lines are traced by colour; the solid and dotted line of each colour are
separated by tracking two continuous paths from the left edge (solid starts low, dotted high) with linear
prediction, assigning each column's pixel runs to the nearer predicted path. Legend boxes are masked.
Right: y axis calibrated on its labelled ticks (-2..4); each box's median line (red), box edges (quartiles) and
whisker ends are read from the box outline pixels.

    python tools/digitize_jeon2022/fig7i.py --src <dir with fig_p11-000.jpg>
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

import common as C

PAGE = 11
PANEL = (1500, 1500, 1790, 1790)
LEGEND = (1568, 1495, 1790, 1580)
TESTS = {'Gria1': lambda r, g, b: (r > 150) & (g < 100) & (b < 100),
         'Gabra1': lambda r, g, b: (b > 200) & (r < 160) & (b - r > 60)}


def candidates(im, x, y0, y1, test):
    col = im[y0:y1, x]
    ok = test(col[:, 0], col[:, 1], col[:, 2])
    return [y0 + (a + b - 1) / 2 for a, b in C.runs(ok)]


def track_two(im, xs, y0, y1, test, jump=5):
    """Each colour has two lines (PV+ solid, PV- dotted) that cross once. Per column, the upper and lower line are the
    extreme candidate runs (a single run is given to the nearer of the two previous values). The crossing column is
    where the two are closest in the middle of the axis; before it the solid (PV+) line is the lower one, after it the
    upper one. Stray picks are removed against a running median; gaps are filled later."""
    up, lo = [], []
    for x in xs:
        cand = candidates(im, x, y0, y1, test)
        pu = next((v for v in reversed(up) if np.isfinite(v)), None)
        pl = next((v for v in reversed(lo) if np.isfinite(v)), None)
        if len(cand) >= 2:
            u, l = min(cand), max(cand)
        elif len(cand) == 1:
            c = cand[0]
            if pu is None or pl is None or abs(c - pu) <= abs(c - pl):
                u, l = c, np.nan
            else:
                u, l = np.nan, c
            if pu is not None and pl is not None and abs(pu - pl) < 4:  # lines touch: both at the run
                u, l = c, c
        else:
            u, l = np.nan, np.nan
        up.append(u)
        lo.append(l)
    up, lo = despike(np.array(up, float), jump), despike(np.array(lo, float), jump)
    n = len(xs)
    mid = np.arange(int(0.3 * n), int(0.65 * n))
    sep = C.fill_gaps(lo)[mid] - C.fill_gaps(up)[mid]
    xc = mid[int(np.nanargmin(sep))]
    solid = np.where(np.arange(n) < xc, lo, up)
    dotted = np.where(np.arange(n) < xc, up, lo)
    return {'solid': solid, 'dotted': dotted, 'crossing_index': xc}


def despike(y, tol, win=15):
    # drop points more than tol px from the running median (stray picks in bands or on the other line at crossings)
    out = y.copy()
    for i in range(len(y)):
        w = y[max(0, i - win):i + win + 1]
        w = w[np.isfinite(w)]
        if np.isfinite(y[i]) and len(w) >= 5 and abs(y[i] - np.median(w)) > tol:
            out[i] = np.nan
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    args = ap.parse_args()
    raw = C.load(args.src, PAGE)
    im = raw.copy()
    im[LEGEND[1]:LEGEND[3], LEGEND[0]:LEGEND[2]] = 255
    s = im.sum(2)
    axis_x = int(np.argmax([(s[1560:1760, x] < 300).sum() for x in range(1500, 1525)])) + 1500
    # x axis = the gradient arrow bar (7 px thick); baseline = centre row of the bar (the axis-label text lies below it)
    bar = [y for y in range(1740, 1776) if (s[y, axis_x + 5:axis_x + 200] < 600).sum() > 40]
    axis_y = int(round(np.mean(bar)))
    xs = np.arange(axis_x + 2, PANEL[2])
    out, pts, crossings = {}, [], {}
    for gene, test in TESTS.items():
        tr = track_two(im, xs, PANEL[1], axis_y - 1, test)
        crossings[gene] = float((xs[tr.pop('crossing_index')] - axis_x))
        for k, v in tr.items():
            out[(gene, k)] = C.fill_gaps(v, max_gap=40)
    last = max(int(np.flatnonzero(np.isfinite(v))[-1]) for v in out.values())
    x_end = xs[last]
    rows = []
    for (gene, style), v in out.items():
        for i, (x, y) in enumerate(zip(xs[:last + 1], v[:last + 1])):
            rows.append(dict(gene=gene, type='PV+' if style == 'solid' else 'PV-', axis1=(x - axis_x) / (x_end - axis_x),
                             height_px=axis_y - y if np.isfinite(y) else np.nan))
            if i % 3 == 0:
                pts.append((x, y, (255, 0, 255) if gene == 'Gria1' else (0, 0, 0)) if style == 'solid' else
                           (x, y, (255, 160, 0) if gene == 'Gria1' else (0, 160, 0)))
    df = pd.DataFrame(rows)
    df['expression_au'] = df.height_px / df.height_px.mean()
    os.makedirs(C.OUT, exist_ok=True)
    df.to_csv(os.path.join(C.OUT, 'fig7I_center_expression.csv'), index=False, float_format='%.4f')
    C.overlay(raw, PANEL, pts, os.path.join(C.OUT, 'overlay_fig7I_center.png'))
    meta = dict(panel='Fig. 7I center', page=PAGE, axis_x_px=axis_x, axis_y_px=axis_y, x_end_px=int(x_end),
                y_unit='A.U.: pixel height above the x axis / mean of all four curves (no y ticks in the panel)',
                axis_per_px=1 / (x_end - axis_x),
                crossing_axis1={g: v / (x_end - axis_x) for g, v in crossings.items()},
                missing_before_fill={f'{g}_{s}': int(np.isnan(v[:last + 1]).sum()) for (g, s), v in out.items()},
                reading_error='line centre +-1 px; baseline +-3 px (arrow bar thickness), i.e. ~2-6% of a curve height of 50-170 px; axis position +-2 px')
    json.dump(meta, open(os.path.join(C.OUT, 'fig7I_center_meta.json'), 'w'), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == '__main__':
    main()
