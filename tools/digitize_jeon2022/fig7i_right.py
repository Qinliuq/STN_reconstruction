"""Fig. 7I right: per-cell log2(Gria1:Gabra1) box plots for PV- and PV+ (PDF page 11, native image).

y axis: tick marks right of the y-axis line, labelled 4, 3, ..., -2 from the top (linear fit). Each box: median =
mean row of the red median line inside the box; quartiles = the box's top / bottom black edges, read in a column
offset from the whisker; whisker ends = extent of the black vertical whisker line through the box centre.

    python tools/digitize_jeon2022/fig7i_right.py --src <dir with fig_p11-000.jpg>
"""
import argparse
import json
import os

import numpy as np

import common as C

PAGE = 11
PANEL = (1880, 1490, 2038, 1780)
BOXES = {'PV-': 1930, 'PV+': 1990}  # approximate centre columns (refined below)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    args = ap.parse_args()
    im = C.load(args.src, PAGE)
    s = im.sum(2)
    sat = im.max(2) - im.min(2)
    y0, y1 = PANEL[1], PANEL[3]
    axis_x = int(np.argmax([(s[y0:y1, x] < 200).sum() for x in range(1885, 1900)])) + 1885
    rows = [y for y in range(y0, y1) if (s[y, axis_x + 2:axis_x + 6] < 300).all() and sat[y, axis_x + 2:axis_x + 6].mean() < 40]
    groups = np.split(np.array(rows), np.flatnonzero(np.diff(rows) > 2) + 1)
    ticks = np.array([g.mean() for g in groups if len(g)])[:7]  # labelled 4 ... -2; the 8th dark row is the x-axis line
    vals = 4 - np.arange(len(ticks))
    a, b = np.polyfit(ticks, vals, 1)
    resid = float(np.abs(vals - (a * ticks + b)).max())
    to_val = lambda y: a * y + b
    out, pts = {}, []
    for name, cx0 in BOXES.items():
        # whisker column: the darkest column near the guessed centre
        cx = int(np.argmax([(s[y0:y1, x] < 150).sum() for x in range(cx0 - 8, cx0 + 9)])) + cx0 - 8
        red = (im[..., 0] > 180) & (im[..., 1] < 90) & (im[..., 2] < 90)
        ry, rx = np.nonzero(red[y0:y1, cx - 25:cx + 26])
        med = y0 + ry.mean()
        col = s[y0:y1, cx + 12] < 150
        dark_runs = [(y0 + p, y0 + q) for p, q in C.runs(col)]
        near = [r for r in dark_runs if abs((r[0] + r[1]) / 2 - med) < 60]
        q3_row = (near[0][0] + near[0][1] - 1) / 2   # top edge
        q1_row = (near[-1][0] + near[-1][1] - 1) / 2  # bottom edge
        wcol = s[y0:y1, cx - 1:cx + 2].min(1) < 150  # whisker line (1-3 px wide) around the centre column
        rr = [(y0 + p, y0 + q - 1) for p, q in C.runs(wcol)]
        above = [r for r in rr if r[0] < q3_row - 2 <= r[1]]
        below = [r for r in rr if r[0] <= q1_row + 2 < r[1]]
        w_top = above[0][0] if above else np.nan
        w_bot = below[0][1] if below else np.nan
        out[name] = dict(median=to_val(med), q1=to_val(q1_row), q3=to_val(q3_row), whisker_low=to_val(w_bot),
                         whisker_high=to_val(w_top), centre_px=cx,
                         rows_px=dict(median=med, q1=q1_row, q3=q3_row, whisker_low=w_bot, whisker_high=w_top))
        for y, c in ((med, (255, 0, 255)), (q1_row, (0, 160, 255)), (q3_row, (0, 160, 255)), (w_bot, (255, 140, 0)), (w_top, (255, 140, 0))):
            for dx in (-6, 0, 6):
                pts.append((cx + dx, y, c))
    meta = dict(panel='Fig. 7I right', page=PAGE, axis_x_px=axis_x, ticks_px=ticks.round(2).tolist(), tick_values=vals.tolist(),
                per_px=abs(a), tick_fit_max_residual=resid, boxes=out,
                reading_error=f'+-1 px = +-{abs(a):.3f} log2 units')
    os.makedirs(C.OUT, exist_ok=True)
    json.dump(meta, open(os.path.join(C.OUT, 'fig7I_right_boxes.json'), 'w'), indent=1, default=float)
    C.overlay(im, PANEL, pts, os.path.join(C.OUT, 'overlay_fig7I_right.png'))
    print(json.dumps(meta, indent=1, default=float))


if __name__ == '__main__':
    main()
