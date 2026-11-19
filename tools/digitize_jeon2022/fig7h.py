"""Fig. 7H (middle): fraction of positive cells along STN axis 1 for Gria1, Gabra1, Pvalb (mean +- SE bands).

Method: colour masks on the native 300 ppi image (PDF page 11). Per pixel column, the line centre is the middle of the
longest run of line-coloured pixels; gaps where another line occludes it are interpolated (<= 15 columns). The y axis
is calibrated on the tick marks (detected right of the y-axis line; labelled 0.6 and 0.3, spacing 0.1). The x axis
maps the y-axis line (axis 1 = 0) to the last data column (axis 1 = 1); the panel has no x ticks, so this assumes the
curves span the full normalized axis, as the arrow does.
SE band: from each line centre, the contiguous run of tinted (non-white) pixels up and down; columns where the run
reaches another curve's line are flagged as overlapping and excluded from the band-width summary.

    python tools/digitize_jeon2022/fig7h.py --src <dir with fig_p11-000.jpg>
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

import common as C

PAGE = 11
PANEL = (440, 1500, 760, 1775)  # x0, y0, x1, y1 in the native image
LEGEND = (638, 1495, 760, 1577)  # legend box (coloured dashes and names), masked before line tracing
TESTS = {
    'Gria1': lambda r, g, b: (r > 150) & (g < 90) & (b < 90),
    'Gabra1': lambda r, g, b: (b > 200) & (r < 150) & (b - r > 70),
    'Pvalb': lambda r, g, b: (g > 110) & (r < 110) & (b < 90) & (g - r > 40),
}
DRAW = {'Gria1': (255, 0, 255), 'Gabra1': (0, 0, 0), 'Pvalb': (255, 140, 0)}


def calibrate_y(im, axis_x):
    s = im.sum(2)
    seg = im[:, axis_x + 2:axis_x + 7]
    grey = (seg.max(2) - seg.min(2)).mean(1) < 40
    rows = [y for y in range(PANEL[1], 1760) if s[y, axis_x + 2:axis_x + 7].mean() < 400 and grey[y]]
    ticks = [np.mean(g) for g in np.split(np.array(rows), np.flatnonzero(np.diff(rows) > 2) + 1) if len(g)]
    ticks = np.array(ticks)
    vals = 0.6 - 0.1 * np.arange(len(ticks))  # first tick is labelled 0.6, then every 0.1
    a, b = np.polyfit(ticks, vals, 1)
    resid = vals - (a * ticks + b)
    return a, b, ticks, vals, float(np.abs(resid).max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    args = ap.parse_args()
    im = C.load(args.src, PAGE)
    raw = im.copy()
    im = im.copy()
    im[LEGEND[1]:LEGEND[3], LEGEND[0]:LEGEND[2]] = 255
    s = im.sum(2)
    axis_x = int(np.argmax([(s[PANEL[1]:1760, x] < 300).sum() for x in range(450, 470)])) + 450
    a, b, ticks, vals, tick_err = calibrate_y(im, axis_x)
    to_val = lambda y: a * y + b
    y0, y1 = PANEL[1], 1760
    cols = np.arange(axis_x + 1, PANEL[2])
    centres = {k: np.array([C.line_centre(im, x, y0, y1, t) for x in cols]) for k, t in TESTS.items()}
    last = max(int(np.flatnonzero(np.isfinite(centres[k]))[-1]) for k in ('Gria1', 'Pvalb'))
    x_end = cols[last]
    cols, centres = cols[:last + 1], {k: v[:last + 1] for k, v in centres.items()}
    to_axis = lambda x: (x - axis_x) / (x_end - axis_x)
    rows, points, bands = [], [], {}
    white = s > 735
    for k, yc in centres.items():
        yf = C.fill_gaps(yc)
        lo, hi, ov = [], [], []
        for x, y in zip(cols, yf):
            if not np.isfinite(y):
                lo.append(np.nan); hi.append(np.nan); ov.append(True)
                continue
            yi = int(round(y))
            up = yi
            while up > y0 and not white[up - 1, x]:
                up -= 1
            dn = yi
            while dn < y1 - 1 and not white[dn + 1, x]:
                dn += 1
            other = any(np.isfinite(centres[o][x - cols[0]]) and up <= centres[o][x - cols[0]] <= dn for o in centres if o != k)
            lo.append(dn); hi.append(up); ov.append(other)
        bands[k] = dict(lo=np.array(lo, float), hi=np.array(hi, float), overlap=np.array(ov))
        for i, (x, y) in enumerate(zip(cols, yf)):
            rows.append(dict(gene=k, axis1=to_axis(x), fraction=to_val(y) if np.isfinite(y) else np.nan,
                             se_lower=to_val(lo[i]) if np.isfinite(lo[i]) else np.nan,
                             se_upper=to_val(hi[i]) if np.isfinite(hi[i]) else np.nan,
                             band_overlaps_other=ov[i], interpolated=not np.isfinite(yc[i]) and np.isfinite(y)))
            if i % 3 == 0:
                points.append((x, y, DRAW[k]))
    df = pd.DataFrame(rows)
    os.makedirs(C.OUT, exist_ok=True)
    df.to_csv(os.path.join(C.OUT, 'fig7H_fraction_positive.csv'), index=False, float_format='%.4f')
    C.overlay(raw, PANEL, points, os.path.join(C.OUT, 'overlay_fig7H.png'))
    px_y, px_x = abs(a), 1 / (x_end - axis_x)
    meta = dict(panel='Fig. 7H middle', page=PAGE, image='native 300 ppi JPEG (pdfimages)', axis_x_px=int(axis_x),
                x_end_px=int(x_end), y_ticks_px=ticks.round(2).tolist(), y_tick_values=vals.round(2).tolist(),
                y_per_px=px_y, axis_per_px=px_x, tick_fit_max_residual=tick_err,
                se_halfwidth_median={k: float(np.nanmedian((to_val(v['hi']) - to_val(v['lo']))[~v['overlap']] / 2))
                                     for k, v in bands.items()},
                interpolated_columns={k: int(df[(df.gene == k)].interpolated.sum()) for k in centres},
                reading_error=f'line centre +-1 px = +-{px_y:.4f} fraction; axis position +-2 px = +-{2 * px_x:.4f}')
    json.dump(meta, open(os.path.join(C.OUT, 'fig7H_meta.json'), 'w'), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == '__main__':
    main()
