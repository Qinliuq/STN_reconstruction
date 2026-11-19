"""Fig. 1H right (top): GPe input position vs STN projection centroid along STN axis 1 (PDF page 3).

The panel has no tick labels: both axes are read in box-normalized units (0 at the axis line, 1 at the axis end).
Points are open red circles; centres are found by matched filtering of the red mask with a ring template (radius
searched 3.5-5.5 px), then greedy non-maximum suppression (min. distance 3 px), keeping peaks above a fraction of the
best response. Overlapping circles can merge, so the count is compared with the 55 valid injections of Fig. 1B.
Fit: STN = a + b * GPe (OLS), residual SD, Spearman rho (paper: 0.886).

    python tools/digitize_jeon2022/fig1h.py --src <dir with fig_p3-000.jpg>
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import fftconvolve

import common as C

PAGE = 3
TEXT = (1836, 1350, 1955, 1375)  # x0, y0, x1, y1 of the in-panel label, masked


def ring(r, w=1.3):
    n = int(np.ceil(r + 3))
    yy, xx = np.mgrid[-n:n + 1, -n:n + 1]
    d = np.sqrt(xx ** 2 + yy ** 2)
    t = np.exp(-(d - r) ** 2 / (2 * (w / 2) ** 2))
    return t - t.mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    ap.add_argument('--keep', type=float, default=0.45, help='peak threshold as a fraction of the best response')
    args = ap.parse_args()
    im = C.load(args.src, PAGE)
    raw = im.copy()
    s = im.sum(2)
    ax_x = int(np.argmax([(s[1330:1540, x] < 200).sum() for x in range(1825, 1840)])) + 1825
    ax_y = int(np.argmax([(s[y, 1830:2030] < 200).sum() for y in range(1525, 1550)])) + 1525
    col = s[1300:ax_y, ax_x] < 200
    top = 1300 + int(np.flatnonzero(col)[0])
    row = s[ax_y, ax_x:2038] < 200
    right = ax_x + int(np.flatnonzero(row)[-1])
    x0, y0, x1, y1 = ax_x + 2, top, right, ax_y - 2
    im = im.copy()
    im[TEXT[1]:TEXT[3], TEXT[0]:TEXT[2]] = 255  # red 'STN axis1' label inside the plot box (top left)
    sub = im[y0:y1, x0:x1]
    red = ((sub[..., 0] > 120) & (sub[..., 1] < 110) & (sub[..., 2] < 120)).astype(float)
    best_r, best = None, None
    for r in np.arange(3.5, 5.6, 0.25):
        resp = fftconvolve(red, ring(r)[::-1, ::-1], mode='same')
        if best is None or resp.max() > best.max():
            best_r, best = r, resp
    thr = args.keep * best.max()
    cand = np.argwhere(best > thr)
    order = np.argsort(-best[cand[:, 0], cand[:, 1]])
    peaks = []
    for i in order:
        y, x = cand[i]
        if all((y - py) ** 2 + (x - px) ** 2 >= 9 for py, px in peaks):
            peaks.append((y, x))
    P = np.array(peaks, float)
    gx = (P[:, 1] + x0 - ax_x) / (right - ax_x)
    gy = (ax_y - (P[:, 0] + y0)) / (ax_y - top)
    df = pd.DataFrame(dict(gpe_input_box=gx, stn_projection_box=gy, x_px=P[:, 1] + x0, y_px=P[:, 0] + y0))
    fit = stats.linregress(gx, gy)
    resid_sd = float(np.std(gy - (fit.intercept + fit.slope * gx), ddof=2))
    rho = stats.spearmanr(gx, gy).correlation
    os.makedirs(C.OUT, exist_ok=True)
    df.to_csv(os.path.join(C.OUT, 'fig1H_gpe_vs_stn_axis1.csv'), index=False, float_format='%.4f')
    C.overlay(raw, (ax_x - 5, top - 5, right + 5, ax_y + 5), [(x, y, (0, 160, 255)) for x, y in zip(df.x_px, df.y_px)],
              os.path.join(C.OUT, 'overlay_fig1H.png'))
    meta = dict(panel='Fig. 1H right top', page=PAGE, axis_box_px=[ax_x, top, right, ax_y], ring_radius_px=float(best_r),
                n_points=int(len(df)), n_expected=55, spearman_rho=float(rho), paper_rho=0.886,
                ols_box_units=dict(slope=float(fit.slope), intercept=float(fit.intercept), r=float(fit.rvalue),
                                   slope_se=float(fit.stderr), resid_sd=resid_sd),
                note='no tick labels: units are fractions of the axis box; data range in box units: GPe %.2f-%.2f, STN %.2f-%.2f'
                     % (gx.min(), gx.max(), gy.min(), gy.max()),
                reading_error='circle centre +-1 px = +-%.4f (x), +-%.4f (y) box units' % (1 / (right - ax_x), 1 / (ax_y - top)))
    json.dump(meta, open(os.path.join(C.OUT, 'fig1H_meta.json'), 'w'), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == '__main__':
    main()
