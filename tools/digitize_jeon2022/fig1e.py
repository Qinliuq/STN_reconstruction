"""Fig. 1E (left): per-injection STN projection profiles along STN axis 1 (heatmap rows), PDF page 3.

Method: colormap inversion on the native 300 ppi image. The colorbar under the heatmap (Low -> High) is sampled along
its middle rows into a lookup table (value 0..1 by position); each heatmap pixel gets the value of the nearest LUT
colour (RGB distance). Rows: 55 valid injections (Fig. 1B), equal height; each row's profile is the mean of its
central pixel rows. x is normalized 0..1 across the heatmap width (STN axis 1). Per row, fit
v(x) = b + A exp(-(x - mu)^2 / (2 sigma^2)) (least squares, bounds 0.02 <= sigma <= 1); report sigma, mu, fit R^2.
Caveat: the published profiles are density estimates (smoothed, normalized per row), so sigma includes that smoothing
and is an upper bound for the axonal spread.

    python tools/digitize_jeon2022/fig1e.py --src <dir with fig_p3-000.jpg>
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

import common as C

PAGE = 3
N_ROWS = 55


def bbox(im, x_probe, y_probe, x_range, y_range):
    s = im.sum(2)
    col = s[y_range[0]:y_range[1], x_probe] < 600
    ys = np.flatnonzero(col) + y_range[0]
    row = s[y_probe, x_range[0]:x_range[1]] < 700
    xs = np.flatnonzero(row) + x_range[0]
    return xs.min(), xs.max(), ys


def gauss(x, b, A, mu, sg):
    return b + A * np.exp(-(x - mu) ** 2 / (2 * sg ** 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    args = ap.parse_args()
    im = C.load(args.src, PAGE)
    s = im.sum(2)
    # heatmap edges: the outermost columns are dark blue (low) over the whole height
    med = np.median(s[470:770, :], axis=0)
    x0 = next(x for x in range(1255, 1300) if med[x] < 300 and med[x + 1] < 300 and med[x + 2] < 300)
    x1 = next(x for x in range(1465, 1400, -1) if med[x] < 600 and med[x - 1] < 600)
    colmask = s[440:810, x0 + 3] < 600
    runs_ = [(440 + a, 440 + b - 1) for a, b in C.runs(colmask)]
    heat = max(runs_, key=lambda r: r[1] - r[0])
    bar = [r for r in runs_ if r[0] > heat[1]][0]
    hy0, hy1 = heat
    by = (bar[0] + bar[1]) // 2
    lut = im[by - 1:by + 2, x0:x1 + 1].mean(0)  # colours Low -> High
    lutv = np.linspace(0, 1, len(lut))
    H = im[hy0:hy1 + 1, x0:x1 + 1].reshape(-1, 3).astype(float)
    d = ((H[:, None, :] - lut[None, :, :]) ** 2).sum(2)
    V = lutv[np.argmin(d, 1)].reshape(hy1 - hy0 + 1, x1 - x0 + 1)
    inv_err = np.sqrt(d.min(1)).reshape(V.shape)
    edges = np.linspace(0, V.shape[0], N_ROWS + 1)
    xs = (np.arange(V.shape[1]) + 0.5) / V.shape[1]
    rows, pts = [], []
    for k in range(N_ROWS):
        a, b = edges[k], edges[k + 1]
        c0, c1 = int(np.ceil(a + 0.25 * (b - a))), int(np.floor(b - 0.25 * (b - a)))
        prof = V[c0:c1 + 1].mean(0)
        try:
            p, _ = curve_fit(gauss, xs, prof, p0=[prof.min(), np.ptp(prof), xs[np.argmax(prof)], 0.15],
                             bounds=([0, 0, -0.2, 0.02], [1, 1.5, 1.2, 1.0]))
            r2 = 1 - ((prof - gauss(xs, *p)) ** 2).sum() / ((prof - prof.mean()) ** 2).sum()
        except RuntimeError:
            p, r2 = [np.nan] * 4, np.nan
        rows.append(dict(row=k, baseline=p[0], amplitude=p[1], mu=p[2], sigma=p[3], r2=r2,
                         peak_axis=xs[np.argmax(prof)], **{f'v{j:03d}': v for j, v in enumerate(prof)}))
        yc = hy0 + (a + b) / 2
        pts.append((x0 + p[2] * V.shape[1], yc, (0, 0, 0)))
        for sgn in (-1, 1):
            pts.append((x0 + (p[2] + sgn * p[3]) * V.shape[1], yc, (255, 0, 255)))
    df = pd.DataFrame(rows)
    os.makedirs(C.OUT, exist_ok=True)
    df.to_csv(os.path.join(C.OUT, 'fig1E_projection_profiles.csv'), index=False, float_format='%.4f')
    C.overlay(im, (x0 - 5, hy0 - 5, x1 + 6, bar[1] + 6), pts, os.path.join(C.OUT, 'overlay_fig1E.png'))
    good = df[df.r2 > 0.5]
    meta = dict(panel='Fig. 1E left (STN axis 1)', page=PAGE, heatmap_px=[int(x0), int(hy0), int(x1), int(hy1)],
                colorbar_row_px=int(by), n_rows=N_ROWS, row_height_px=float(V.shape[0] / N_ROWS),
                px_per_axis_unit=int(V.shape[1]), colormap_inversion_rms_rgb=float(inv_err.mean()),
                sigma_all=dict(median=float(df.sigma.median()), q25=float(df.sigma.quantile(.25)), q75=float(df.sigma.quantile(.75))),
                sigma_r2_gt_0_5=dict(n=int(len(good)), median=float(good.sigma.median()), q25=float(good.sigma.quantile(.25)),
                                     q75=float(good.sigma.quantile(.75))),
                fit_r2_median=float(df.r2.median()),
                reading_error='axis +-1 px = +-%.4f; value: LUT step %.4f; row assignment assumes 55 equal rows' % (1 / V.shape[1], 1 / len(lut)))
    json.dump(meta, open(os.path.join(C.OUT, 'fig1E_meta.json'), 'w'), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == '__main__':
    main()
