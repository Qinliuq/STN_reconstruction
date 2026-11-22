"""Fig. 3C (PDF page 7): hyperdirect projection (HP_STN) along STN axis 1 per cortical injection, with the cortical
group colour stripe (green sensorimotor, orange associative, red limbic), and HP_GPe along the GPe gradient axis 1.

Method: colormap inversion against the Fig. 3B colorbar (same colormap; 3C has none of its own), per pixel row of the
heatmaps (rows are ~1.8 px high; each pixel row is treated as one sample and the group is read from the stripe pixel
of that row). Per row the centroid along the axis is the Gaussian-fit mean (as Fig. 1E; fallback: weighted mean of the
values above the row minimum). Somatomotor locus = median centroid of the sensorimotor (green) rows, uncertainty =
bootstrap 95% CI of that median over rows and the IQR (spread of the territory). The same for HP_GPe, mapped to the
STN with the Fig. 1H mapping in params (cross-check).

    python tools/digitize_jeon2022/fig3c.py --src <dir with fig_p7-000.jpg>
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

import common as C

PAGE = 7
BAR = (858, 1000, 380)  # x0, x1, row of the Fig. 3B colorbar
STRIPE_X = (1547, 1554)
GROUPS = {'sensorimotor': lambda r, g, b: (g > r + 40) & (g > 100),
          'limbic': lambda r, g, b: (r > 140) & (g < 70) & (b < 80),
          'associative': lambda r, g, b: (r > 200) & (g > 110) & (g < 210) & (b < 90)}


def gauss(x, b, A, mu, sg):
    return b + A * np.exp(-(x - mu) ** 2 / (2 * sg ** 2))


def heat_edges(s, lo, hi):
    med = np.median(s[60:360, :], axis=0)
    xs = [x for x in range(lo, hi) if med[x] < 700]
    return min(xs), max(xs)


def centroids(im, x0, x1, y0, y1, lut, lutv):
    H = im[y0:y1 + 1, x0:x1 + 1].reshape(-1, 3).astype(float)
    d = ((H[:, None, :] - lut[None]) ** 2).sum(2)
    V = lutv[np.argmin(d, 1)].reshape(y1 - y0 + 1, x1 - x0 + 1)
    V[:, :2] = V[:, 2:3]  # first columns carry anti-aliased stripe colour (orange/red read as 'high'); copy column 3
    xs = (np.arange(V.shape[1]) + 0.5) / V.shape[1]
    out = []
    for prof in V:
        try:
            p, _ = curve_fit(gauss, xs, prof, p0=[prof.min(), np.ptp(prof) + 1e-3, xs[np.argmax(prof)], 0.15],
                             bounds=([0, 0, -0.2, 0.02], [1, 1.5, 1.2, 1.0]))
            mu, sg = p[2], p[3]
        except (RuntimeError, ValueError):
            w = prof - prof.min()
            mu, sg = (xs * w).sum() / w.sum() if w.sum() > 0 else np.nan, np.nan
        out.append((mu, sg))
    return np.array(out), float(np.sqrt(d.min(1)).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    args = ap.parse_args()
    im = C.load(args.src, PAGE)
    s = im.sum(2)
    bx0, bx1, by = BAR
    row = s[by, bx0 - 10:bx1 + 10] < 700
    bx0, bx1 = bx0 - 10 + int(np.flatnonzero(row)[0]), bx0 - 10 + int(np.flatnonzero(row)[-1])
    lut = im[by - 2:by + 3, bx0:bx1 + 1].mean(0)
    lutv = np.linspace(0, 1, len(lut))
    hx0, hx1 = heat_edges(s, STRIPE_X[1] + 1, 1700)
    col = s[30:400, (hx0 + hx1) // 2] < 760
    rr = max(C.runs(col), key=lambda r: r[1] - r[0])
    y0, y1 = 30 + rr[0], 30 + rr[1] - 1
    gx0, gx1 = heat_edges(s, 1880, 2038)
    stn, err1 = centroids(im, hx0, hx1, y0, y1, lut, lutv)
    gpe, err2 = centroids(im, gx0, gx1, y0, y1, lut, lutv)
    stripe = im[y0:y1 + 1, STRIPE_X[0] + 2:STRIPE_X[1] - 1].mean(1)
    grp = np.array(['unknown'] * len(stripe), dtype=object)
    for g, t in GROUPS.items():
        grp[t(stripe[:, 0], stripe[:, 1], stripe[:, 2])] = g
    df = pd.DataFrame(dict(y_px=np.arange(y0, y1 + 1), group=grp, stn_centroid=stn[:, 0], stn_sigma=stn[:, 1],
                           gpe_centroid=gpe[:, 0]))
    os.makedirs(C.OUT, exist_ok=True)
    df.to_csv(os.path.join(C.OUT, 'fig3C_hyperdirect_rows.csv'), index=False, float_format='%.4f')
    P = json.load(open(os.path.join(C.ROOT, 'params', 'jeon2022_profiles.json')))
    slope, intercept = P['gpe_topography']['group_center_map']
    rng = np.random.default_rng(0)
    summ = {}
    for g in ('sensorimotor', 'associative', 'limbic'):
        sub = df[df.group == g].dropna(subset=['stn_centroid'])
        v = sub.stn_centroid.values
        boot = np.median(rng.choice(v, (5000, len(v))), axis=1)
        gp = sub.gpe_centroid.dropna().values
        summ[g] = dict(n_pixel_rows=int(len(v)), stn_median=float(np.median(v)), stn_ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                       stn_iqr=[float(np.percentile(v, 25)), float(np.percentile(v, 75))],
                       gpe_median=float(np.median(gp)), stn_via_gpe_mapping=float(intercept + slope * np.median(gp)))
    pts = []
    for _, r in df.iterrows():
        if np.isfinite(r.stn_centroid) and r.y_px % 3 == 0:
            c = {'sensorimotor': (0, 0, 0), 'associative': (255, 140, 0), 'limbic': (255, 0, 255)}.get(r.group, (128, 128, 128))
            pts.append((hx0 + r.stn_centroid * (hx1 - hx0 + 1), r.y_px, c))
    C.overlay(im, (STRIPE_X[0] - 4, y0 - 4, hx1 + 5, y1 + 5), pts, os.path.join(C.OUT, 'overlay_fig3C.png'))
    meta = dict(panel='Fig. 3C', page=PAGE, colorbar_px=[bx0, bx1, by], stn_heatmap_px=[hx0, y0, hx1, y1], gpe_heatmap_px=[gx0, y0, gx1, y1],
                colormap_inversion_rms_rgb=[err1, err2], n_pixel_rows=int(len(df)),
                group_counts={g: int((df.group == g).sum()) for g in list(GROUPS) + ['unknown']}, groups=summ,
                reading_error='axis +-1 px = +-%.4f; rows are ~1.8 px high, so neighbouring pixel rows are not independent '
                              '(CI of the median is optimistic by ~sqrt(1.8))' % (1 / (hx1 - hx0 + 1)))
    json.dump(meta, open(os.path.join(C.OUT, 'fig3C_meta.json'), 'w'), indent=1, default=float)
    print(json.dumps(meta, indent=1, default=float))


if __name__ == '__main__':
    main()
