"""Fig. S5D (PDF page 33, 220 ppi image): human 7T-MRI hyperdirect profile along STN axis 1 for the motor territory
(green mean line). Supplementary only: human tractography, an axis defined on the human STN, and a bimodal profile,
so it is reported next to the mouse locus (Fig. 3C), not used to set it.

Method: per pixel column, the dark-green mean line (thick) centre; height above the panel baseline; centroid and
peak positions over the curve's x span (0..1).

    python tools/digitize_jeon2022/figs5d.py --src <dir with fig_p33-000.ppm>
"""
import argparse
import json
import os

import numpy as np
from PIL import Image

import common as C

PANEL = (1040, 1015, 1125, 1058)  # motor sub-panel: x0, y0, x1, baseline row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    args = ap.parse_args()
    im = np.asarray(Image.open(os.path.join(args.src, 'fig_p33-000.ppm')).convert('RGB')).astype(int)
    x0, y0, x1, base = PANEL
    s = im.sum(2)
    base = int(np.argmax([(s[y, x0:x1] < 300).sum() for y in range(base - 4, base + 6)])) + base - 4
    test = lambda r, g, b: (g > 120) & (r < 110) & (b < 110) & (g - r > 50)
    xs = np.arange(x0, x1)
    yc = np.array([C.line_centre(im, x, y0, base, test) for x in xs])
    ok = np.isfinite(yc)
    xs, h = xs[ok], base - yc[ok]
    a = (xs - xs[0]) / (xs[-1] - xs[0])
    cen = float((a * h).sum() / h.sum())
    peaks = [float(a[i]) for i in range(1, len(h) - 1) if h[i] >= h[i - 1] and h[i] >= h[i + 1] and h[i] > 0.6 * h.max()]
    C.overlay(im, (x0 - 3, y0 - 3, x1 + 3, base + 3), [(x, base - v, (255, 0, 255)) for x, v in zip(xs, h)],
              os.path.join(C.OUT, 'overlay_figS5D_motor.png'), scale=6)
    meta = dict(panel='Fig. S5D motor (human)', page=33, baseline_row=int(base), x_span_px=[int(xs[0]), int(xs[-1])],
                centroid_axis1=cen, peaks_axis1_gt60pct=sorted(set(round(p, 2) for p in peaks)),
                reading_error='+-1 px = +-%.3f axis (85 px panel); y has no ticks' % (1 / (xs[-1] - xs[0])),
                note='supplementary: human, bimodal; not used for the model locus')
    json.dump(meta, open(os.path.join(C.OUT, 'figS5D_meta.json'), 'w'), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == '__main__':
    main()
