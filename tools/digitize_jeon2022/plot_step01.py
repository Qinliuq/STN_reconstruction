"""Figure for Steps 0-1: digitized Jeon 2022 profiles, fits, and the per-cell model scales they imply.

    python tools/digitize_jeon2022/plot_step01.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
os.chdir(ROOT)
from cfg import cfg  # noqa: E402
from tools import generate_locs  # noqa: E402
from stn_inputs import axis_position, profile_scale, gpe_cell_scale  # noqa: E402

DIG = os.path.join(ROOT, 'results', 'emergent_peak', 'digitized')
P = json.load(open(os.path.join(ROOT, 'params', 'jeon2022_profiles.json')))
COL = {'Gria1': '#c0392b', 'Gabra1': '#2e86c1', 'Pvalb': '#27ae60'}

pvp, pvn = generate_locs(cfg.pop_order, 0.2, seed=cfg.seeds['loc'], plot=False, pvp_placement='gradient')
ax = np.array([axis_position(l) for l in pvp + pvn])
isP = np.r_[np.ones(len(pvp), bool), np.zeros(len(pvn), bool)]
o = np.argsort(ax)

fig, axs = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
h = pd.read_csv(os.path.join(DIG, 'fig7H_fraction_positive.csv'))
a = axs[0, 0]
for g, s in h.groupby('gene'):
    a.fill_between(s.axis1, s.se_lower, s.se_upper, color=COL[g], alpha=0.12, lw=0)
    a.plot(s.axis1, s.fraction, color=COL[g], lw=1, label=f'{g} (digitized)')
    p = P['profiles'][g]
    a.plot(p['x'], p['y'], color=COL[g], ls='--', lw=1.5)
a.set(xlabel='STN axis 1', ylabel='fraction of positive cells', title='Fig. 7H: digitized (solid), cubic fit (dashed)')
a.legend(fontsize=7)

a = axs[0, 1]
I = pd.read_csv(os.path.join(DIG, 'fig7I_center_expression.csv'))
for (g, t), s in I.groupby(['gene', 'type']):
    a.plot(s.axis1, s.expression_au, color=COL[g], ls='-' if t == 'PV+' else ':', label=f'{g} {t}')
a.set(xlabel='STN axis 1', ylabel='expression (A.U.)', title='Fig. 7I center (PV+ solid, PV- dotted)')
a.legend(fontsize=7)

a = axs[0, 2]
e = pd.read_csv(os.path.join(DIG, 'fig1E_projection_profiles.csv'))
a.hist(e.sigma, bins=15, color='0.5')
a.axvline(P['gpe_topography']['topo_sigma'], color='k', label=f"median {P['gpe_topography']['topo_sigma']:.3f}")
a.axvline(0.1, color='tab:red', ls='--', label='config D 0.10')
a.set(xlabel='projection width sigma (axis units)', ylabel='injections', title='Fig. 1E: Gaussian width per injection')
a.legend(fontsize=7)

a = axs[1, 0]
a.plot(ax[o], profile_scale(ax, P['ctx_gscale'])[o], color=COL['Gria1'], label='Gria1 conductance scale (new)')
lin = lambda i, s: (i + s * ax) / (i + s * ax).mean()
a.plot(ax[o], lin(0.65, -0.5)[o], 'k--', label='config D rate scale 0.65 - 0.5a')
a.plot(ax[o], lin(0.55, -0.3)[o], 'k:', label='0.55 - 0.3a (cfg comment)')
a.set(xlabel='model axis', ylabel='per-cell factor (mean 1)', title='Glutamate')
a.legend(fontsize=7)

a = axs[1, 1]
for nm, spec, c in (('Gabra1 + ratio (F)', P['gpe_gscale'], 'tab:blue'), ('Gabra1, no ratio (abs)', P['gpe_gscale_abs'], 'tab:cyan'),
                    ('flat + ratio (noRec)', P['gpe_gscale_flat_with_ratio'], 'tab:gray')):
    sc = gpe_cell_scale(ax, isP, spec)
    a.scatter(ax[isP], sc[isP], marker='^', color=c, s=18)
    a.scatter(ax[~isP], sc[~isP], marker='o', color=c, s=10, label=nm)
cd = np.where(isP, 1.13, 0.87)
a.scatter(ax, cd / cd.mean(), marker='x', color='k', s=10, label='config D (1.13 / 0.87)')
a.set(xlabel='model axis', ylabel='GPe->STN factor (mean 1)', title='GABA per cell (triangles PV+)')
a.legend(fontsize=6)

a = axs[1, 2]
s1 = pd.read_csv(os.path.join(DIG, 'fig1H_gpe_vs_stn_axis1.csv'))
g = P['gpe_topography']
lo, hi = g['centroid_range']
u = (s1.gpe_input_box - s1.gpe_input_box.min()) / (s1.gpe_input_box.max() - s1.gpe_input_box.min())
c = lo + (s1.stn_projection_box - s1.stn_projection_box.min()) / (s1.stn_projection_box.max() - s1.stn_projection_box.min()) * (hi - lo)
a.scatter(u, c, s=12, color='#a93226')
xx = np.linspace(0, 1, 2)
a.plot(xx, g['group_center_map'][1] + g['group_center_map'][0] * xx, 'k-', label='fit: %.3f + %.3f u' % (g['group_center_map'][1], g['group_center_map'][0]))
a.plot(xx, 1 / 24 + xx * 22 / 24, 'r--', label='config D homes (even)')
a.set(xlabel='GPe input position (0..1)', ylabel='STN projection centroid (axis 1)', title='Fig. 1H + 1E: GPe -> STN mapping')
a.legend(fontsize=7)
fig.savefig(os.path.join(ROOT, 'results', 'emergent_peak', 'step01_profiles.png'), dpi=130)
print('saved')
