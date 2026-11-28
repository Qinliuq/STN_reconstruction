"""Check, before any run, that the advanced-set depth scale M (0.8825, params/drive_topography.json) keeps every GPe
group depth <= 1 under the fixes (F1 positions, F3 48 groups); and report realised positions under F1.
    python tools/edge_artefacts/check_depth.py
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
sys.path.insert(0, os.path.join(ROOT, 'tools', 'drive_topography'))
import gpe_feedback_v3 as G  # noqa: E402
import stn_drive_topography as D  # noqa: E402
import precompute as PC  # noqa: E402
from tools import generate_locs  # noqa: E402

P = json.load(open(os.path.join(ROOT, 'params', 'drive_topography.json')))
M = P['sets']['advanced']['M']
F3 = {'n_gpe_groups': 48, 'home_range': [1 / 24, 23 / 24]}


def positions(pl, grid):
    pc = PC.placement_cfg(pl)
    pvp, pvn = generate_locs(60, 0.2, seed=PC.LOC_SEED, plot=False, pvp_placement=pc.get('pvp_placement', 'gradient'),
                             pvp_profile=pc.get('pvp_profile'), grid=grid)
    return np.clip([(c['xnorm'] + c['ynorm']) / 2 for c in pvp + pvn], 0, 1), len(pvp)


if __name__ == '__main__':
    out = {}
    for grid in ('original', 'centred'):
        for pl in ('gradient', 'uniform', 'monotonic', 'pvalb'):
            ax, npv = positions(pl, grid)
            print(f'{grid:8s} {pl:9s} range {ax.min():.3f}-{ax.max():.3f}  exact0 {(ax == 0).sum()}  PV+ {np.round(np.sort(ax[:npv]), 2).tolist()}')
            out[f'{grid}/{pl}'] = dict(range=[float(ax.min()), float(ax.max())], pvp=np.sort(ax[:npv]).round(3).tolist())
    for grid, gp in (('original', {}), ('centred', {}), ('centred', F3), ('original', F3)):
        for cond in ('2_territory_adv', '2b_territory_adv_pvalb', '3_pv_uniform', '6_pv_monotonic'):
            c = D.CONDITIONS[cond]
            ax, _ = positions(c['pv_placement'], grid)
            mx = []
            for s in D.SEEDS:
                g = G.GPeInputGenerator(n_stn=len(ax), dt_s=0.001, stn_axis=ax, seed=100 + s, feedback_delay_s=0.004,
                                        beta_depth_shape=lambda x, kw=c['drive']: D.drive_scaling(np.asarray(x), **kw),
                                        beta_depth_mean=M, beta_depth_base=0.1, **{**G.STATES['parkinsonian'], **gp})
                mx.append(g.beta_depth.max())
            print(f'{grid:8s} groups {gp.get("n_gpe_groups", 12):2d} {cond:24s} max depth {max(mx):.4f}')
    json.dump(out, open(os.path.join(ROOT, 'results', 'edge_artefacts', 'f1_positions.json'), 'w'), indent=1)
