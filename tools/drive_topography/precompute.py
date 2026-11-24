"""Depth scale M per comparison set for the data-derived drive topography (results/drive_topography).

Committed before any run. For every condition of src/stn_drive_topography.CONDITIONS and every input seed
(gpe_seed = 100 + s), rebuild the realised STN positions (tools.generate_locs with the fixed network seed
cfg.seeds['loc'] and the condition's PV+ placement, in the network's cell order PV+ then PV-) and the realised GPe->STN
connections (gpe_feedback_v3.GPeInputGenerator with the parkinsonian state, i.e. the same random draws as the runs),
and read the normalized group shape s_g (connection-weighted mean over cells = 1). Per comparison set:
    M = base + (1 - base) / max(s_g over all conditions and seeds of the set), base = 0.1 (healthy depth),
rounded DOWN to 4 decimals, so that every group depth base + (M - base) * s_g <= 1. One M per set, control included.
Also reports the per-cell synapse-scale range 1 + (2 - 1) * s(x_i) and the PV+ positions per placement.

    python tools/drive_topography/precompute.py   -> params/drive_topography.json
"""
import datetime
import json
import math
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
import gpe_feedback_v3 as G  # noqa: E402
import stn_drive_topography as D  # noqa: E402
from tools import generate_locs  # noqa: E402

BASE = 0.1
SEEDS = D.SEEDS
LOC_SEED = 1  # cfg.seeds['loc']
GRID = np.round(np.arange(0, 1.0001, 0.05), 2)


def placement_cfg(name):
    if name == 'gradient':
        return {}
    if name == 'uniform':
        return {'pvp_placement': 'uniform'}
    if name == 'monotonic':
        return {'pvp_placement': 'profile', 'pvp_profile': {'x': GRID.tolist(), 'y': np.round(D.pv_shape_monotonic(GRID), 4).tolist()}}
    if name == 'pvalb':
        P = json.load(open(os.path.join(ROOT, 'params', 'jeon2022_profiles.json')))
        return {'pvp_placement': 'profile', 'pvp_profile': P['pvp_profile']}
    raise ValueError(name)


def positions(name):
    pc = placement_cfg(name)
    pvp, pvn = generate_locs(60, 0.2, seed=LOC_SEED, plot=False, pvp_placement=pc.get('pvp_placement', 'gradient'),
                             pvp_profile=pc.get('pvp_profile'))
    ax = np.clip([(c['xnorm'] + c['ynorm']) / 2 for c in pvp + pvn], 0, 1)
    return ax, len(pvp)


def main():
    import cfg as _  # noqa: F401  (only to fail early if src/cfg.py is broken)
    from cfg import cfg
    st = dict(G.STATES['parkinsonian'])
    out = dict(created=datetime.date.today().isoformat(),
               git_at_creation=subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True).stdout.strip(),
               depth_base=BASE, input_seeds=SEEDS, loc_seed=LOC_SEED, placements={}, conditions={}, sets={})
    for pl in ('gradient', 'uniform', 'monotonic', 'pvalb'):
        ax, npv = positions(pl)
        out['placements'][pl] = dict(cfg=placement_cfg(pl), n_pvp=int(npv), n_cells=int(len(ax)),
                                     pvp_axis=sorted(np.round(ax[:npv], 3).tolist()), axis_range=[float(ax.min()), float(ax.max())])
    for name, c in D.CONDITIONS.items():
        ax, _ = positions(c['pv_placement'])
        s_cell = D.drive_scaling(ax, **c['drive'])
        s_cell = s_cell / s_cell.mean()
        smax = []
        for seed in SEEDS:
            gen = G.GPeInputGenerator(n_stn=len(ax), dt_s=cfg.gpe_dt, stn_axis=ax, seed=100 + seed,
                                      feedback_delay_s=max(0.0, cfg.gpe_feedback_delay - cfg.gpe_dt),
                                      beta_depth_shape=lambda x, kw=c['drive']: D.drive_scaling(np.asarray(x), **kw),
                                      beta_depth_mean=BASE, beta_depth_base=BASE, **st)
            smax.append(float(gen.beta_depth_shape_norm.max()))
        out['conditions'][name] = dict(depth_set=c['depth_set'], pv_placement=c['pv_placement'], drive=c['drive'],
                                       s_cell_range=[float(s_cell.min()), float(s_cell.max())],
                                       synapse_scale_range=[float(1 + s_cell.min()), float(1 + s_cell.max())],
                                       group_s_max_per_seed=smax)
    for setname in sorted({c['depth_set'] for c in D.CONDITIONS.values()}):
        members = [n for n, c in out['conditions'].items() if c['depth_set'] == setname]
        smax = max(max(out['conditions'][n]['group_s_max_per_seed']) for n in members)
        M = math.floor((BASE + (1 - BASE) / smax) * 1e4) / 1e4
        out['sets'][setname] = dict(members=members, binding_max_s=smax, M=M,
                                    binding_condition=max(members, key=lambda n: max(out['conditions'][n]['group_s_max_per_seed'])),
                                    total_beta_drive_vs_configD_uniform=M)
    os.makedirs(os.path.join(ROOT, 'params'), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, 'params', 'drive_topography.json'), 'w'), indent=1)
    for k, v in out['sets'].items():
        print(k, v)
    for n, c in out['conditions'].items():
        print(f"{n:26s} s_cell {c['s_cell_range'][0]:.3f}-{c['s_cell_range'][1]:.3f}  synapse {c['synapse_scale_range'][0]:.3f}-"
              f"{c['synapse_scale_range'][1]:.3f}  group s max {min(c['group_s_max_per_seed']):.3f}-{max(c['group_s_max_per_seed']):.3f}")
    for pl, p in out['placements'].items():
        print(pl, p['n_pvp'], p['n_cells'], p['pvp_axis'])


if __name__ == '__main__':
    main()
