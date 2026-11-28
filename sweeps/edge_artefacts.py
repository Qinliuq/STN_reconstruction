# Edge / coverage artefact check (results/edge_artefacts): bias fixes F1-F3 on configuration D, variant b (beta depth
# only), input seeds 1-5 on the fixed network. Matched healthy references per fix set and per PV+ placement.
#   F1 cfg.loc_grid = 'centred'; F2 cfg.gpe_coverage_norm = True;
#   F3 cfg.gpe_params = {'n_gpe_groups': 48, 'home_range': [1/24, 23/24]} + cfg.gpe_drive_ref = {'n_groups': 12}
# R0 = existing drive_topography runs (reused) + its geometry control here; R1 F1; R2 F2; R3 F1+F2; R4 F1+F2+F3.
# Geometry control 'geo': uniform beta depth (condition 1's drive) + uniform STN PV+ placement.
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drive_topography_common import *  # noqa: F401,F403,E402
from drive_topography_common import D, HEALTHY, P, PD, placement, pd_condition  # noqa: E402

outdir = 'sweeps/out/edge_artefacts'
F1 = {'loc_grid': 'centred'}
F2 = {'gpe_coverage_norm': True}
F3 = {'gpe_params': {'n_gpe_groups': 48, 'home_range': [1 / 24, 23 / 24]}, 'gpe_drive_ref': {'n_groups': 12}}
RUNS = {'R1': ({**F1}, ['1_uniform_control', '2_territory_adv', '3_pv_uniform', '6_pv_monotonic', '2b_territory_adv_pvalb', 'geo']),
        'R2': ({**F2}, ['1_uniform_control', '2_territory_adv', '3_pv_uniform', '6_pv_monotonic', '2b_territory_adv_pvalb']),
        'R3': ({**F1, **F2}, ['1_uniform_control', '2_territory_adv', '3_pv_uniform', '6_pv_monotonic', '2b_territory_adv_pvalb', 'geo']),
        # R4 (F1+F2+F3) HELD: its seed-1 smoke run (sweeps/edge_artefacts_check.py) shifted PV+ / PV- rates by -17% / -30%
        # and the burst fraction by -53% vs configuration D (limits 10% / 25%); not recalibrated, see results/edge_artefacts
        # 'R4': ({**F1, **F2, **F3}, ['1_uniform_control', '2_territory_adv', '3_pv_uniform']),
        'R0': ({}, ['geo'])}
# R1 also gets the geometry control (cheap; the request asks for R0 and R3, R1 isolates F1's effect on it)


def pd_run(run, cond):
    fix = RUNS[run][0]
    if cond == 'geo':
        c = pd_condition('1_uniform_control', 'b')
        c['cfg'] = {**c['cfg'], **placement('uniform')}
    else:
        c = pd_condition(cond, 'b')
    return {'label': f'{run}_{cond}_pd', 'group': 'parkinsonian', 'x': 0, 'cfg': {**c['cfg'], **fix}}


def placement_of(cond):
    return 'uniform' if cond == 'geo' else D.CONDITIONS[cond]['pv_placement']


conditions = []
for run, (fix, conds) in RUNS.items():
    conditions += [pd_run(run, c) for c in conds]
    if run != 'R0':  # R0 healthy references exist (drive_topography / emergent_peak_followup)
        for pl in sorted({placement_of(c) for c in conds}):
            conditions.append({'label': f'{run}_H_{pl}_healthy', 'group': 'healthy', 'x': -1,
                               'cfg': {**HEALTHY, **placement(pl), **fix}})
