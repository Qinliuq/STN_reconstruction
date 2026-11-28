# Cortical->STN topography test (results/ctx_topography): variants T0-T3 (src/ctx_topography.py, cfg.ctx_topo_variant)
# x pallidal beta depth {uniform: gpe_params beta_width 1e9 (as Dref_uni_pd), configuration D graded (default)}, with
# bias fixes F1 (cfg.loc_grid 'centred') + F2 (cfg.gpe_coverage_norm). Configuration D otherwise; input seeds 1-5
# (ctx_seed s, gpe_seed 100 + s) on one fixed network (loc 1, conn 4).
# Healthy references: one per variant. The healthy GPe depth is 0.1 in every group whatever beta_width is
# (gpe_feedback_v3.STATES['healthy']: beta_floor = beta_modulation = 0.1), so the healthy state does not depend on the
# pallidal-depth factor and each variant's healthy set is the matched reference for both PD depths.
# CTX_TOPO_SEEDS=1 restricts every condition to seed 1 (operating-point smoke check before the full set).
import os

duration = 3000
seeds = [int(s) for s in os.environ.get('CTX_TOPO_SEEDS', '1,2,3,4,5').split(',')]
x_label = 'condition'
reference_group = 'healthy'
outdir = 'sweeps/out/ctx_topography'

FIX = {'loc_grid': 'centred', 'gpe_coverage_norm': True}
PD = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
HEALTHY = {'ctx_state': 'healthy', 'gpe_state': 'healthy'}
DEPTH = {'graded': {}, 'uniform': {'gpe_params': {'beta_width': 1e9}}}
VARIANTS = ['T0', 'T1', 'T2', 'T3']

conditions = []
for i, v in enumerate(VARIANTS):
    conditions.append({'label': f'{v}_healthy', 'group': 'healthy', 'x': i,
                       'cfg': {**HEALTHY, **FIX, 'ctx_topo_variant': v}})
    for dep, dcfg in DEPTH.items():
        conditions.append({'label': f'{v}_{dep}_pd', 'group': 'parkinsonian', 'x': i,
                           'cfg': {**PD, **FIX, **dcfg, 'ctx_topo_variant': v}})
