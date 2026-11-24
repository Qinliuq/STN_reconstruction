# Shared definitions for the data-derived PD drive topography runs (results/drive_topography; Priority 1 item 3).
# Imported by sweeps/drive_topography_c.py (variant c + controls + healthy references) and
# sweeps/drive_topography_ab.py (variants a and b). Configuration D otherwise unchanged; cfg.pd_drive is off by default.
#
# Seeds: INPUT seeds 1-5 (ctx_seed = s, gpe_seed = 100 + s) on one fixed network (cfg.seeds['loc'] = 1,
# ['conn'] = 4), as in Table 3. Peak positions are conditional on this one wiring.
# Drive: the PD CHANGE is distributed by s(x) (src/stn_drive_topography.py, mean 1):
#   GPe->STN synapse scale of cell i = 1 + (2.0 - 1) * s(x_i)                       (hooks 'both' / 'gpe')
#   GPe group beta depth = 0.1 + (M - 0.1) * s(c_g), connection-weighted mean = M    (hooks 'both' / 'beta')
#   with the other hook uniform (synapse 2.0 / depth M). M from params/drive_topography.json, one per comparison set
#   (advanced 0.8825, early 0.6757), the uniform control of each set included.
# Healthy references: the healthy state does not depend on the drive, so one healthy set per PV+ placement:
#   gradient -> sweeps/out/emergent_peak_followup/Dref_healthy, pvalb -> .../FU1_healthy (configuration D with the
#   Pvalb placement), uniform / monotonic -> H_<placement>_healthy here.
import json
import os
import sys

sys.path.insert(0, 'src')
import stn_drive_topography as D  # noqa: E402

P = json.load(open(os.path.join('params', 'drive_topography.json')))
duration = 3000
seeds = D.SEEDS
x_label = 'condition'
reference_group = 'healthy'
outdir = 'sweeps/out/drive_topography'

PD = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
HEALTHY = {'ctx_state': 'healthy', 'gpe_state': 'healthy'}
HOOKS = {'c': 'both', 'a': 'gpe', 'b': 'beta'}
CONTROLS = ('1_uniform_control', '5c_uniform_control_early')


def placement(name):
    return dict(P['placements'][name]['cfg'])


def pd_condition(name, variant):
    c = D.CONDITIONS[name]
    M = P['sets'][c['depth_set']]['M']
    drv = {'drive': dict(c['drive']), 'hooks': HOOKS[variant], 'depth_mean': M, 'depth_base': P['depth_base']}
    label = f'ctl_{name}_pd' if name in CONTROLS else f'{variant}_{name}_pd'
    return {'label': label, 'group': 'parkinsonian', 'x': list(D.CONDITIONS).index(name),
            'cfg': {**PD, **placement(c['pv_placement']), 'pd_drive': drv}}


def healthy_conditions():
    return [{'label': f'H_{pl}_healthy', 'group': 'healthy', 'x': -1, 'cfg': {**HEALTHY, **placement(pl)}}
            for pl in ('uniform', 'monotonic')]
