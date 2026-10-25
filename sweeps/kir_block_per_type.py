# Per-type Kir block on the topographic parkinsonian network with equal Kir in PV+ and PV- (0.3 each):
# Kir block in PV+ only, PV- only, or both (cfg.channel_block keys PVP_gKir / PVN_gKir / gKir).
# The healthy and untreated runs are the ctx0.8 runs of sweeps/kir_equal_calib.py (identical settings).
duration = 3000
seeds = [1, 2, 3]
x_label = 'Kir factor (1 = default, 0 = full block)'
reference_group = 'healthy'
baseline_label = 'untreated'
band_groups = ['healthy', 'untreated']
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
conditions = [{'label': 'healthy', 'group': 'healthy', 'x': 1.0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
              {'label': 'untreated', 'group': 'untreated', 'x': 1.0, 'cfg': {**pd}}]
for name, key in [('PV+ only', 'PVP_gKir'), ('PV- only', 'PVN_gKir'), ('both', 'gKir')]:
    for f in [0.5, 0.0]:
        conditions.append({'label': f"{key.replace('_gKir', '').lower() if key != 'gKir' else 'all'}_kir{f:g}", 'group': name, 'x': f,
                           'cfg': {**pd, 'channel_block': {key: f}}})

# pinned to the parkinsonian GPe->STN synapse scale these runs used (default later moved to 2.0)
for c in conditions:
    if c['cfg'].get('gpe_state') == 'parkinsonian':
        c['cfg']['gpe_params'] = {'synapse_scale': 1.5, **c['cfg'].get('gpe_params', {})}
