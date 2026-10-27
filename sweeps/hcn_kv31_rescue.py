# HCN block and Kv3.1 enhancement (opener-like: factor > 1 on the Kv3.1 scale in all STN cells) in the topographic
# parkinsonian network with equal Kir (0.3 / 0.3) and parkinsonian GPe->STN synapse scale 2.0 (current defaults).
# The healthy and untreated runs are the healthy / pd_syn2 runs of sweeps/pd_synapse_calib.py (identical settings).
duration = 3000
seeds = [1, 2, 3]
x_label = 'conductance factor (1 = default)'
reference_group = 'healthy'
baseline_label = 'untreated'
band_groups = ['healthy', 'untreated']
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
conditions = [{'label': 'healthy', 'group': 'healthy', 'x': 1.0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
              {'label': 'untreated', 'group': 'untreated', 'x': 1.0, 'cfg': {**pd}}]
conditions += [{'label': f'hcn{f:g}', 'group': 'HCN block', 'x': f, 'cfg': {**pd, 'channel_block': {'gHCN': f}}}
               for f in [0.75, 0.5, 0.25, 0.0]]
conditions += [{'label': f'kv31x{f:g}', 'group': 'Kv3.1 increase', 'x': f, 'cfg': {**pd, 'channel_block': {'gKv31': f}}}
               for f in [1.1, 1.25, 1.5, 2.0]]

# pinned to equal GPe->STN conductance in PV+ and PV- (default later moved to cfg.gpe_gaba_scale 1.13 / 0.87)
for c in conditions:
    c['cfg'] = {'gpe_gaba_scale': {'PVP': 1.0, 'PVN': 1.0}, **c['cfg']}
