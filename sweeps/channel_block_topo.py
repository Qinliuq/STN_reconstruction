# HCN, CaT and Kir block, alone and combined with Kir block (same factor for both channels), in the
# topographic parkinsonian network (input_version 3, current cfg.py defaults; see results/topo_pd_state/).
# Same design as sweeps/channel_block_gpe_beta.py; the combinations are run at 0.5 and 0 only.
duration = 3000
seeds = [1, 2, 3]
x_label = 'conductance factor (1 = default, 0 = full block)'
reference_group = 'healthy'
baseline_label = 'untreated'
band_groups = ['healthy', 'untreated']
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
conditions = [{'label': 'healthy', 'group': 'healthy', 'x': 1.0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
              {'label': 'untreated', 'group': 'untreated', 'x': 1.0, 'cfg': {**pd}}]
blocks = {'Kir': (['gKir'], [0.75, 0.5, 0.25, 0.0]), 'HCN': (['gHCN'], [0.75, 0.5, 0.25, 0.0]),
          'CaT': (['gCaT'], [0.75, 0.5, 0.25, 0.0]), 'HCN+Kir': (['gHCN', 'gKir'], [0.5, 0.0]),
          'CaT+Kir': (['gCaT', 'gKir'], [0.5, 0.0])}
for name, (channels, factors) in blocks.items():
    for f in factors:
        conditions.append({'label': f"{name.lower().replace('+', '_')}{f:g}", 'group': name, 'x': f,
                           'cfg': {**pd, 'channel_block': {ch: f for ch in channels}}})

# pinned to the Kir scales and cortical conductance these runs used (Kir defaults later moved to 0.3 in both types)
for c in conditions:
    c['cfg'] = {'PVP_gKir_scale': 0.5, 'PVN_gKir_scale': 0.1, 'ctx_overrides': {'conductance_nS': 0.8}, **c['cfg']}
