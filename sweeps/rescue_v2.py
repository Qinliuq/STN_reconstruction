# Channel block, Kv3.1 enhancement and combinations on the current default model: configuration D (graded GPe beta,
# topographic GPe<->STN), equal Kir 0.3, Kv3.1 PV+ 1.2 / PV- 1.15, parkinsonian GPe->STN synapse scale 2.0 and the
# per-type GABA:glutamate ratio (gpe_gaba_scale PV+ 1.13 / PV- 0.87). All factors act on all STN cells
# (cfg.channel_block); Kv3.1 factors > 1 model a Kv3.1 opener. The healthy and untreated runs are those of
# sweeps/gaba_ratio.py (seeds 1-5).
duration = 3000
seeds = [1, 2, 3, 4, 5]
x_label = 'factor (1 = default)'
reference_group = 'healthy'
baseline_label = 'untreated'
band_groups = ['healthy', 'untreated']
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
conditions = [{'label': 'healthy', 'group': 'healthy', 'x': 1.0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
              {'label': 'untreated', 'group': 'untreated', 'x': 1.0, 'cfg': {**pd}}]
singles = {'CaT': ('gCaT', [0.9, 0.8, 0.5]), 'Kir': ('gKir', [0.8, 0.5, 0.0]), 'HCN': ('gHCN', [0.75, 0.5, 0.0]),
           'Kv3.1': ('gKv31', [1.1, 1.2, 1.25, 1.3, 1.5, 2.0])}
for name, (key, factors) in singles.items():
    conditions += [{'label': f"{name.lower().replace('.', '')}{'x' if f > 1 else ''}{f:g}", 'group': name, 'x': f,
                    'cfg': {**pd, 'channel_block': {key: f}}} for f in factors]
combos = [('gCaT', 0.9), ('gCaT', 0.8), ('gKir', 0.9), ('gKir', 0.8), ('gHCN', 0.5)]
conditions += [{'label': f"{key[1:].lower()}{f:g}_kv31x1.25", 'group': 'with Kv3.1 x1.25', 'x': f,
                'cfg': {**pd, 'channel_block': {key: f, 'gKv31': 1.25}}} for key, f in combos]
# stage 2: the same combinations at Kv3.1 x1.3 (the most consistent Kv3.1 dose of stage 1)
conditions += [{'label': f"{key[1:].lower()}{f:g}_kv31x1.3", 'group': 'with Kv3.1 x1.3', 'x': f,
                'cfg': {**pd, 'channel_block': {key: f, 'gKv31': 1.3}}} for key, f in combos]
