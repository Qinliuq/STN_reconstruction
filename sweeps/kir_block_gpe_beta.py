# Kir block in the parkinsonian network with input_version 2 (GPe is the primary beta source, cortical beta a
# supplement). Healthy reference: same network, healthy GPe and cortical input (not swept).
duration = 3000
seeds = [1, 2, 3]
x_label = 'Kir factor (1 = default, 0 = full block)'
reference_group = 'healthy'
baseline_label = 'pd_kir1'
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
conditions = [{'label': 'healthy', 'group': 'healthy', 'x': 1.0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}}]
conditions += [{'label': f'pd_kir{f:g}', 'group': 'parkinsonian', 'x': f, 'cfg': {**pd, 'channel_block': {'gKir': f}}}
               for f in [1.0, 0.75, 0.5, 0.25, 0.0]]

# pinned to the model version these runs used (defaults later moved to input_version 3, pop_order 60)
for c in conditions:
    c['cfg'] = {'input_version': 2, 'pop_order': 30, 'gpe_conductance_nS': 0.5, 'ctx_overrides': {}, **c['cfg']}
