# Tuning so that GPe is the primary beta source (input_version 2): GPe conductance x GABA decay,
# healthy vs parkinsonian GPe only vs parkinsonian cortex only (2 s, seed 1)
duration = 2000
seeds = [1]
x_label = 'condition (0 healthy, 1 pd_ctx_only, 2 pd_gpe_only)'
conditions = []
for g in [1.0, 2.0]:
    for tau2 in [15.0, 8.0]:
        net = {'gpe_conductance_nS': g, 'gpe_syn': {'tau2': tau2}}
        tag = f'g{g:g}_tau{tau2:g}'
        conditions += [
            {'label': f'{tag}_healthy', 'group': tag, 'x': 0, 'cfg': {**net, 'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
            {'label': f'{tag}_ctx', 'group': tag, 'x': 1, 'cfg': {**net, 'ctx_state': 'parkinsonian', 'gpe_state': 'healthy'}},
            {'label': f'{tag}_gpe', 'group': tag, 'x': 2, 'cfg': {**net, 'ctx_state': 'healthy', 'gpe_state': 'parkinsonian'}},
        ]

# pinned to the model version these runs used (defaults later moved to input_version 3, pop_order 60)
for c in conditions:
    c['cfg'] = {'input_version': 2, 'pop_order': 30, 'gpe_conductance_nS': 0.5, 'ctx_overrides': {}, **c['cfg']}
