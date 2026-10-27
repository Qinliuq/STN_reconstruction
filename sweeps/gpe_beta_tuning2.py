# Round 2: all GPe groups project to every STN cell (target_probability 1.0), parkinsonian beta depth 1.0,
# GABA decay 8 ms; GPe conductance 0.5 vs 1.0 nS. Healthy vs parkinsonian GPe only vs parkinsonian cortex only.
duration = 2000
seeds = [1, 2]
x_label = 'condition (0 healthy, 1 pd_ctx_only, 2 pd_gpe_only)'
conditions = []
for g in [0.5, 1.0]:
    net = {'gpe_conductance_nS': g, 'gpe_syn': {'tau2': 8.0}, 'gpe_params': {'target_probability': 1.0}}
    pdGpe = {'gpe_params': {'target_probability': 1.0, 'beta_modulation': 1.0}}
    tag = f'g{g:g}'
    conditions += [
        {'label': f'{tag}_healthy', 'group': tag, 'x': 0, 'cfg': {**net, 'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
        {'label': f'{tag}_ctx', 'group': tag, 'x': 1, 'cfg': {**net, 'ctx_state': 'parkinsonian', 'gpe_state': 'healthy'}},
        {'label': f'{tag}_gpe', 'group': tag, 'x': 2, 'cfg': {**net, **pdGpe, 'ctx_state': 'healthy', 'gpe_state': 'parkinsonian'}},
    ]

# pinned to the model version these runs used (defaults later moved to input_version 3, pop_order 60, Kv3.1 1.2 / 1.15)
for c in conditions:
    c['cfg'] = {'input_version': 2, 'pop_order': 30, 'gpe_conductance_nS': 0.5, 'ctx_overrides': {},
                'PVP_gKv31_scale': 1.0, 'PVN_gKv31_scale': 1.0, **c['cfg']}

# pinned to the Kir scales these runs used (defaults later moved to 0.3 in both types)
for c in conditions:
    c['cfg'] = {'PVP_gKir_scale': 0.5, 'PVN_gKir_scale': 0.1, **c['cfg']}

# pinned to equal GPe->STN conductance in PV+ and PV- (default later moved to cfg.gpe_gaba_scale 1.13 / 0.87)
for c in conditions:
    c['cfg'] = {'gpe_gaba_scale': {'PVP': 1.0, 'PVN': 1.0}, **c['cfg']}
