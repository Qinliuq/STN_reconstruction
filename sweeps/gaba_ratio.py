# Per-type GABA:glutamate ratio (cfg.gpe_gaba_scale PV+ 1.13 / PV- 0.87, default) in configuration D (GPe beta
# graded along the axis, topographic GPe<->STN, 1 nS base): healthy vs parkinsonian, and the uniform-beta-depth
# parkinsonian network with the PV+ gradient. Compare with the same conditions without the ratio in
# sweeps/beta_pvp_layout.py (healthy, untreated, pd_unibeta; same seeds).
duration = 3000
seeds = [1, 2, 3, 4, 5]
x_label = 'condition'
reference_group = 'healthy'
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
conditions = [
    {'label': 'healthy', 'group': 'healthy', 'x': 0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
    {'label': 'pd', 'group': 'parkinsonian', 'x': 1, 'cfg': {**pd}},
    {'label': 'pd_unibeta', 'group': 'parkinsonian', 'x': 2, 'cfg': {**pd, 'gpe_params': {'beta_width': 1e9}}},
]
