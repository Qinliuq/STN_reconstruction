# Parkinsonian network with uniform vs graded GPe beta depth, and with the PV+ gradient (PV+ concentrated around
# axis ~0.67) vs uniformly placed PV+ cells (tools.generate_locs pvp_placement; same cells per grid). Uniform beta
# depth: gpe_feedback_v3 beta_width very large, so every GPe group has the full depth (beta_modulation 1.0).
# Healthy and graded-beta PD with the PV+ gradient are the healthy / untreated runs of sweeps/kv31_combos.py.
duration = 3000
seeds = [1, 2, 3, 4, 5]
x_label = 'condition'
reference_group = 'healthy'
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
uni = {'gpe_params': {'beta_width': 1e9}}
conditions = [
    {'label': 'healthy', 'group': 'healthy', 'x': 0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
    {'label': 'untreated', 'group': 'PV+ gradient', 'x': 1, 'cfg': {**pd}},
    {'label': 'pd_unibeta', 'group': 'PV+ gradient', 'x': 2, 'cfg': {**pd, **uni}},
    {'label': 'healthy_pvpuni', 'group': 'healthy', 'x': 0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy', 'pvp_placement': 'uniform'}},
    {'label': 'pd_pvpuni', 'group': 'PV+ uniform', 'x': 1, 'cfg': {**pd, 'pvp_placement': 'uniform'}},
    {'label': 'pd_unibeta_pvpuni', 'group': 'PV+ uniform', 'x': 2, 'cfg': {**pd, **uni, 'pvp_placement': 'uniform'}},
]

# pinned to equal GPe->STN conductance in PV+ and PV- (default later moved to cfg.gpe_gaba_scale 1.13 / 0.87)
for c in conditions:
    c['cfg'] = {'gpe_gaba_scale': {'PVP': 1.0, 'PVN': 1.0}, **c['cfg']}
