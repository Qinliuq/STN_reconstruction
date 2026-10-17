# Which input carries the parkinsonian beta (input_version 2)? Healthy vs parkinsonian, and each source alone:
# pd_gpe_only = parkinsonian GPe + healthy cortex, pd_ctx_only = healthy GPe + parkinsonian cortex.
duration = 3000
seeds = [1, 2, 3]
x_label = 'condition'
conditions = [
    {'label': 'healthy', 'group': 'healthy', 'x': 0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
    {'label': 'pd_ctx_only', 'group': 'parkinsonian', 'x': 1, 'cfg': {'ctx_state': 'parkinsonian', 'gpe_state': 'healthy'}},
    {'label': 'pd_gpe_only', 'group': 'parkinsonian', 'x': 2, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'parkinsonian'}},
    {'label': 'pd', 'group': 'parkinsonian', 'x': 3, 'cfg': {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}},
]
