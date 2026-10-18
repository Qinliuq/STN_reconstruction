# quick check of input_version 3 at the default network size (pop_order 60)
duration = 600
seeds = [1]
conditions = [
    {'label': 'healthy', 'group': 'healthy', 'x': 0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
    {'label': 'pd', 'group': 'pd', 'x': 1, 'cfg': {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}},
]
