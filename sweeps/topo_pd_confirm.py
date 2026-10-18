# Confirmation of the best calibration setting (topo_sigma 0.1, base GPe->STN 1 nS, parkinsonian synapse
# scale 1.5 = gpe_feedback_v3 default) with more seeds, and a lower cortical conductance to bring rates down.
duration = 3000
seeds = [1, 2, 3]
x_label = 'condition'
reference_group = 'healthy'
base = {'gpe_conductance_nS': 1.0, 'gpe_params': {'topo_sigma': 0.10}}
conditions = []
for ctx in [1.0, 0.7]:
    net = {**base, 'ctx_overrides': {'conductance_nS': ctx}}
    conditions += [
        {'label': f'ctx{ctx:g}_healthy', 'group': f'ctx{ctx:g}', 'x': 0, 'cfg': {**net, 'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
        {'label': f'ctx{ctx:g}_pd', 'group': f'ctx{ctx:g}', 'x': 1, 'cfg': {**net, 'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}},
    ]
