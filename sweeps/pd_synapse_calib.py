# Restore the parkinsonian state with equal Kir (0.3 / 0.3) by strengthening the parkinsonian GPe->STN synapse
# (gpe_feedback_v3 STATES['parkinsonian']['synapse_scale'], 1.5 before). Healthy runs (synapse_scale 1.0) and
# the untreated PD at 1.5 are the ctx0.8 runs of sweeps/kir_equal_calib.py (identical settings).
duration = 3000
seeds = [1, 2, 3]
x_label = 'parkinsonian GPe->STN synapse scale'
reference_group = 'healthy'
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
conditions = [{'label': 'healthy', 'group': 'healthy', 'x': 1.0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}}]
conditions += [{'label': f'pd_syn{s:g}', 'group': 'parkinsonian', 'x': s, 'cfg': {**pd, 'gpe_params': {'synapse_scale': s}}}
               for s in [1.5, 1.75, 2.0, 2.5]]
