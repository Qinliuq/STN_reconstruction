# H3 held-out run: parkinsonian configuration D + anesthesia-like cortical slow oscillation (priors.json)
import json, os
_p = json.load(open(os.path.join('validation', 'heldout', 'priors.json')))['H3_slow_oscillation']
duration = _p['duration_ms']
seeds = [1, 2, 3, 4, 5]
outdir = 'sweeps/out/heldout/H3_slowosc'
conditions = [{'label': 'pd_so', 'cfg': {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian', 'ctx_rate_mod': _p['ctx_rate_mod']}}]
