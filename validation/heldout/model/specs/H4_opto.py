# H4 held-out runs: virtual D2- / D1-MSN stimulation in the healthy state, compressed protocol (priors.json)
import json, os
_p = json.load(open(os.path.join('validation', 'heldout', 'priors.json'), encoding='utf-8'))['H4_optogenetics']
_proto = {k: _p['protocol'][k] for k in ('pre_s', 'on_s', 'post_s', 't_start_s')}
_sc = _p['scenarios']
duration = _p['duration_ms']
seeds = [1, 2, 3, 4, 5]
outdir = 'sweeps/out/heldout/H4_opto'
healthy = {'ctx_state': 'healthy', 'gpe_state': 'healthy'}
conditions = [
    {'label': 'D2_gpe', 'cfg': {**healthy, 'gpe_rate_mod': {'type': 'blocks', **_proto, 'factor': _sc['D2_gpe']['gpe_rate_mod_factor']}}},
    {'label': 'D2_ctx', 'cfg': {**healthy, 'ctx_rate_mod': {'type': 'blocks', **_proto, 'factor': _sc['D2_ctx']['ctx_rate_mod_factor']}}},
    {'label': 'D1_ctx', 'cfg': {**healthy, 'ctx_rate_mod': {'type': 'blocks', **_proto, 'factor': _sc['D1_ctx']['ctx_rate_mod_factor']}}},
]
