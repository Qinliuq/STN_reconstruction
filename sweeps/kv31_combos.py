# Follow-up to sweeps/hcn_kv31_rescue.py on the same network (current defaults):
# - Kv3.1 x1.25 and x1.5 with 8 seeds (healthy and untreated also 8 seeds);
# - Kv3.1 x1.25 combined with HCN 0.5, with partial CaT block (0.9, 0.8) or with slight Kir block (0.9, 0.8; all
#   STN cells, equal Kir), with the single drugs as references (5 seeds each).
# Seeds 1-3 of healthy / untreated / hcn0.5 / kv31x1.25 / kv31x1.5 are the runs of sweeps/hcn_kv31_rescue.py.
duration = 3000
seeds = [1, 2, 3, 4, 5, 6, 7, 8]
few = [1, 2, 3, 4, 5]
x_label = 'condition'
reference_group = 'healthy'
baseline_label = 'untreated'
band_groups = ['healthy', 'untreated']
pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
conditions = [{'label': 'healthy', 'group': 'healthy', 'x': 0, 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}},
              {'label': 'untreated', 'group': 'untreated', 'x': 0, 'cfg': {**pd}},
              {'label': 'kv31x1.25', 'group': 'Kv3.1', 'x': 1.25, 'cfg': {**pd, 'channel_block': {'gKv31': 1.25}}},
              {'label': 'kv31x1.5', 'group': 'Kv3.1', 'x': 1.5, 'cfg': {**pd, 'channel_block': {'gKv31': 1.5}}},
              {'label': 'hcn0.5', 'group': 'single', 'x': 1, 'seeds': few, 'cfg': {**pd, 'channel_block': {'gHCN': 0.5}}},
              {'label': 'cat0.9', 'group': 'single', 'x': 2, 'seeds': few, 'cfg': {**pd, 'channel_block': {'gCaT': 0.9}}},
              {'label': 'cat0.8', 'group': 'single', 'x': 3, 'seeds': few, 'cfg': {**pd, 'channel_block': {'gCaT': 0.8}}},
              {'label': 'hcn0.5_kv31x1.25', 'group': 'combined', 'x': 1, 'seeds': few,
               'cfg': {**pd, 'channel_block': {'gHCN': 0.5, 'gKv31': 1.25}}},
              {'label': 'cat0.9_kv31x1.25', 'group': 'combined', 'x': 2, 'seeds': few,
               'cfg': {**pd, 'channel_block': {'gCaT': 0.9, 'gKv31': 1.25}}},
              {'label': 'cat0.8_kv31x1.25', 'group': 'combined', 'x': 3, 'seeds': few,
               'cfg': {**pd, 'channel_block': {'gCaT': 0.8, 'gKv31': 1.25}}},
              {'label': 'kir0.9', 'group': 'single', 'x': 4, 'seeds': few, 'cfg': {**pd, 'channel_block': {'gKir': 0.9}}},
              {'label': 'kir0.8', 'group': 'single', 'x': 5, 'seeds': few, 'cfg': {**pd, 'channel_block': {'gKir': 0.8}}},
              {'label': 'kir0.9_kv31x1.25', 'group': 'combined', 'x': 4, 'seeds': few,
               'cfg': {**pd, 'channel_block': {'gKir': 0.9, 'gKv31': 1.25}}},
              {'label': 'kir0.8_kv31x1.25', 'group': 'combined', 'x': 5, 'seeds': few,
               'cfg': {**pd, 'channel_block': {'gKir': 0.8, 'gKv31': 1.25}}}]

# pinned to equal GPe->STN conductance in PV+ and PV- (default later moved to cfg.gpe_gaba_scale 1.13 / 0.87)
for c in conditions:
    c['cfg'] = {'gpe_gaba_scale': {'PVP': 1.0, 'PVN': 1.0}, **c['cfg']}
