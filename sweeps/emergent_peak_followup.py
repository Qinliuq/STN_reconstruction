# POST-HOC follow-ups to the emergent-peak test (results/emergent_peak/followup.md), requested after the Step 2 results.
# Same seeds, duration and criteria as sweeps/emergent_peak.py; parameters from params/jeon2022_profiles.json.
#
# FU1      configuration D (graded PD beta drive centre 0.65 / width 0.15, sigma 0.1, config-D receptors) with the
#          measured Pvalb PV+ placement. Tests Table 3's "the PV+ gradient is necessary" in its original form.
# FU2_F    F (measured placement, Gria1, Gabra1 + per-type ratio, uniform drive) with configuration D's GPe->STN
#          topography (sigma 0.1, evenly spaced homes) instead of sigma 0.219 + mapping.
# FU2_M*   M (graded drive at the somatomotor locus 0.572 / 0.634 / 0.697, drive width 0.219) with configuration D's
#          GPe->STN topography. Healthy reference = FU2_F healthy (healthy beta depth is uniform: modulation = floor).
# Dref_*   configuration D unchanged (healthy; PD graded; PD uniform drive): reruns of the Table 3/4 rows, for the
#          PV+/PV- pairwise shared-beta ratio (not stored in the earlier CSVs) and a reproducibility check.
import json
import os

P = json.load(open(os.path.join('params', 'jeon2022_profiles.json')))
duration = 3000
seeds = [1, 2, 3, 4, 5]
x_label = 'condition'
reference_group = 'healthy'
outdir = 'sweeps/out/emergent_peak_followup'

PLACE_DATA = {'pvp_placement': 'profile', 'pvp_profile': P['pvp_profile']}
F_REC = {**PLACE_DATA, 'ctx_topography': None, 'ctx_gscale': P['ctx_gscale'], 'gpe_gscale': P['gpe_gscale']}
M_LOCI = P['somatomotor_locus']['run_at']
M_WIDTH = P['somatomotor_locus']['width']

pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
healthy = {'ctx_state': 'healthy', 'gpe_state': 'healthy'}
COND = {
    'FU1': ({**PLACE_DATA}, ('healthy', 'pd')),
    'FU2_F': ({**F_REC, 'gpe_params': {'beta_width': 1e9}}, ('healthy', 'pd')),
    'Dref': ({}, ('healthy', 'pd')),
    'Dref_uni': ({'gpe_params': {'beta_width': 1e9}}, ('pd',)),
}
for j, loc in enumerate(M_LOCI):
    COND[f'FU2_M{j}'] = ({**F_REC, 'gpe_params': {'beta_center': loc, 'beta_width': M_WIDTH}}, ('pd',))

conditions = []
for i, (name, (c, states)) in enumerate(COND.items()):
    for st in states:
        conditions.append({'label': f'{name}_{st}', 'group': 'healthy' if st == 'healthy' else 'parkinsonian', 'x': i,
                           'cfg': {**(healthy if st == 'healthy' else pd), **c}})
