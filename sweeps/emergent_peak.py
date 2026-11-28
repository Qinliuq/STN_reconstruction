# Can the parkinsonian beta peak at STN axis ~0.6-0.7 emerge from measured distributions alone (results/emergent_peak)?
# PD state with a UNIFORM pallidal beta drive (full depth in all GPe groups, as in Table 3; beta_width 1e9); healthy
# runs per condition as reference. All spatial parameters come from params/jeon2022_profiles.json (Jeon et al. 2022),
# committed before any of these runs. Configuration D itself is unchanged (every option here is off by default).
#
# F      PV+ placement from the Pvalb profile, Gria1 glutamate conductance (rate topography off), Gabra1 GABA with the
#        per-type ratio (no double counting), data GPe->STN sigma + mapping
# F_cfgD F with configuration D's (skew-normal, exaggerated) PV+ placement
# noPV   F with uniform PV+ placement
# noRec  F with flat glutamate and flat GABA (per-type ratio kept)
# abs    F without the per-type ratio (absolute Gria1 / Gabra1 gradients only)
# sig    measured placement, configuration-D glutamate (0.65 - 0.5 axis, rate) and flat GABA with config-D per-type scale,
#        data GPe->STN sigma + mapping
# ref    as sig with configuration D's GPe->STN topography (sigma 0.1, even homes)
# M_*    as F, but graded PD beta drive centred on the somatomotor hyperdirect locus (Fig. 3C) +- its reading error,
#        width = data sigma; labelled "assumes motor-territory origin of beta", not part of the emergent-peak verdict.
#        Healthy reference = F healthy (the healthy beta depth is uniform, so the centre/width do not change it).
import json
import os

P = json.load(open(os.path.join('params', 'jeon2022_profiles.json')))
duration = 3000
seeds = [1, 2, 3, 4, 5]
x_label = 'condition'
reference_group = 'healthy'
outdir = 'sweeps/out/emergent_peak'

UNIFORM_BETA = {'beta_width': 1e9}
DATA_TOPO = {'topo_sigma': P['gpe_topography']['topo_sigma'], 'group_center_map': P['gpe_topography']['group_center_map']}
PLACE_DATA = {'pvp_placement': 'profile', 'pvp_profile': P['pvp_profile']}
GLUT_DATA = {'ctx_topography': None, 'ctx_gscale': P['ctx_gscale']}
GLUT_FLAT = {'ctx_topography': None}

COND = {
    'F': {**PLACE_DATA, **GLUT_DATA, 'gpe_gscale': P['gpe_gscale'], 'gpe_params': {**DATA_TOPO, **UNIFORM_BETA}},
    'F_cfgD': {'pvp_placement': 'gradient', **GLUT_DATA, 'gpe_gscale': P['gpe_gscale'], 'gpe_params': {**DATA_TOPO, **UNIFORM_BETA}},
    'noPV': {'pvp_placement': 'uniform', **GLUT_DATA, 'gpe_gscale': P['gpe_gscale'], 'gpe_params': {**DATA_TOPO, **UNIFORM_BETA}},
    'noRec': {**PLACE_DATA, **GLUT_FLAT, 'gpe_gscale': P['gpe_gscale_flat_with_ratio'], 'gpe_params': {**DATA_TOPO, **UNIFORM_BETA}},
    'abs': {**PLACE_DATA, **GLUT_DATA, 'gpe_gscale': P['gpe_gscale_abs'], 'gpe_params': {**DATA_TOPO, **UNIFORM_BETA}},
    'sig': {**PLACE_DATA, 'gpe_params': {**DATA_TOPO, **UNIFORM_BETA}},
    'ref': {**PLACE_DATA, 'gpe_params': {**UNIFORM_BETA}},
}
M_LOCI = P['somatomotor_locus']['run_at']
M_WIDTH = P['somatomotor_locus']['width']

pd = {'ctx_state': 'parkinsonian', 'gpe_state': 'parkinsonian'}
healthy = {'ctx_state': 'healthy', 'gpe_state': 'healthy'}
conditions = []
for i, (name, c) in enumerate(COND.items()):
    conditions.append({'label': f'{name}_healthy', 'group': 'healthy', 'x': i, 'cfg': {**healthy, **c}})
    conditions.append({'label': f'{name}_pd', 'group': 'parkinsonian', 'x': i, 'cfg': {**pd, **c}})
for j, loc in enumerate(M_LOCI):
    c = dict(COND['F'])
    c['gpe_params'] = {**DATA_TOPO, 'beta_center': loc, 'beta_width': M_WIDTH}  # state's beta_modulation 1.0 / floor 0.2
    conditions.append({'label': f'M{j}_pd', 'group': 'parkinsonian', 'x': len(COND) + j, 'cfg': {**pd, **c}})
