# Pre-run checks for sweeps/edge_artefacts.py: identity (all fix switches at their defaults reproduce
# drive_topography/ctl_1_uniform_control_pd_s1 exactly) and seed-1 smoke runs of R1, R2, R4 (operating point).
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runpy  # noqa: E402

E = runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'edge_artefacts.py'))
duration, x_label, reference_group = E['duration'], 'condition', 'healthy'
seeds = [1]
outdir = 'sweeps/out/edge_artefacts_check'
ident = E['pd_condition']('1_uniform_control', 'b')
ident = {**ident, 'label': 'ident_ctl_1_uniform_control_pd',
         'cfg': {**ident['cfg'], 'loc_grid': 'original', 'gpe_coverage_norm': False, 'gpe_drive_ref': None}}
conditions = [ident] + [E['pd_run'](r, '1_uniform_control') for r in ('R1', 'R2', 'R4')]
