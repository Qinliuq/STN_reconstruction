# H5 held-out run: healthy configuration D, 5 seeds x 20.5 s (validation/heldout/priors.json, H5_primate_change).
# The parkinsonian side reuses the H1 runs (H1_pd.py).
duration = 20500
seeds = [1, 2, 3, 4, 5]
outdir = 'sweeps/out/heldout/H5_healthy'
conditions = [{'label': 'healthy', 'cfg': {'ctx_state': 'healthy', 'gpe_state': 'healthy'}}]
