# STN reconstruction

Requirements: NetPyNE and NEURON 8.2.x (NEURON 9 is not yet compatible with NetPyNE's LFP recording), e.g.
```
pip install "neuron==8.2.6" netpyne
```

First, compile .mod files (this should be done only once):
```
nrnivmodl mod
```

Launch the simulation:
```
python src/init.py
```


STN inputs (settings in `src/cfg.py`, glue code in `src/stn_inputs.py`). `cfg.input_version` selects the input models:
- `2` (default): `src/cortical_arrivals_v2.py` + `src/gpe_feedback_v2.py`: GPe is the primary beta source (parkinsonian pallidal beta drive), cortical beta is a supplement
- `1`: `src/cortical_arrivals.py` + `src/gpe_feedback.py`: beta comes from cortex (use with `cfg.ctx_overrides = {'conductance_nS': 1.0}`)

`cfg.ctx_state` and `cfg.gpe_state` (`'healthy'` / `'parkinsonian'`) set each input's state.
- cortex: arrival events for every STN cell are generated before the simulation; `cfg.ctx_topography` scales each cell's arrival rate along the STN axis
- GPe: closed-loop; `init.py` runs the simulation in steps of `cfg.gpe_dt`, passing the STN spikes to the GPe model and delivering its inhibitory arrivals. GPe group rates are saved as `<filename>_GPe_rates.png` and in the saved data (`GPe_rate_hz`, `GPe_arrival_t`, `GPe_arrival_stn`)

Each run also saves network metrics (rates, bursting, beta power, synchrony; see `src/network_metrics.py`) to `<filename>_metrics.json`.
`cfg.channel_block` applies a drug-like block to all STN cells, e.g. `{'gKir': 0.5}`.

Parameter sweeps (spec files in `sweeps/`, output in `sweeps/out/`, curated results in `results/`):
```
python src/sweep.py run sweeps/kir_block_gpe_beta.py --jobs 4
python src/sweep.py summary sweeps/kir_block_gpe_beta.py
```

Check single-cell behaviour of PV+ and PV- cells (no synaptic input; firing/burst stats and `single_cell.png`):
```
python src/single_cell.py
```
To sweep one conductance scale for one cell type (other scales as in `src/cfg.py`), e.g. CaT in PV- cells:
```
python src/single_cell.py --cell PVN --param gCaT --values 0.8 1.0 1.2 1.4
```

Ion channel conductance scales for PV+ and PV- cells are set in `src/cfg.py` (`cfg.PVP_*` / `cfg.PVN_*`), so they can be varied in `src/batch.py`.

Note: cells/sample.hoc and cells/tools.hoc are not used. They are here just for the reference to orig. model