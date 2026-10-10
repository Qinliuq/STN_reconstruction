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


STN inputs (settings in `src/cfg.py`, glue code in `src/stn_inputs.py`):
- cortex: `src/cortical_arrivals.py` generates arrival events for every STN cell before the simulation (`cfg.ctx_state` = `'healthy'` or `'parkinsonian'`)
- GPe: `src/gpe_feedback.py` is closed-loop; `init.py` runs the simulation in steps of `cfg.gpe_dt`, passing the STN spikes to the GPe model and delivering its inhibitory arrivals. GPe group rates are saved as `<filename>_GPe_rates.png` and in the saved data (`GPe_rate_hz`, `GPe_arrival_t`, `GPe_arrival_stn`)

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