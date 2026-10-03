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