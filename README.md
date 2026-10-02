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


Ion channel conductance scales for PV+ and PV- cells are set in `src/cfg.py` (`cfg.PVP_*` / `cfg.PVN_*`), so they can be varied in `src/batch.py`.

Note: cells/sample.hoc and cells/tools.hoc are not used. They are here just for the reference to orig. model