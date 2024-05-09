from netpyne import specs

cfg = specs.SimConfig()
cfg.hParams['celsius'] = 37
cfg.hParams['v_init'] = -65.0
cfg.duration = 50
cfg.recordCells = [0]
cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}

cfg.analysis['plotTraces'] = {'saveFig': True}
cfg.analysis['plotRaster'] = {'saveFig': True}