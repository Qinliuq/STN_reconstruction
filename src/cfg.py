from netpyne import specs

cfg = specs.SimConfig()

cfg.duration = 100
cfg.recordCells = [0]
cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}

cfg.analysis['plotTraces'] = {'saveFig': True}
cfg.analysis['plotRaster'] = {'saveFig': True}