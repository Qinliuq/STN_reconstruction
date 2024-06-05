from netpyne import specs

cfg = specs.SimConfig()
cfg.hParams['celsius'] = 37
cfg.hParams['v_init'] = -60.0
cfg.duration = 6000
cfg.recordCells = [0]
cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}

cfg.savePickle = True
tRange = [1000, cfg.duration]
cfg.analysis['plotTraces'] = {'saveFig': True, 'timeRange': tRange}
cfg.analysis['plotRaster'] = {'saveFig': True}