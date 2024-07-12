from netpyne import specs

cfg = specs.SimConfig()
cfg.hParams['celsius'] = 37 # to reproduce fast bursting, set to 36 (along with other parameters described in Fig 6A)
cfg.hParams['v_init'] = -60.0
cfg.duration = 750
cfg.seeds['conn'] = 4
cfg.dt = 0.05
cfg.recordStim = True
# cfg.recordCellsSpikes = [0, 1]
cfg.recordCells = [('PVP_pop', 0), ('PVN_pop', 0)]
cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'},
                    # 'gHCN': {'sec':'soma', 'loc':0.5, 'mech': 'Ih', 'var': 'gk'},
                    # 'gMaxHCN': {'sec':'soma', 'loc':0.5, 'mech': 'Ih', 'var': 'gmax_k'}
                    }

cfg.savePickle = True
# cfg.verbose = True
tRange = [0, cfg.duration]
cfg.analysis['plotTraces'] = {'saveFig': True, 'timeRange': tRange, 'oneFigPer': 'trace'}
cfg.analysis['plot2Dnet'] = {'include': ['PVP_pop', 'PVN_pop'], 'saveFig': True}
cfg.analysis['plotRaster'] = {'saveFig': True}