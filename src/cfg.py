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
sizeX = 8.5 * 1e3 # TODO: use var here and in netParams
sizeZ = 3 * 1e3
cfg.recordLFP = [[0.2*sizeX, 0.2*sizeX, sizeZ/2], 
                 [0.4*sizeX, 0.4*sizeX, sizeZ/2], 
                 [0.6*sizeX, 0.6*sizeX, sizeZ/2], 
                 [0.8*sizeX, 0.8*sizeX, sizeZ/2]]

cfg.savePickle = True
# cfg.verbose = True
tRange = [0, cfg.duration]
cfg.analysis['plotTraces'] = {'saveFig': True, 'timeRange': tRange, 'oneFigPer': 'trace'}
cfg.analysis['plot2Dnet'] = {'include': ['PVP_pop', 'PVN_pop'], 'saveFig': True, 'showConns': False}
cfg.analysis['plotRaster'] = {'saveFig': True}
cfg.analysis['plotLFP'] = {'plots': ['timeSeries',  'locations', 'PSD', 'spectrogram'], 'figSize': (7.5,13.5), 'saveFig': True}