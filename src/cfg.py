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
# network dimensions (um), also used by netParams.py
cfg.sizeX = 850 # 8.5 * 1e3
cfg.sizeY = 850 # 8.5 * 1e3
cfg.sizeZ = 100 # 3 * 1e3
# LFP electrodes placed along the STN (diagonal) axis
cfg.recordLFP = [[f*cfg.sizeX, f*cfg.sizeY, cfg.sizeZ/2] for f in (0.2, 0.4, 0.6, 0.8)]

# ion channel conductance scaling per cell type (relative to default values in cells/sth-data)
# PV+ "best candidate" (bursting); this combo gives burst-pause like behavior during 20Hz input: cat, cal, hcn, sk = 1.4, 0.86, 0.65, 0.9
cfg.PVP_gCaT_scale = 1.4
cfg.PVP_gCaL_scale = 0.85
cfg.PVP_gHCN_scale = 0.65
cfg.PVP_gSK_scale = 0.9
cfg.PVP_gKir_scale = 0.5
# PV- (single-spiking)
cfg.PVN_gCaT_scale = 1.2 # tonic firing at rest (silent at <= 1.0, bursts at -0.16 nA at 1.4)
cfg.PVN_gCaL_scale = 0.8 # short rebound (1.2 gives a ~250 ms rebound plateau)
cfg.PVN_gHCN_scale = 0.6
cfg.PVN_gSK_scale = 1.0
cfg.PVN_gKir_scale = 0.1

cfg.savePickle = True
# cfg.verbose = True
tRange = [0, cfg.duration]
cfg.analysis['plotTraces'] = {'saveFig': True, 'timeRange': tRange, 'oneFigPer': 'trace'}
cfg.analysis['plot2Dnet'] = {'include': ['PVP_pop', 'PVN_pop'], 'saveFig': True, 'showConns': False}
cfg.analysis['plotRaster'] = {'saveFig': True}
cfg.analysis['plotLFP'] = {'plots': ['timeSeries',  'locations', 'PSD', 'spectrogram'], 'figSize': (7.5,13.5), 'saveFig': True}