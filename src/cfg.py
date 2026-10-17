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

# STN->STN collateral excitation (uS per connection)
cfg.stn_stn_weight = 0.0002 # 0.001 gave a loop-driven 6-8 Hz bursting state that masked the input differences

# STN inputs
# cortex (src/cortical_arrivals.py): open-loop arrival events, generated before the simulation
cfg.ctx_state = 'healthy' # 'healthy' or 'parkinsonian' (STATES in cortical_arrivals.py)
# override fields of the state's InputParams (applied to both states); 1 nS per arrival (script default 0.2 nS) lets the
# parkinsonian beta modulation reach STN: ~2-3x healthy beta power and ~4x synchrony
cfg.ctx_overrides = {'conductance_nS': 1.0}
# topography: per-cell rate scale = intercept + slope * position along the STN axis (0..1), normalized to mean 1;
# None for uniform input. Jeon et al. Fig 7H.2 (Gria): 0.55 - 0.3x; default is the previous, exaggerated rule 0.65 - 0.5x
cfg.ctx_topography = (0.65, -0.5)
cfg.ctx_cells_per_assembly = 31 # STN cells sharing one cortical assembly (2000 cells / 64 assemblies in cortical_arrivals.py)
cfg.ctx_dt = 0.001 # s, time bin of the arrival generator
cfg.ctx_delay = 1 # ms, synaptic delay
cfg.ctx_seed = 1
# GPe (src/gpe_feedback.py): closed loop, stepped during the simulation with the STN spikes of the previous step
cfg.gpe_dt = 0.001 # s, GPe update step
cfg.gpe_feedback_delay = 0.005 # s, total STN->GPe feedback delay (includes the one-step lag of the coupling)
cfg.gpe_conductance_nS = 0.5 # per inhibitory arrival
cfg.gpe_params = {} # other GPeInputGenerator arguments, e.g. {'n_gpe_groups': 12, 'target_probability': 0.35, 'stn_to_gpe_gain': 0.5}
cfg.gpe_seed = 2

cfg.savePickle = True
# cfg.verbose = True
tRange = [0, cfg.duration]
cfg.analysis['plotTraces'] = {'saveFig': True, 'timeRange': tRange, 'oneFigPer': 'trace'}
cfg.analysis['plot2Dnet'] = {'include': ['PVP_pop', 'PVN_pop'], 'saveFig': True, 'showConns': False}
cfg.analysis['plotRaster'] = {'saveFig': True}
cfg.analysis['plotLFP'] = {'plots': ['timeSeries',  'locations', 'PSD', 'spectrogram'], 'figSize': (7.5,13.5), 'saveFig': True}