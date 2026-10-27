"""
Connect the STN input models to the NetPyNE network:

- cortical_arrivals.py (open loop): arrival events for every STN cell are generated
  before the simulation; each STN cell gets one VecStim with its own event times,
  connected 1:1 to a 'glut' synapse on the soma. Optionally (cfg.ctx_topography) the
  arrival rate is scaled per cell by its position along the STN axis. Called from netParams.py.
- gpe_feedback.py (closed loop): GPeInputGenerator is stepped during the simulation
  (every cfg.gpe_dt) with the STN spikes of the previous step; its inhibitory arrivals
  are delivered to 'GABA' synapses on the STN cells. Called from init.py.

cfg.input_version selects the input models: 1 = cortical_arrivals.py + gpe_feedback.py (cortical beta),
2 = cortical_arrivals_v2.py + gpe_feedback_v2.py (GPe is the primary beta source, cortical beta a
supplement), 3 = cortical_arrivals_v2.py + gpe_feedback_v3.py (as 2, with topographic GPe-STN connections
and parkinsonian beta graded along the STN axis). cfg.ctx_state / cfg.gpe_state choose each model's 'healthy' or 'parkinsonian' state
(gpe_state applies to v2 only).

Both generators index STN cells 0..nSTN-1 in gid order: PV+ cells first, then PV-.
Conductances (nS per arrival) are used as NetCon weights in uS.
Single-process runs only (the GPe coupling reads spikes on the local node).
"""
import math
import dataclasses

import numpy as np
from neuron import h

import cortical_arrivals
import cortical_arrivals_v2
import gpe_feedback
import gpe_feedback_v2
import gpe_feedback_v3

STN_POPS = ['PVP_pop', 'PVN_pop']


def input_models(cfg):
    if cfg.input_version == 1:
        return cortical_arrivals, gpe_feedback
    if cfg.input_version == 2:
        return cortical_arrivals_v2, gpe_feedback_v2
    if cfg.input_version == 3:
        return cortical_arrivals_v2, gpe_feedback_v3
    raise ValueError(f'unknown cfg.input_version {cfg.input_version}')


# ------------------------------------------------ CORTEX ------------------------------------------------
def axis_position(loc):
    # position along the STN axis (the 45-deg diagonal), 0..1, from normalized cell coordinates
    return float(np.clip((loc['xnorm'] + loc['ynorm']) / 2, 0, 1))

def topography_scale(cellLocs, cfg):
    # per-cell cortical rate scale from the position along the STN axis: scale = intercept + slope * axis,
    # normalized to mean 1
    axis = np.array([axis_position(loc) for loc in cellLocs])
    if cfg.ctx_topography is None:
        return axis, np.ones(len(axis))
    intercept, slope = cfg.ctx_topography
    scale = intercept + slope * axis
    if np.any(scale <= 0):
        raise ValueError(f'cfg.ctx_topography {cfg.ctx_topography} gives a non-positive rate scale along the STN axis')
    return axis, scale / scale.mean()

def separate_coincident(times_ms, step=1e-6):
    # VecStim replays coincident events back to back; with float rounding the 2nd can fall a hair before t and
    # NEURON aborts (net_send td-t < 0). Spread coincident events by `step` ms (1 ns) so times strictly increase.
    t = np.asarray(times_ms, dtype=float)
    if len(t) < 2:
        return t
    runId = np.cumsum(np.r_[True, np.diff(t) > step / 2])
    rank = np.arange(len(t)) - np.searchsorted(runId, runId)
    return t + rank * step

def add_cortical_input(netParams, cfg, cellLocsPerPop):
    # cellLocsPerPop: {pop: cellsList with 'xnorm'/'ynorm'}, in STN index order (PVP_pop, PVN_pop)
    ctxModel, _ = input_models(cfg)
    params = dataclasses.replace(ctxModel.STATES[cfg.ctx_state], **cfg.ctx_overrides)
    cellLocs = [loc for locs in cellLocsPerPop.values() for loc in locs]
    nSTN = len(cellLocs)
    nAssemblies = max(1, math.ceil(nSTN / cfg.ctx_cells_per_assembly))
    axis, scale = topography_scale(cellLocs, cfg)

    # generate at the highest per-cell rate, then thin each cell's events to its own rate
    # (thinning keeps the arrivals Poisson and keeps the shared fraction of each cell's input)
    time_s, neuron_id, _ = ctxModel.generate_arrivals(n_neurons=nSTN, duration_s=cfg.duration / 1000, dt_s=cfg.ctx_dt,
                                             params=dataclasses.replace(params, rate_hz=params.rate_hz * scale.max()),
                                             n_assemblies=nAssemblies, seed=cfg.ctx_seed)
    keep = np.random.default_rng([cfg.ctx_seed, 1]).random(len(time_s)) < (scale / scale.max())[neuron_id]
    time_ms, neuron_id = time_s[keep] * 1000, neuron_id[keep]

    first = 0
    for pop, locs in cellLocsPerPop.items():
        n = len(locs)
        spkTimes = [separate_coincident(time_ms[neuron_id == first + i]).tolist() for i in range(n)] # sorted by generate_arrivals
        netParams.popParams[f'Ctx_{pop}'] = {'cellModel': 'VecStim', 'numCells': n, 'spkTimes': spkTimes}
        netParams.connParams[f'Ctx->{pop}'] = {'preConds': {'pop': f'Ctx_{pop}'}, 'postConds': {'pop': pop},
                                               'connList': [[i, i] for i in range(n)], 'synMech': 'glut',
                                               'weight': params.conductance_nS * 1e-3, 'delay': cfg.ctx_delay,
                                               'sec': 'soma', 'loc': 0.5}
        first += n
    rates = np.bincount(neuron_id, minlength=nSTN) / (cfg.duration / 1000)
    print(f"Cortical input (v{cfg.input_version}, {cfg.ctx_state}): {len(time_ms)} arrivals to {nSTN} STN cells "
          f"({rates.mean():.1f} /s per cell, {nAssemblies} assemblies)")
    if cfg.ctx_topography is not None:
        print(f"  topography {cfg.ctx_topography}: rate scale {scale.min():.2f}-{scale.max():.2f} along the STN axis "
              f"(arrivals {rates.min():.0f}-{rates.max():.0f} /s; corr(axis, arrivals) = {np.corrcoef(axis, rates)[0, 1]:.2f})")


# ------------------------------------------------ GPe ------------------------------------------------
class GPeFeedback:
    def __init__(self, sim, cfg, synMech):
        self.sim = sim
        self.dt = cfg.gpe_dt * 1000 # ms
        gids = [gid for pop in STN_POPS for gid in sim.net.pops[pop].cellGids]
        self.stnIndex = {gid: i for i, gid in enumerate(gids)}
        # arrivals for a step are computed from the previous step's spikes, which adds one step to the feedback delay
        _, gpeModel = input_models(cfg)
        stateParams = {**(gpeModel.STATES[cfg.gpe_state] if hasattr(gpeModel, 'STATES') else {}), **cfg.gpe_params}
        synapseScale = stateParams.get('synapse_scale', 1.0) # GPe->STN conductance multiplier of the state (v3)
        cellsByGid = {cell.gid: cell for cell in sim.net.cells}
        if cfg.input_version >= 3:
            stateParams['stn_axis'] = [axis_position(cellsByGid[gid].tags) for gid in gids]
        self.generator = gpeModel.GPeInputGenerator(n_stn=len(gids), dt_s=cfg.gpe_dt,
                                                    feedback_delay_s=max(0.0, cfg.gpe_feedback_delay - cfg.gpe_dt),
                                                    seed=cfg.gpe_seed, **stateParams)
        self.rng = np.random.default_rng(cfg.gpe_seed)

        # GABA synapses: on all dendritic sections for PV+ cells, on the soma for PV- cells (as in the previous GPe->STN rules);
        # each arrival goes to one randomly chosen synapse of the target cell
        synMech = {**synMech, **cfg.gpe_syn} # cfg.gpe_syn overrides the GABA synapse ('e', 'tau1', 'tau2')
        # per-type GPe->STN conductance factor (GABA:glutamate balance per cell type), cfg.gpe_gaba_scale
        gabaScale = {'PVP_cell': cfg.gpe_gaba_scale.get('PVP', 1.0), 'PVN_cell': cfg.gpe_gaba_scale.get('PVN', 1.0)}
        self.syns, self.netcons = [], []
        for gid in gids:
            cell = cellsByGid[gid]
            secNames = [s for s in cell.secs if s != 'soma'] if cell.tags['cellType'] == 'PVP_cell' else ['soma']
            ncs = []
            for secName in secNames:
                syn = h.Exp2Syn(cell.secs[secName]['hObj'](0.5))
                syn.e, syn.tau1, syn.tau2 = synMech['e'], synMech['tau1'], synMech['tau2']
                nc = h.NetCon(None, syn)
                nc.weight[0] = cfg.gpe_conductance_nS * synapseScale * gabaScale[cell.tags['cellType']] * 1e-3
                self.syns.append(syn)
                ncs.append(nc)
            self.netcons.append(ncs)

        self.nSpikesRead = 0
        self.rates, self.arrivalTimes, self.arrivalCells = [], [], []

    def step(self, t):
        # STN spikes recorded since the last call (spike vectors are filled during the run)
        spkid = self.sim.simData['spkid']
        n = int(spkid.size())
        newGids = spkid.c(self.nSpikesRead, n - 1).to_python() if n > self.nSpikesRead else []
        newIds = [self.stnIndex[int(g)] for g in newGids if int(g) in self.stnIndex]
        self.nSpikesRead = n

        arrivals = self.generator.step(newIds)
        for i in arrivals:
            ncs = self.netcons[i]
            ncs[self.rng.integers(len(ncs))].event(t)
        # emitted rate (includes the v2 beta modulation; equals rate_hz in v1)
        self.rates.append(getattr(self.generator, 'emitted_rate_hz', self.generator.rate_hz).tolist())
        self.arrivalTimes.extend([t] * len(arrivals))
        self.arrivalCells.extend(arrivals.tolist())

    def run(self):
        # like sim.runSimWithIntervalFunc, but also steps the GPe model at t = 0
        from netpyne.sim.run import prepareSimWithIntervalFunc, runForInterval, postRun
        stopTime, _ = prepareSimWithIntervalFunc()
        self.step(h.t)
        while round(h.t) < stopTime:
            runForInterval(self.dt, self.step)
        postRun(stopTime)

    def save(self):
        # call after sim.gatherData(); stored in the saved simulation data
        self.sim.allSimData['GPe_rate_hz'] = self.rates # one row (GPe groups) per step of cfg.gpe_dt, from t = 0
        self.sim.allSimData['GPe_arrival_t'] = self.arrivalTimes # ms
        self.sim.allSimData['GPe_arrival_stn'] = self.arrivalCells # STN index (PV+ first, then PV-)

    def plot(self, fileName):
        import matplotlib.pyplot as plt
        rates = np.array(self.rates)
        t = np.arange(len(rates)) * self.dt
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.plot(t, rates, lw=0.8)
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('GPe group rate (Hz)')
        ax.set_title(f'GPe groups ({rates.shape[1]}), {len(self.arrivalTimes)} arrivals to STN', fontsize=9)
        fig.tight_layout()
        fig.savefig(fileName, dpi=110)
        plt.close(fig)
