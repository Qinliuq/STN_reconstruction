"""
Single-cell check of PV+ and PV- STN cells (no synaptic input).

Each cell type is run through a set of current-clamp protocols, using the
conductance scales from cfg.py and the same setup as init.py (set_conductances,
ion styles, Bevan aCSF). Prints firing/burst statistics and saves voltage traces
to single_cell.png.

Run from the repo root (after `nrnivmodl mod`):
    python src/single_cell.py
"""
from netpyne import specs, sim
from neuron import h
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from cfg import cfg as netCfg # network config: conductance scales, temperature, v_init
from tools import set_conductances, set_ion_styles, apply_CSF_Bevan

duration = 4000
transient = 1000 # all cells fire a high-frequency transient right after init; exclude from stats

# name: list of IClamp pulses (del, dur, amp)
protocols = {
    'spontaneous (0 nA)': [],
    'depolarizing +0.05 nA': [(0, duration, 0.05)],
    'hold -0.16 nA': [(0, duration, -0.16)], # slow rhythmic bursting protocol used for PV+ 'best candidate'
    'hold -0.25 nA': [(0, duration, -0.25)],
    'hold -0.35 nA': [(0, duration, -0.35)],
    'rebound (-0.25 nA, 1000-1500 ms)': [(1000, 500, -0.25)],
}
rebound = {'protocol': 'rebound (-0.25 nA, 1000-1500 ms)', 'onset': 1000, 'release': 1500}

cellTypes = {'PVP_cell': 'PV+', 'PVN_cell': 'PV-'}

# ------------------------------------------------ SIMULATION ------------------------------------------------
cfg = specs.SimConfig()
cfg.hParams['celsius'] = netCfg.hParams['celsius']
cfg.hParams['v_init'] = netCfg.hParams['v_init']
cfg.duration = duration
cfg.dt = netCfg.dt
cfg.recordCells = ['all']
cfg.recordTraces = {'V_soma': {'sec': 'soma', 'loc': 0.5, 'var': 'v'}}
cfg.recordStep = 0.1

netParams = specs.NetParams()
resultCode = h.load_file('cells/SThprotocell.hoc')
if not resultCode:
    raise RuntimeError("Couldn't load prototype STN cell from `cells/STHprotocell.hoc")
protoCell = h.SThproto()

names = list(protocols)
for cellType in cellTypes:
    cellParams = netParams.importCellParams(label=cellType, fileName='cells/SThprotocell.hoc', cellName='SThcell', cellArgs=[0, protoCell])
    cellParams.secs.soma['threshold'] = -30
    pop = f'{cellType}_pop'
    netParams.popParams[pop] = {'cellType': cellType, 'numCells': len(names)} # one unconnected cell per protocol
    for i, name in enumerate(names):
        for j, (delay, dur, amp) in enumerate(protocols[name]):
            label = f'{pop}_{i}_{j}'
            netParams.stimSourceParams[label] = {'type': 'IClamp', 'del': delay, 'dur': dur, 'amp': amp}
            netParams.stimTargetParams[label] = {'source': label, 'conds': {'pop': pop, 'cellList': [i]}, 'sec': 'soma', 'loc': 0.5}

sim.create(netParams, cfg)
for cell in sim.net.cells:
    set_conductances(cell, netCfg)
set_ion_styles()
apply_CSF_Bevan()
sim.simulate()

# ------------------------------------------------ ANALYSIS ------------------------------------------------
def find_bursts(spikes, max_isi=20.0, min_spikes=3, pause_ratio=3.0):
    # burst = run of >= min_spikes spikes with ISI < max_isi, preceded and followed by a pause
    # longer than pause_ratio * mean intra-burst ISI (so regular fast firing is not a burst)
    isi = np.diff(spikes)
    runs, run = [], [0]
    for i, x in enumerate(isi):
        if x < max_isi:
            run.append(i + 1)
        else:
            runs.append(run)
            run = [i + 1]
    runs.append(run)

    bursts = []
    for run in runs:
        if len(run) < min_spikes or run[0] == 0 or run[-1] == len(spikes) - 1:
            continue # too short, or not bounded by spikes on both sides
        intra = np.mean(np.diff(spikes[run]))
        if min(isi[run[0] - 1], isi[run[-1]]) > pause_ratio * intra:
            bursts.append(spikes[run])
    return bursts

def firing_stats(spikes, t0, t1):
    inWin = spikes[(spikes >= t0) & (spikes < t1)]
    stats = {'rate': 1000 * len(inWin) / (t1 - t0), 'cv': np.nan, 'bursts': 0, 'fracInBursts': 0.0, 'spikesPerBurst': np.nan}
    if len(inWin) >= 3:
        isi = np.diff(inWin)
        stats['cv'] = isi.std() / isi.mean()
    bursts = [b for b in find_bursts(spikes) if t0 <= b[0] < t1]
    if bursts:
        stats['bursts'] = len(bursts)
        stats['fracInBursts'] = sum(len(b) for b in bursts) / len(inWin)
        stats['spikesPerBurst'] = np.mean([len(b) for b in bursts])
    return stats

def rebound_stats(spikes, release, max_isi=30.0):
    # spikes after release of hyperpolarization, until the first ISI > max_isi
    post = spikes[spikes >= release]
    if len(post) < 2:
        return len(post), 0.0
    isi = np.diff(post)
    end = int(np.argmax(isi > max_isi)) if (isi > max_isi).any() else len(isi)
    return end + 1, post[end] - post[0]

spikeTimes = {}
for gid, t in zip(sim.allSimData['spkid'], sim.allSimData['spkt']):
    spikeTimes.setdefault(int(gid), []).append(t)

def cell_spikes(cellType, i):
    gid = sim.net.pops[f'{cellType}_pop'].cellGids[i]
    return gid, np.array(sorted(spikeTimes.get(gid, [])))

print(f"\nCells without synaptic input; stats over {transient}-{duration} ms")
print("burst = >=3 spikes with ISI < 20 ms, flanked by pauses > 3x the intra-burst ISI\n")
for name in names:
    print(f"=== {name}")
    for cellType, lab in cellTypes.items():
        _, spikes = cell_spikes(cellType, names.index(name))
        if name == rebound['protocol']:
            pre = firing_stats(spikes, rebound['onset'] - 500, rebound['onset'])
            n, dur = rebound_stats(spikes, rebound['release'])
            late = firing_stats(spikes, rebound['release'] + 500, duration)
            print(f"  {lab}: before step {pre['rate']:5.1f} Hz | rebound {n:3d} spikes over {dur:4.0f} ms | late {late['rate']:5.1f} Hz")
        else:
            s = firing_stats(spikes, transient, duration)
            print(f"  {lab}: rate {s['rate']:5.1f} Hz | ISI CV {s['cv']:4.2f} | bursts {s['bursts']:3d} | "
                  f"spikes in bursts {s['fracInBursts']:4.2f} | spikes/burst {s['spikesPerBurst']:4.1f}")

# ------------------------------------------------ FIGURE ------------------------------------------------
t = np.array(sim.allSimData['t'])
fig, axs = plt.subplots(len(names), len(cellTypes), figsize=(14, 2.1 * len(names)), sharex=True, sharey=True, squeeze=False)
for c, (cellType, lab) in enumerate(cellTypes.items()):
    for i, name in enumerate(names):
        gid, _ = cell_spikes(cellType, i)
        axs[i, c].plot(t, sim.allSimData['V_soma'][f'cell_{gid}'], color=['tab:green', 'tab:blue'][c], lw=0.6)
        axs[i, c].set_title(f'{lab}: {name}', fontsize=9, loc='left')
        axs[i, c].set_ylim(-90, 40)
for ax in axs[:, 0]:
    ax.set_ylabel('Vm (mV)')
for ax in axs[-1]:
    ax.set_xlabel('Time (ms)')
fig.tight_layout()
fig.savefig('single_cell.png', dpi=110)
print("\nSaved single_cell.png")
