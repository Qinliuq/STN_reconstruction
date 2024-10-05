import numpy as np
import matplotlib.pyplot as plt
from netpyne import sim

def load_from_file(name):
    sim.loadSimData(name)

def rate_per_cell(sim, duration, transient=0):
    num_cells = len(sim.net.cells)
    rates = np.zeros(num_cells)

    spks = np.array([sim.allSimData['spkid'], sim.allSimData['spkt']])
    for cellId in range(num_cells):
        this_cell_spiketimes = spks[1,spks[0]==cellId]
        this_cell_spiketimes = this_cell_spiketimes[this_cell_spiketimes > transient]
        rate = 1000 * len(this_cell_spiketimes) / (duration - transient)
        print(rate)
        rates[cellId] = rate
    return rates

def plot_all():
    fig, axs = plt.subplots(2, 2)

    def plot(ax, x, y, y_ax_lab):
        ax.plot(x, y)
        ax.scatter(x, y)
        ax.set_xlabel(y_ax_lab)
        ax.set_ylabel('Frequency (Hz)')
        ax.set_ylim(bottom=0)

    def num_cells():
        return len(sim.net.cells)

    load_from_file('iapps_data.pkl')
    iapps = np.linspace(0, 10, num_cells()) / 100
    rates = rate_per_cell(sim, 1500, 500)
    plot(axs[0, 0], iapps, rates, 'Iapp')

    load_from_file('cats_data.pkl')
    cats = np.linspace(15, 40, num_cells()) / 20
    rates = rate_per_cell(sim, 1500, 500)
    plot(axs[0, 1], cats, rates, 'gCaT')

    load_from_file('hcns_data.pkl')
    hcns = np.linspace(0, 4, num_cells()) / 2
    rates = rate_per_cell(sim, 1500, 500)
    plot(axs[1, 0], hcns, rates, 'gHCN')

    load_from_file('cals_data.pkl')
    cals = np.linspace(0, 40, num_cells()) / 5
    rates = rate_per_cell(sim, 1500, 500)
    plot(axs[1, 1], cals, rates, 'gCaL')

    fig.savefig('tonic_repr.png')

def extractBursts(icell):
    # duration = 2000
    transient = 500

    spks = np.array([sim.allSimData['spkid'], sim.allSimData['spkt']])

    this_cell_spiketimes = spks[1,spks[0]==icell]
    this_cell_spiketimes = this_cell_spiketimes[this_cell_spiketimes > transient]

    if len(this_cell_spiketimes) < 4: # to have at very least two bursts 
        print(f'Not enough spikes (N = {len(this_cell_spiketimes)})')
        return None

    # all_inds = np.arange(transient, duration, 1)
    # spk_inds = np.digitize(this_cell_spiketimes, all)

    # ts = np.zeros_like(all_inds)
    # ts[spk_inds] = 1

    ISI = np.diff(this_cell_spiketimes)

    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score, silhouette_samples

    kmeans = KMeans(n_clusters=2, n_init=10, random_state=123)

    ISI_ = ISI.reshape(-1, 1)
    kmeans.fit(ISI_)

    cluster_labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_.flatten()

    # Silhouette coefficients - mean and by clusters
    score = silhouette_score(ISI_, cluster_labels)
    dist = cluster_centers[1] - cluster_centers[0]
    # # TODO: check how good are thu clusts
    # if bad score or bad dist:
    #     return None

    bursts_clust_i = 0 if dist > 0 else 1
    burst_isi = ISI[cluster_labels == bursts_clust_i]

    thresh = burst_isi.max()
    
    # thresh_isi = cl

    bursts = []
    current_burst = []

    for i, sp in enumerate(this_cell_spiketimes[:-1]): # exclude last spike
        if ISI[i] <= thresh: # if next spike is close enough to this one, we have an ongoing burst
            current_burst.append(sp)                        
        else:
            if len(current_burst) > 0 and ISI[i-1] <= thresh:
                # this was the last spike in burst
                current_burst.append(sp)
                bursts.append(current_burst)
                current_burst = []
            else: # this was a single spike
                pass

    # handle last spike
    sp = this_cell_spiketimes[-1]
    if len(current_burst) > 0 and ISI[i-1] <= thresh:
        # this was the last spike in burst
        current_burst.append(sp)
        bursts.append(current_burst)

    bnum = len(bursts)
    if bnum == 0:
        return None

    bdur = 0
    bscount = 0
    period = 0
    for iburst, spikes in enumerate(bursts):
        bdur += spikes[-1] - spikes[0]
        bscount += len(spikes)
        if iburst+1 < len(bursts): # time from first spike uuntil next burst's first spike
            period += bursts[iburst+1][0] - spikes[0]

    avg_busrt_dur = bdur / bnum
    avg_sp_count = bscount / bnum
    avg_period = period / (bnum-1) if bnum > 1 else -1
    print(f'{bnum} bursts, avg. dur {round(avg_busrt_dur, 2)} sp. count {round(avg_sp_count, 2)}')
    return [avg_period, avg_busrt_dur, avg_sp_count]
