import numpy as np
import matplotlib.pyplot as plt
from netpyne import sim

def load_from_file(name):
    sim.initialize()
    sim.loadAll(name)

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

import os
os.chdir('bursts_37_cal0')
load_from_file('cat28_cal0.9_hcn1.3_sk0.9_iclamp_-0.25_data.pkl')
plot_all()