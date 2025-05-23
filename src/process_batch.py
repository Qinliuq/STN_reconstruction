import os, sys, shutil
import numpy as np
import matplotlib.pyplot as plt

def build(sim_label):
    from featureExtr import load_from_file, extractSpikes, extractBursts

    os.chdir(sim_label)

    data = np.zeros((6, 20, 20, 3)) # shape (num-vals-CaL, num-vals-CaT, num-vals-HCN, num_features). There are 3 features: 0 - burst period (ms), 1 - duration of burst (ms), 2 - number of spikes in burst
    for ihcn in range(20):
        for icat in range(20):
            fname = f'{sim_label}_{ihcn}_{icat}_data.pkl'
            load_from_file(fname)
            for CaL in range(6):
                dat = extractBursts(CaL, transient=500)
                # alternatively (rare case), to use just spikes, not bursts:
                # dat = [len(extractSpikes(CaL, transient=500)), -1, -1] # last dimensions are dummy
                if dat == None:
                    dat = [-1, -1, -1]
                data[CaL, ihcn, icat] = dat
                print(dat)
    np.save('batch_summary.npy', data)

def plot_bursts_stats(sim_label):
    os.chdir(sim_label)
    data = np.load('batch_summary.npy') # shape (num-vals-CaL, num-vals-CaT, num-vals-HCN, num_features). There are 3 features: 0 - burst period (ms), 1 - duration of burst (ms), 2 - number of spiikes in burst 

    CaLs_to_use = [2, 3, 4] # range(data.shape[0])
    first_valid_HCN = 8 # 0

    fig, axs = plt.subplots(len(CaLs_to_use), 3)
    fig.set_size_inches(15, 15)

    # all CaLs: np.linspace(0.6, 1.6, cfg.num_vals)
    # all HCN (scale): np.linspace(0.5, 1.5, 20)
    # all CaT (scale): np.linspace(0.65, 2, 20)

    low_high = {}
    # to provide boundaries per each feature:
    # low_high = {
    #     'period': (200, 600),
    #     'duration': (0, 130),
    #     'num spikes': (0, 20)
    # }

    for CaL in CaLs_to_use:
        for feature, feat_name in zip([0, 1, 2], ['period', 'duration', 'num spikes']):
            this_dat = data[CaL,:,first_valid_HCN:,feature].T # transposed to get HCN-CaT, not CaT-HCN

            masked = np.ma.masked_where(this_dat == -1, this_dat)

            ax = axs[CaL-min(CaLs_to_use), feature]

            cmap = plt.cm.viridis
            cmap.set_bad(color='gray')

            vmin, vmax = low_high.get(feat_name, (None, None))
            img = ax.imshow(masked, cmap=cmap, interpolation='nearest', aspect=1.75, vmin=vmin, vmax=vmax)
            ax.invert_yaxis()
            ax.set_xlabel('CaT')
            ax.set_ylabel('HCN')
            ax.set_title(feat_name)
            plt.colorbar(img, ax=ax, orientation = 'vertical') # , label=feat_name

    plt.savefig('plot_burst_stats.png')


if __name__ == '__main__':
    sim_label = 'gCaT-gHCN' # 'kir-gCaT-gHCN'
    # build(sim_label)
    plot_bursts_stats(sim_label)


