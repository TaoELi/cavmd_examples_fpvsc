import numpy as np
import columnplots as clp
from scipy import signal
import os, sys, time
import math
from itertools import islice
from scipy import fftpack
import glob
import json
import matplotlib.cm as cm
import matplotlib.pyplot as plt


def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def obtain_spectrum(path, n=1, pattern="simu_*.phsp", subspectrum=False):
    data = obtain_avg_data(path=path, pattern=pattern)
    freq = data[:,0]
    sp = data[:,1:n+1]
    if subspectrum:
        # reset frequency starting point as 2000 cm-1, ending point as 3000 cm-1
        x = freq
        dx = x[2] - x[1]
        nstart, nend = int(1800 / dx), int(3000 / dx)
        freq = freq[nstart:nend]
        sp = sp[nstart:nend,:]
    return freq, sp

def plot_IR_all():
    '''
    A function to plot Setup, Rabi splitting, and Avoid crossing at the same time
    '''
    fig, axes = clp.initialize(2, 3, width=4.3*3, height=4.3*0.7*2,  LaTeX=True, fontsize=12, return_fig_args=True, labelthem=True, labelsize=13, labelthemPosition=[-0.05, 1.06])

    def subplot(path, idx=0, text=r"cav-\#1.1.1"):
        N = 36
        E0 = "5e-5"
        freq, sp = obtain_spectrum(path=path, n=N)

        x = freq
        dx = x[2] - x[1]
        nstart, nend = int(2200 / dx), int(3000 / dx)
        x = x[nstart:nend]
        sp = np.abs(sp[nstart:nend, :])
        sp /= np.max(np.max(sp))
        sp = sp[::-1, :]
        freq_cav_min = (2320.0 ** 2 + 50 ** 2) ** 0.5
        freq_cav_max = (2320.0 ** 2 + (50 * N) ** 2) ** 0.5
        freq_cav_inplane_min = 50
        freq_cav_inplane_max = 50 * N
        extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0], x[-1]]

        from matplotlib.colors import LogNorm
        vmax = np.max(np.max(sp))
        vmin = vmax * 1e-3
        pos = axes[idx].imshow(sp, aspect='auto', extent=extent,
                             cmap=cm.hot,
                             interpolation='nearest',
                             norm=LogNorm(vmin=vmin, vmax=vmax)
                             )

        freq_cav_inplane = np.linspace(freq_cav_inplane_min, freq_cav_inplane_max, x.size)

        xs = [freq_cav_inplane] * 2
        ys = [np.ones(len(freq_cav_inplane)) * 2327, (2320.0 ** 2 + freq_cav_inplane ** 2) ** 0.5]
        clp.plotone(xs, ys, axes[idx], showlegend=False, colors=["w--", "g--"], lw=1.2,
                    xlabel="cavity in-plane frequency [cm$^{-1}$]" if idx[0]==1 else None,
                    ylabel="IR frequency [cm$^{-1}$]" if idx[1]==0 else None)
        axes[idx].text(1100, 2484, "cavity photon", color='g', fontsize=12)
        axes[idx].text(1000, 2355, "C=O asym. stretch", color='w', fontsize=12)
        axes[idx].text(200, 2900, text, color='w', fontsize=12)
        axes[idx].tick_params(color='c', labelsize='medium', width=2)
        return pos

    # Figure a: Cavity detuning dependence of the IR spectrum
    pos = subplot(path="../co2_liquid_phase_systemsize/mode_times_4/", idx=(0,0), text=r"cav-\#1.4.4")

    # Figure b: Cavity detuning dependence of the IR spectrum
    pos = subplot(path="../co2_liquid_phase_systemsize/size_times_16/", idx=(0,1), text=r"cav-\#1.16.16")

    # Figure c: Cavity detuning dependence of the IR spectrum
    pos = subplot(path="../co2_liquid_phase_systemsize/size_times_64/", idx=(0, 2), text=r"cav-\#1.64.64")

    # Figure c: Cavity detuning dependence of the IR spectrum
    pos = subplot(path="../co2_liquid_phase_systemsize/mode_times_1/", idx=(1,0), text=r"cav-\#1.4.1")

    # Figure d: Cavity detuning dependence of the IR spectrum
    pos = subplot(path="../co2_liquid_phase_systemsize/mode_times_6/", idx=(1,1), text=r"cav-\#1.4.6")

    # Figure d: Cavity detuning dependence of the IR spectrum
    pos = subplot(path="../co2_liquid_phase_systemsize/mode_times_8/", idx=(1,2), text=r"cav-\#1.4.8")

    # Additonal scaling of spectrum
    fig.subplots_adjust(right=0.95)
    cbar_ax = fig.add_axes([0.96, 0.15, 0.01, 0.70])

    cbar = fig.colorbar(pos, cax=cbar_ax)
    cbar.set_label('spectrum intensity [arb. units]')

    clp.adjust(savefile="IR_CO2_mode_dep.pdf")


if __name__ == "__main__":
    plot_IR_all()
