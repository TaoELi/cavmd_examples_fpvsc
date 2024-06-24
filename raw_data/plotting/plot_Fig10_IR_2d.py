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

def obtain_spectrum_2d(path, idx_y_slice=1, pattern="simu_*.phsp", subspectrum=False):
    N = 18 # in one dimension, there are 18 cavity modes
    data = obtain_avg_data(path=path, pattern=pattern)
    freq = data[:,0]
    nstart = 1 + idx_y_slice*N
    sp = data[:,nstart:N+nstart]
    if subspectrum:
        # reset frequency starting point as 2000 cm-1, ending point as 3000 cm-1
        x = freq
        dx = x[2] - x[1]
        nstart, nend = int(1800 / dx), int(3000 / dx)
        freq = freq[nstart:nend]
        sp = sp[nstart:nend,:]

    # post processing
    x = freq
    dx = x[2] - x[1]
    nstart, nend = int(2200 / dx), int(2700 / dx)
    x = x[nstart:nend]
    sp = np.abs(sp[nstart:nend,:])
    print(np.max(np.max(sp)))
    #sp /= np.max(np.max(sp))
    sp /= 4.0e31
    sp = sp[::-1, :]
    freq_cav_inplane_min = 50.0
    freq_cav_inplane_max = 50.0 * N

    return freq, sp, freq_cav_inplane_min, freq_cav_inplane_max, x

def plot_IR_all():
    '''
    A function to plot Setup, Rabi splitting, and Avoid crossing at the same time
    '''
    fig, axes = clp.initialize(1, 3, width=4.3*2, height=4.3*0.7,  LaTeX=True, fontsize=12, return_fig_args=True, labelthem=True, labelsize=13, labelthemPosition=[-0.05, 1.06])

    # Fig one
    freq, sp, freq_cav_inplane_min, freq_cav_inplane_max, x = obtain_spectrum_2d(path="../co2_liquid_phase/E0_8.3e-6_2d", idx_y_slice=0)

    extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]

    from matplotlib.colors import LogNorm
    vmax = 1.0 #np.max(np.max(sp))
    vmin = vmax * 0.003
    pos = axes[0].imshow(sp, aspect='auto', extent=extent,
            cmap=cm.hot,
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )

    freq_cav_inplane = np.linspace(freq_cav_inplane_min, freq_cav_inplane_max, x.size)

    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * 2327, (2320.0**2 + freq_cav_inplane**2 + 50.0**2)**0.5]
    clp.plotone(xs, ys, axes[0], showlegend=False, colors=["w--", "g--"], lw=1.2,
            xlabel="$\omega_x$ [cm$^{-1}$]",
            ylabel="IR frequency [cm$^{-1}$]"),
    axes[0].text(519, 2351, "cavity photon", color='g', fontsize=10)
    axes[0].text(300, 2300, "C=O asym. stretch", color='w', fontsize=10)
    axes[0].tick_params(color='c', labelsize='medium', width=2)
    axes[0].text(150, 2600, "$\omega_y$ = 50 cm$^{-1}$", color='w', fontsize=10)

    # Fig two
    freq, sp, freq_cav_inplane_min, freq_cav_inplane_max, x = obtain_spectrum_2d(path="../co2_liquid_phase/E0_8.3e-6_2d", idx_y_slice=10)

    extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]

    from matplotlib.colors import LogNorm
    vmax = 1.0 #np.max(np.max(sp))
    vmin = vmax * 0.003
    pos = axes[1].imshow(sp, aspect='auto', extent=extent,
            cmap=cm.hot,
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )

    freq_cav_inplane = np.linspace(freq_cav_inplane_min, freq_cav_inplane_max, x.size)

    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * 2327, (2320.0**2 + freq_cav_inplane**2 + (50.0 * 11)**2)**0.5]
    clp.plotone(xs, ys, axes[1], showlegend=False, colors=["w--", "g--"], lw=1.2,
            xlabel="$\omega_x$ [cm$^{-1}$]")
            #ylabel="IR frequency [cm$^{-1}$]")
    axes[1].tick_params(color='c', labelsize='medium', width=2)
    axes[1].text(150, 2600, "$\omega_y$ = 550 cm$^{-1}$", color='w', fontsize=10)

    # Fig three
    freq, sp, freq_cav_inplane_min, freq_cav_inplane_max, x = obtain_spectrum_2d(path="../co2_liquid_phase/E0_8.3e-6_2d", idx_y_slice=17)

    extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]

    from matplotlib.colors import LogNorm
    vmax = 1.0 #np.max(np.max(sp))
    vmin = vmax * 0.003
    pos = axes[2].imshow(sp, aspect='auto', extent=extent,
            cmap=cm.hot,
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )

    freq_cav_inplane = np.linspace(freq_cav_inplane_min, freq_cav_inplane_max, x.size)

    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * 2327, (2320.0**2 + freq_cav_inplane**2 + (50.0 * 18)**2)**0.5]
    clp.plotone(xs, ys, axes[2], showlegend=False, colors=["w--", "g--"], lw=1.2,
            xlabel="$\omega_x$ [cm$^{-1}$]")
            #ylabel="IR frequency [cm$^{-1}$]")
    axes[2].tick_params(color='c', labelsize='medium', width=2)
    axes[2].text(150, 2600, "$\omega_y$ = 900 cm$^{-1}$", color='w', fontsize=10)

    # Additonal scaling of spectrum
    fig.subplots_adjust(right=0.95)
    cbar_ax = fig.add_axes([0.96, 0.15, 0.02, 0.70])

    cbar = fig.colorbar(pos, cax=cbar_ax)
    cbar.set_label('spectrum intensity [arb. units]')

    clp.adjust(savefile="IR_CO2_2D.pdf")


if __name__ == "__main__":
    plot_IR_all()
