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
    fig, axes = clp.initialize(2, 2, width=4.3*2, height=4.3*0.7*2,  LaTeX=True, fontsize=12, return_fig_args=True, labelthem=True, labelsize=13, labelthemPosition=[-0.05, 1.06])

    # Figure a: Linear IR spectroscopy outside vs inside the cavity
    xs, ys = [], []
    freq, sp = obtain_spectrum(path="../co2_liquid_phase/E0_0e-4", pattern="simu_*.sp",  subspectrum=True)
    xs.append(freq)
    ys.append(sp / np.max(sp))
    E0 = "5e-5"
    paths = ["../co2_liquid_phase/E0_%s_singlemode" %E0, "../co2_liquid_phase/E0_%s" %E0]
    for i, path in enumerate(paths):
        freq, sp = obtain_spectrum(path=path, subspectrum=True)
        xs.append(freq)
        ys.append(sp / np.max(sp) +  (i+1)*1.)
    colors = ["k"] * len(path)
    clp.plotone(xs, ys, axes[0,0], colors=colors, showlegend=False, lw=1, xlim=[1900, 2650], ylim=[-0.01, 3.05],
                #xlabel="frequency [cm$^{-1}$]", 
                ylabel="spectrum [arb. units]")
    colormap = plt.cm.hot
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(axes[0,0].lines))]
    for i,j in enumerate(axes[0,0].lines):
        j.set_color(colors[i])

    axes[0,0].text(1910, 0.1, "liquid CO$_2$\noutside cavity", fontsize=11, color=colors[0])
    axes[0,0].text(1910, 1.1, "single-mode cavity\n36 CO$_2$ subsystems",fontsize=11, color=colors[1])
    axes[0,0].text(1910, 2.1, "36-mode cavity\n36 CO$_2$ subsystems",fontsize=11, color=colors[2])
    axes[0,0].text(1910, 2.7, "Classical",fontsize=12)
    axes[0,0].axvline(x=2320, lw=3, color=clp.dark_green, alpha=0.5)

    # Figure b: Cavity detuning dependence of the IR spectrum
    N = 36
    E0 = "5e-5"
    freq, sp = obtain_spectrum(path="../co2_liquid_phase/E0_%s" %E0, n=36)

    x = freq
    dx = x[2] - x[1]
    nstart, nend = int(2200 / dx), int(3000 / dx)
    x = x[nstart:nend]
    sp = np.abs(sp[nstart:nend,:])
    sp /= np.max(np.max(sp)) 
    sp = sp[::-1, :]
    freq_cav_min = (2320.0**2 + 50.0**2)**0.5
    freq_cav_max = (2320.0**2 + (50.0*N)**2)**0.5
    freq_cav_inplane_min = 50.0
    freq_cav_inplane_max = 50.0 * N
    extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]

    from matplotlib.colors import LogNorm
    vmax = np.max(np.max(sp))
    vmin = vmax * 0.003
    pos = axes[0,1].imshow(sp, aspect='auto', extent=extent,  
            cmap=cm.hot,
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )
    
    freq_cav_inplane = np.linspace(freq_cav_inplane_min, freq_cav_inplane_max, x.size)
    
    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * 2327, (2320.0**2 + freq_cav_inplane**2)**0.5]
    clp.plotone(xs, ys, axes[0,1], showlegend=False, colors=["w--", "g--"], lw=1.2,
            #xlabel="cavity in-plane frequency [cm$^{-1}$]",
            ylabel="IR frequency [cm$^{-1}$]")
    axes[0,1].text(1250, 2550, "cavity photon", color='g', fontsize=12)
    axes[0,1].text(1000, 2370, "C=O asym. stretch", color='w', fontsize=12)
    axes[0,1].tick_params(color='c', labelsize='medium', width=2)

    # Figure c: Linear IR spectroscopy outside vs inside the cavity
    xs, ys = [], []
    freq, sp = obtain_spectrum(path="../co2_liquid_phase/pirerun_E0_0e-4", pattern="simu_*.sp",  subspectrum=True)
    xs.append(freq)
    ys.append(sp / np.max(sp))
    E0 = "5e-5"
    paths = ["../co2_liquid_phase/pi_E0_%s_singlemode" %E0, "../co2_liquid_phase/pi_E0_%s" %E0]
    for i, path in enumerate(paths):
        freq, sp = obtain_spectrum(path=path, subspectrum=True)
        xs.append(freq)
        ys.append(sp / np.max(sp) +  (i+1)*1.)
    colors = ["k"] * len(path)
    clp.plotone(xs, ys, axes[1,0], colors=colors, showlegend=False, lw=1, xlim=[1900, 2650], ylim=[-0.01, 3.05],
                xlabel="frequency [cm$^{-1}$]", ylabel="spectrum [arb. units]")
    colormap = plt.cm.hot
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(axes[1,0].lines))]
    for i,j in enumerate(axes[1,0].lines):
        j.set_color(colors[i])

    #axes[1,0].text(1910, 0.1, "liquid CO$_2$\noutside cavity", fontsize=11, color=colors[0])
    #axes[1,0].text(1910, 1.1, "single-mode cavity",fontsize=11, color=colors[1])
    #axes[1,0].text(1910, 2.1, "multimode cavity",fontsize=11, color=colors[2])
    axes[1,0].text(1910, 2.7, "TRPMD",fontsize=12)
    axes[1,0].axvline(x=2320, lw=3, color=clp.dark_green, alpha=0.5)

    # Figure d: Cavity detuning dependence of the IR spectrum
    N = 36
    E0 = "5e-5"
    freq, sp = obtain_spectrum(path="../co2_liquid_phase/pi_E0_%s" %E0, n=36)

    x = freq
    dx = x[2] - x[1]
    nstart, nend = int(2200 / dx), int(3000 / dx)
    x = x[nstart:nend]
    sp = np.abs(sp[nstart:nend,:])
    sp /= np.max(np.max(sp)) 
    sp = sp[::-1, :]
    freq_cav_min = (2320.0**2 + 50.0**2)**0.5
    freq_cav_max = (2320.0**2 + (50.0*N)**2)**0.5
    freq_cav_inplane_min = 50.0
    freq_cav_inplane_max = 50.0 * N
    extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]

    from matplotlib.colors import LogNorm
    #vmax = np.max(np.max(sp))
    #vmin = vmax * 0.001
    pos = axes[1,1].imshow(sp, aspect='auto', extent=extent,  
            cmap=cm.hot,
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )
    
    freq_cav_inplane = np.linspace(freq_cav_inplane_min, freq_cav_inplane_max, x.size)
    
    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * 2309, (2320.0**2 + freq_cav_inplane**2)**0.5]
    clp.plotone(xs, ys, axes[1,1], showlegend=False, colors=["w--", "g--"], lw=1.2,
            xlabel="cavity in-plane frequency $\omega_{\parallel}$ [cm$^{-1}$]",
            ylabel="IR frequency [cm$^{-1}$]")
    #axes[1,1].text(1250, 2550, "cavity photon", color='g', fontsize=12)
    #axes[1,1].text(1000, 2370, "C=O asym. stretch", color='w', fontsize=12)
    axes[1,1].tick_params(color='c', labelsize='medium', width=2)

    # Additonal scaling of spectrum
    fig.subplots_adjust(right=0.95)
    cbar_ax = fig.add_axes([0.96, 0.15, 0.02, 0.70])

    cbar = fig.colorbar(pos, cax=cbar_ax)
    cbar.set_label('spectrum intensity [arb. units]')

    clp.adjust(savefile="IR_CO2.pdf")


if __name__ == "__main__":
    plot_IR_all()
