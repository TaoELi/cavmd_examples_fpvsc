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

def calc_ph_energy(data, n=36):
    ps = data[:, 7:113:3]
    xs = data[:, 115::3]
    #print("ps size", np.shape(ps))
    #print("xs size", np.size(xs))
    omega_c = np.zeros(np.shape(ps))
    for idx in range(n):
        omega_c[:, idx] = (2320**2 + (50*(idx+1))**2)**0.5 / 219474.63 # cm-1 to au
    #print("omega_c is", omega_c[0,:]*219474.63)
    es = 0.5 * ps**2 + 0.5*omega_c**2 * xs**2 - 0.00095003725571098 # minus kT
    return es


def obtain_ph_dynamics(path, pattern="simu_*.out", n=36):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    # calculate energy
    es = calc_ph_energy(data=data, n=n)
    for filename in filenames[1:]:
        data = np.loadtxt(filename)
        es += calc_ph_energy(data=data, n=n)
    es /= float(len(filenames))
    t = data[:,1]
    phs = []
    for idx in range(n):
        phs.append(es[:,idx])
    ts = [t]*len(phs)
    return ts, phs

def plot_ph_dynamics(n=36):

    ts, phs = obtain_ph_dynamics(path="../co2_liquid_phase/E0_5e-5_excited2/", n=n)
    
    # calculate the excited and unexcited photon energy
    e_excited, e_rest = np.zeros(ts[0].size), np.zeros(ts[0].size)
    for idx in range(n):
        if idx == 11:
            e_excited = phs[idx]
        else:
            e_rest += phs[idx]
    # calculate photon frequency for all modes
    omega_lst = np.array([(2320**2 + (50*(idx+1))**2)**0.5 for idx in range(n)])

    axes = clp.initialize(3, 1, width=4.3, height=4.3*0.618*3, LaTeX=True, fontsize=12,
                          labelthem=True, labelthemPosition=[-0.06, 1.05], labelsize=13,
                          commonY=[-0.03, 0.5, r"$E_{\rm ph} - k_{\rm B}T$ [$\hbar\omega_{\perp}$]"])

    x1s = ts[0:2]
    e_ph_unit = 2320.0 / 219474.63
    y1s = [e_excited / e_ph_unit, e_rest / e_ph_unit]
    colors = [clp.red, clp.sky_blue]
    labels = [r"$\omega_{\rm c} = 2396$ cm$^{-1}$", r"rest photon modes"]
    clp.plotone(x1s, y1s, axes[0], colors=colors, labels=labels, lw=1.0, 
                xlim=[0,20], ylim=[0, 180], 
                )

    # Fig b
    x2s = ts[5:11]
    y2s = phs[5:11] 
    y2s = [y/e_ph_unit for y in y2s]
    labels = [r"$%d$ cm$^{-1}$" %omega for omega in omega_lst[5:11]]
    clp.plotone(x2s, y2s, axes[1], lw=1.0, labels=labels, ylim=[0,5],
                xlim=[0,20], legendFontSize=8, 
                )
    colormap = plt.cm.hot
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(axes[1].lines))]
    leg = axes[1].get_legend()
    for i,j in enumerate(axes[1].lines):
        j.set_color(colors[i])
        leg.legendHandles[i].set_color(colors[i])
    
    # Fig c
    x3s = ts[12:18]
    y3s = phs[12:18] 
    y3s = [y/e_ph_unit for y in y3s]
    labels = [r"$%d$ cm$^{-1}$" %omega for omega in omega_lst[12:18]]
    clp.plotone(x3s, y3s, axes[2], lw=1.0, labels=labels, ylim=[0,5], legendFontSize=8,
                xlim=[0,20], xlabel="time [ps]", 
                )
    colormap = plt.cm.hot
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(axes[2].lines))][::-1]
    leg = axes[2].get_legend()
    for i,j in enumerate(axes[2].lines):
        j.set_color(colors[i])
        leg.legendHandles[i].set_color(colors[i])

    clp.adjust(savefile="ph_dynamics.pdf")


if __name__ == "__main__":
    plot_ph_dynamics()
