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

def smooth(x,window_len=11,window='hamming'):
    """smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal
    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    see also:
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

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

def obtain_ph_dynamics_from_file(path, pattern="simu_*.phe", n=36):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    t = data[:,0]
    ts, phs = [], []
    for idx in range(n):
        phs.append(data[:,idx+1])
        ts.append(t*1e-3)
    return ts, phs


def get_excited_unexcited_energy_traj(path="../co2_liquid_phase_systemsize/excite_ph_size_times_4/", n=36, rescaling=1.0/4.0):

    ts, phs = obtain_ph_dynamics_from_file(path=path, n=n)
    # calculate the excited and unexcited photon energy
    e_excited, e_rest = np.zeros(ts[0].size), np.zeros(ts[0].size)
    for idx in range(n):
        if idx == int(n // 36) * 12 - 1:
            e_excited = phs[idx] * rescaling
        else:
            e_rest += phs[idx] * rescaling

    e_ph_unit = 2320.0 / 219474.63
    ys = [e_excited / e_ph_unit, e_rest / e_ph_unit]

    return ts[0:2], ys

def plot_ph_dynamics(n=36*4):


    axes = clp.initialize(2, 1, width=4.3, height=4.3*0.618*2, LaTeX=True, fontsize=12,
                          labelthem=True, labelthemPosition=[-0.06, 1.05], labelsize=13,
                          commonY=[-0.03, 0.5, r"$(E_{\rm ph} - k_{\rm B}T) / N_{\rm simu}$ [$\hbar\omega_{\perp}$]"])

    Nmol = 36*36

    colors = [clp.red, clp.sky_blue]

    labels_ref = [r"cav-\#1 $\omega_{\rm c} = 2396$ cm$^{-1}$", r"cav-\#1 rest photon modes"]
    x1s, y1s = get_excited_unexcited_energy_traj(path="../co2_liquid_phase_systemsize/excite_ph_size_times_1/", n=36, rescaling=1.0/Nmol)

    clp.plotone(x1s, y1s, axes[0], colors=colors, labels=labels_ref, lw=1.0,
                xlim=[0,20], ylim=[0, 180/Nmol],
                )

    clp.plotone(x1s, y1s, axes[1], colors=colors, labels=labels_ref, lw=1.0,
                xlim=[0, 20], ylim=[0, 180/Nmol],
                )

    x2s, y2s = get_excited_unexcited_energy_traj(path="../co2_liquid_phase_systemsize/excite_ph_size_times_4/",
                                                 n=36 * 4,
                                                 rescaling=1.0 / 4.0 / Nmol)
    labels = [r"cav-\#1.4.4 $\omega_{\rm c} = 2396$ cm$^{-1}$", r"cav-\#1.4.4 rest photon modes"]
    clp.plotone(x2s, y2s, axes[0], colors=colors, labels=labels, lw=1.0,
                xlim=[0, 20], ylim=[0, 180 / Nmol], alpha=0.7,
                )

    x2s, y2s = get_excited_unexcited_energy_traj(path="../co2_liquid_phase_systemsize/excite_ph_size_times_16/", n=36*16,
                                                 rescaling=1.0/16.0/Nmol)
    labels = [r"cav-\#1.16.16 $\omega_{\rm c} = 2396$ cm$^{-1}$", r"cav-\#1.16.16 rest photon modes"]
    clp.plotone(x2s, y2s, axes[0], colors=colors, labels=labels, lw=1.0,
                xlim=[0, 20], ylim=[0, 180/Nmol],  alpha=0.3,
                )

    x3s, y3s = get_excited_unexcited_energy_traj(path="../co2_liquid_phase_systemsize/excite_ph_mode_times_1/", n=36,
                                                 rescaling=1.0/4.0/Nmol)
    labels = [r"cav-\#1.4.1 $\omega_{\rm c} = 2396$ cm$^{-1}$", r"cav-\#1.4.1 rest photon modes"]
    clp.plotone(x3s, y3s, axes[1], colors=colors, labels=labels, lw=1.0,
                xlim=[0, 20], ylim=[0, 180/Nmol], alpha=0.5, xlabel="time [ps]",
                )


    clp.adjust(savefile="ph_dynamics_size_dep.pdf")


if __name__ == "__main__":
    plot_ph_dynamics()
