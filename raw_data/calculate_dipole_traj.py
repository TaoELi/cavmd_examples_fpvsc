import numpy as np
from scipy import signal
from scipy import fftpack
import glob
#import columnplots as clp
import MDAnalysis as mda
#import columnplots as clp
import sys

'''
Note that IR calculations should be run at NVE ensemble!!!
Outside the cavity, a very long simulation of a single Fe(CO)5 is presented
Dipole autocorrelaction function will be calculated!
Inside the cavity, a 1 ps simulation of 16 Fe(CO)5 will be given
Photon autocorrelation function will be calculated!
'''

'''
Universal parameters
'''

au2fs = 0.02418884254
au2eV = 27.211399
eV2cminv = 8065.540106923572
fsinv2eV = 4.135668 # 1fs-1 to 4.13 eV
fsinv2cminv = fsinv2eV * eV2cminv

'''
Parameters for controlling the Pade approximation
'''
sigma = 5e5
w_step = 1e-5
e_cutoff_ev = 0.5
e_cutoff_au = e_cutoff_ev / au2eV
e_start_ev = 0.01
e_start_au = e_start_ev / au2eV


def fft_pade(time,signal,sigma=sigma,max_len=None,w_min=e_start_au,w_max=e_cutoff_au,w_step=w_step,read_freq=None):
    """ Routine to take the Fourier transform of a time signal using the method
          of Pade approximants.
        Inputs:
          time:      (list or Numpy NDArray) signal sampling times
          signal:    (list or Numpy NDArray)
        Optional Inputs:
          sigma:     (float) signal damp factor, yields peaks with
                       FWHM of 2/sigma
          max_len:   (int) maximum number of points to use in Fourier transform
          w_min:     (float) lower returned frequency bound
          w_max:     (float) upper returned frequency bound
          w_step:    (float) returned frequency bin width
        Returns:
          fsignal:   (complex NDArray) transformed signal
          frequency: (NDArray) transformed signal frequencies
        From: Bruner, Adam, Daniel LaMaster, and Kenneth Lopata. "Accelerated
          broadband spectra using transition signal decomposition and Pade
          approximants." Journal of chemical theory and computation 12.8
          (2016): 3741-3750.
    """

    # center signal about zero
    signal = np.asarray(signal) - signal[0]

    stepsize = time[1] - time[0]

    # Damp the signal with an exponential decay.
    damp = np.exp(-(stepsize*np.arange(len(signal)))/float(sigma))
    signal *= damp

    M = len(signal)
    N = int(np.floor(M / 2))

    # Check signal length, and truncate if too long
    if max_len:
        if M > max_len:
            N = int(np.floor(max_len / 2))

    # G and d are (N-1) x (N-1)
    # d[k] = -signal[N+k] for k in range(1,N)
    d = -signal[N+1:2*N]

    try:
        from scipy.linalg import toeplitz, solve_toeplitz
        # Instead, form G = (c,r) as toeplitz
        #c = signal[N:2*N-1]
        #r = np.hstack((signal[1],signal[N-1:1:-1]))
        b = solve_toeplitz((signal[N:2*N-1],\
            np.hstack((signal[1],signal[N-1:1:-1]))),d,check_finite=False)
    except (ImportError,np.linalg.linalg.LinAlgError) as e:
        # OLD CODE: sometimes more stable
        # G[k,m] = signal[N - m + k] for m,k in range(1,N)
        G = signal[N + np.arange(1,N)[:,None] - np.arange(1,N)]
        b = np.linalg.solve(G,d)

    # Now make b Nx1 where b0 = 1
    b = np.hstack((1,b))

    # b[m]*signal[k-m] for k in range(0,N), for m in range(k)
    a = np.dot(np.tril(toeplitz(signal[0:N])),b)
    p = np.poly1d(np.flip(a))
    q = np.poly1d(np.flip(b))

    if read_freq is None:
        # choose frequencies to evaluate over
        frequency = np.arange(w_min,w_max,w_step)
    else:
        frequency = read_freq

    W = np.exp(-1j*frequency*stepsize)

    fsignal = p(W)/q(W)

    return frequency, fsignal

def fft(x, dtfs):
    # Adding zeros to the end of x
    #N = 1500
    #x = np.pad(x, (0, N), 'constant')
    lineshape = fftpack.dct(x, type=1)
    freq_au = np.linspace(0, 0.5/dtfs * 1e15, len(x))
    # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
    freq_cminverse = freq_au / (100.0 * 299792458.0)
    # Calculate spectra
    #field_description =  freq_au**2
    field_description =  freq_au**2
    spectra = lineshape * field_description
    return freq_cminverse, spectra
    #return freq_cminverse[0:spectra.size//2], spectra[0:spectra.size//2]

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

def auto_correlation_function_simple(x):
    n = x.size
    if n % 2 == 0:
        x_shifted = np.zeros(n*2)
    else:
        x_shifted = np.zeros(n*2-1)
    x_shifted[n//2 : n//2+n] = x
    # Convolute the shifted array with the flipped array, which is equivalent to performing a correlation
    autocorr_full = (signal.fftconvolve(x_shifted, x[::-1], mode='same')[-n:]/ np.arange(n, 0, -1))
    # Truncate the autocorrelation array
    autocorr = autocorr_full[0:n//2]
    return autocorr

def get_dipole_spectrum(xyz_filename, save2file=False, n_co2=36*36, dt_fs=2):
    # load xyz file
    traj = mda.Universe(xyz_filename)
    frames = traj.trajectory
    nframes = len(frames)
    print("In %s nframes = %d" %(xyz_filename, nframes))
    partial_charge = np.array([1.0, -0.5, -0.5] * n_co2)
    dx, dy, dz = [], [], []
    for idx, ts in enumerate(frames):
        current_coord = ts._pos.copy()
        dx.append(np.sum(current_coord[:n_co2*3, 0] * partial_charge))
        dy.append(np.sum(current_coord[:n_co2*3, 1] * partial_charge))
        dz.append(np.sum(current_coord[:n_co2*3, 2] * partial_charge))
    dipole_x = np.array(dx)
    dipole_y = np.array(dy)
    dipole_z = np.array(dz)

    dacf_x = auto_correlation_function_simple(dipole_x)
    dacf_y = auto_correlation_function_simple(dipole_y)
    dacf_z = auto_correlation_function_simple(dipole_z)

    dacf_x_freq, dacf_x_sp = fft(dacf_x, dt_fs)
    dacf_y_freq, dacf_y_sp = fft(dacf_y, dt_fs)
    dacf_z_freq, dacf_z_sp = fft(dacf_z, dt_fs)

    sp_tot = smooth(dacf_x_sp) + smooth(dacf_y_sp) + smooth(dacf_z_sp)
    #sp_tot = dacf_x_sp + dacf_y_sp + dacf_z_sp
    #sp_tot /= np.max(sp_tot)

    if save2file:
        print("save dipole spectrum file to disk")
        data = np.zeros((dacf_x_freq.size, 2))
        data[:,0] = dacf_x_freq
        data[:,1] = sp_tot
        np.savetxt(xyz_filename + ".sp", data)
    return dacf_x_freq, sp_tot

pathname = sys.argv[-1]
xyz_filenames = glob.glob(pathname + "/simu_*.xc.xyz")

xs, ys = [], []
for xyz_filename in xyz_filenames:
    print("Checking", xyz_filename)
    print("Checking", xyz_filename)
    try:
        from pathlib import Path
        path = Path(xyz_filename+".sp")
        if path.is_file():
            print(xyz_filename, "done!")
            continue
        freq, sp = get_dipole_spectrum(xyz_filename=xyz_filename, save2file=True)
        xs.append(freq)
        ys.append(sp)
    except:
        print("Error!!!")
