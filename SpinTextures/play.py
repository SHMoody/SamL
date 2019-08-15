import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft as fft


def f_REXS_4(F0, F1, F2, z_n, th, r_n):

    a = F0+F2*z2**2
    b = -1j*F1*(z1*np.cos(th)+z3*np.sin(th))-F2 * \
        z2*(z1*np.sin(th)-z3*np.cos(th))
    return np.array([a, b])


def z_n_helix(tau, r_n):
    dot = np.dot(r_n, tau)
    return np.array([np.cos(dot),
                     np.sin(dot),
                     0.*dot])


t = np.array([0, 0, 1])
rs = np.array([0, 0, 4.2])
m = z_n_helix(t, rs)
print m
