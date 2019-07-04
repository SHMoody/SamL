import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def Theta_dw(x, y, R, w):
    r = np.sqrt(x**2+y**2)
    return 2*np.arctan(np.sinh(R/w)/np.sinh(r/w))


def m(xs, ys, nu=1, gamma=90.*np.pi/180):
    r = np.sqrt((xs)**2+((ys)**2))
    phi = np.arctan2(ys, xs)
    return (np.array([np.sin(Theta_dw(r))*np.cos(nu*phi + gamma),
                      np.sin(Theta_dw(r))*np.sin(nu*phi+gamma), -1.+np.cos(Theta_dw(r))]))
