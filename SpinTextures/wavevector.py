import numpy as plt
import matplotlib.pyplot as plt


def kp_x(q, theta):
    theta = theta*np.pi/180.
    return np.array([
        q[0]/2. - np.linalg.norm(q)/np.tan(theta) *
        np.sin(np.pi/2. - np.atan(q[0]/q[1])),
        q[1]/2. - np.linalg.norm(q)/np.tan(theta) *
        np.cos(np.pi/2. - np.atan(q[0]/q[1]))])
