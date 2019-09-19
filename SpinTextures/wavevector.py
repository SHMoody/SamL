import numpy as np


def k_vectors(q, theta):
    theta = theta*np.pi/180.
    r = np.array([q[1]/q[0], -1]) * (
        np.tan(theta)*np.linalg.norm(q)/(2*np.sqrt(1+(q[1]/q[0])**2)))

    return np.array([-(q/2.+r), q/2.-r])


if __name__ == "__main__":

    q = np.array([8.716, 0.339])
    kv = k_vectors(q, 44.147)
    print kv[1] - kv[0]
