import numpy as np


def satellite(tau):
    a = np.sqrt(3)/2.
    return np.array([[tau, 0, 0],
                     [0.5*tau, a*tau, 0],
                     [0.5*tau, -a*tau, 0],
                     [-tau, 0, 0],
                     [-0.5*tau, -a*tau, 0],
                     [-0.5*tau, a*tau, 0]])


def d_hexgonal(a, c, h, k, l):
    return 1/np.sqrt(4./3. * (h**2 + h*k + k**2)/a**2+l**2/c**2)


def bragg_angle(lam, d):
    return np.arcsin(lam/(2*d))*180/np.pi


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    sats = (satellite(0.138))
    nsats = (satellite(0.138))
    peaks = np.array([[0, 0, 0],
                      [1, 0, 0],
                      [2, 0, 0],
                      [3, 0, 0],
                      [4, 0, 0],
                      [5, 0, 0],
                      [6, 0, 0],
                      [7, 0, 0],
                      [8, 0, 0],
                      [9, 0, 0]])
    for i in peaks:
        nsats = np.append(nsats, sats+i, axis=0)
    sats = np.transpose(nsats)
    print sats
    plt.figure(figsize=(18, 2.3))
    plt.scatter(sats[0], sats[1])
    plt.scatter(np.transpose(peaks)[0], np.transpose(peaks)[1])
    plt.axis('equal')
    plt.xlabel('h /$r.l.u$')
    plt.ylabel('k /$r.l.u$')
    plt.tight_layout()
    plt.show()
    print bragg_angle(1.5634, d_hexgonal(8.8142, 9.5692, 7, 0, 0))
