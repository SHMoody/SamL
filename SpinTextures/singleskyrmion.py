import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def Theta_dw(r, R, w):
    return 2*np.arctan(np.sinh(R/w)/np.sinh(r/w))


def m(xs, ys, nu=1, gamma=90.*np.pi/180):
    r = np.sqrt((xs)**2+((ys)**2))
    phi = np.arctan2(ys, xs)
    return (np.array([np.sin(Theta_dw(r, 10, 5))*np.cos(nu*phi + gamma),
                      np.sin(Theta_dw(r, 10, 5))*np.sin(nu*phi+gamma),
                      -1.+np.cos(Theta_dw(r, 10, 5))]))


def lattice(xs, ys, x, y):
    lat = np.zeros(np.shape(xs))

    for i in range(len(x)):
        for j in range(len(y)):
            lat[x[i]][y[j]] = 1
    return lat


if __name__ == "__main__":
    xs, ys = np.meshgrid(np.linspace(-500, 500, 201),
                         np.linspace(-500, 500, 201))

    m = m(xs, ys)
    la = lattice(xs, ys, np.linspace(0, 200, 55, dtype=int),
                 np.linspace(0, 200, 55, dtype=int))

    plt.figure()
    mx = np.abs(np.fft.ifft2((np.fft.fft2(m[0]))*(np.fft.fft2(la))))
    my = np.abs(np.fft.ifft2((np.fft.fft2(m[1]))*(np.fft.fft2(la))))
    mz = np.abs(np.fft.ifft2((np.fft.fft2(m[2]))*(np.fft.fft2(la))))

    plt.figure()
    plt.imshow(mz)
    plt.show()

    #fig, ax = plt.subplots(3, 2)
    # ax[0][0].imshow(m[0])
    # ax[1][0].imshow(m[1])
    # ax[2][0].imshow(m[2])
#
    # ax[0][1].imshow(np.log(np.abs(np.fft.fftshift(np.fft.fft2(m[0])))))
    # ax[1][1].imshow(np.log(np.abs(np.fft.fftshift(np.fft.fft2(m[1])))))
    # ax[2][1].imshow(np.log(np.abs(np.fft.fftshift(np.fft.fft2(m[2])))))
#
    # for i in range(3):
    #    for j in range(2):
    #        ax[i][j].set_xticks([])
    #        ax[i][j].set_yticks([])
    ## plt.quiver(xs, ys, m[0], m[1], m[2], pivot='middle', units='x')
    # plt.tight_layout()
    # plt.show()
#
