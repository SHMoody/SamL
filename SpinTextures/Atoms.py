from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


def M_helix(r, tau):
    dot = np.dot(r, tau)
    return np.array([np.cos(dot),
                     np.sin(dot),
                     dot*0])


unit_cell = np.transpose(
    np.loadtxt('Cu2OSeO3.txt', usecols=(1, 2, 3), skiprows=0) / 8.92500)


for j in np.arange(5000, 25001, 5000):
    print j
    x1, y1, z1 = 1, 1, j

    num_cells = 3
    y, z, x = np.meshgrid(np.arange(y1), np.arange(
        z1), np.arange(x1))
    x = x.ravel()
    y = y.ravel()
    z = z.ravel()
    for i in range(x1*y1*z1):
        print i, x1*y1*z1, j
        if i == 0:
            lattice = unit_cell
        else:
            lattice = np.append(lattice, np.array([unit_cell[0]+x[i],
                                                   unit_cell[1]+y[i],
                                                   unit_cell[2]+z[i]]), axis=1)

    moments = M_helix(np.transpose(lattice), np.array([0, 0, 2*np.pi/69.468]))

    q = np.linspace(4.01, 4.02, 1500)
    I1 = []
    I2 = []
    I3 = []
    for i in q:
        print i, j
        #   I1.append(np.abs(np.sum(
        #       np.exp(1j*np.dot(np.transpose(lattice),
        #                        2*np.pi*np.array([0, 0, i])))
        #       * moments[0]**2)))  # moments[1]**2)))  # 1j*moments[0])))  # (1j*moments[0]+moments[1]*moments[1]))))

        I2.append(np.abs(np.sum(
            np.exp(1j*np.dot(np.transpose(lattice),
                             2*np.pi*np.array([0, 0, i])))
            * 1j*moments[0])))

        #    I3.append(np.abs(np.sum(
        #        np.exp(1j*np.dot(np.transpose(lattice),
        #                         2*np.pi*np.array([0, 0, i])))
        #        * moments[1]*moments[0])))

    # print np.shape(np.array(
    #    [q, np.array(I1)+np.array(I2)+np.array(I3)]))
    np.savetxt(str(j)+'_datafromatoms.txt', np.transpose(np.array(
        [q, np.array(I2)])))

    #plt.plot(q, I1, C='r')
    #plt.plot(q, I2, c='b')
    #plt.plot(q, I3, c='g')
    # plt.show()
