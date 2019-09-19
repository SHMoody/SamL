import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pyvtk import *


class magnetic_structure:

    def __init__(self):
        self.x1 = np.array([1, 0, 0])
        self.y1 = np.array([0, 1, 0])
        self.z1 = np.array([0, 0, 1])

    def cubespace(self, maxht):
        rvals = np.meshgrid(np.arange(maxht),
                            np.arange(maxht),
                            np.arange(maxht))
        self.rvals = np.reshape(rvals, (3, maxht**3), order='C')

    def unit_cell(self, filename, nx, ny, nz):
        unit_cell = np.transpose(
            np.loadtxt(filename, usecols=(1, 2, 3),
                       skiprows=2, dtype=np.single) / 8.92500)

        r = np.meshgrid(np.arange(nx, dtype=np.single),
                        np.arange(ny, dtype=np.single),
                        np.arange(nz, dtype=np.single))
        r = np.reshape(r, (3, nx*ny*nz), order='F')
        lattice = np.tile(unit_cell, nx*ny*nz)
        r = np.repeat(r, np.shape(unit_cell)[1], axis=1)
        self.lat = lattice + r

    def helix(self, r, tau):
        dot = np.dot(tau, r)
        dot1 = np.dot(tau/np.linalg.norm(tau), self.y1)
        dot2 = np.dot(tau/np.linalg.norm(tau), self.x1)
        self.m = np.array([dot1*np.cos(dot),
                           dot2*np.cos(dot),
                           np.sin(dot)])

    def REXS(self, q, lattice, mags, F0, F1, F2, theta, e_in, e_out):
        edote = np.dot(e_in, e_out)
        ecroe = np.cross(e_out, e_in)
        eoutm = np.tensordot(e_in, mags, axes=1)
        ein_m = np.tensordot(e_out, mags, axes=1)
        return np.abs(np.sum(np.exp(1j*np.dot(np.transpose(lattice),
                                              2*np.pi*q)) *
                             (edote*F0 - 1j*np.tensordot(ecroe, mags, axes=1) * F1
                              - eoutm*ein_m*F2)))

    def makevecvtk(self, fname='vectors'):

        points = np.transpose(self.lat)
        magvec = np.transpose(self.m)
        vtk = VtkData(UnstructuredGrid(points,),
                      PointData(Vectors(magvec),
                                Scalars(np.sum(np.absolute(magvec), axis=1))),
                      'test'
                      )
        vtk.tofile(fname)


if __name__ == "__main__":
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.integrate import trapz

    m = magnetic_structure()
    m.unit_cell('Cu2OSeO3.txt', 100, 100, 1)
    m.helix(m.lat, np.array([1./138.547*np.pi, 0, 0]))
    xs = np.linspace(0.95, 1.05, 1000)

    I = np.zeros(shape=(np.shape(xs)))
    for j in range(len(xs)):
        print j
        q = np.array([xs[j], 0, 0])
        I[j] = m.REXS(q, m.lat, m.m, 1, 1, 1, 45*np.pi/180.,
                      np.array([0, 0, 1]), np.array([0, 1, 0]))

    print trapz(I, xs)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xs, I)
    plt.show()
