import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from Peakshapes2D import adv_gaussian2D


class skyrmion2D():

    def makehexspace(self, xpoints, ypoints, skyrmionsize, spacing):
        '''
        Creates a 2d grid with hexagonal references
        '''
        s = 1./spacing
        self.a = skyrmionsize
        self.b = 4*np.pi/(skyrmionsize*np.sqrt(3))

        self.xs, self.ys = np.meshgrid(
            np.linspace(-xpoints*self.a,
                        xpoints*self.a,
                        xpoints/s)/2.,
            np.sqrt(3)*np.linspace(-ypoints*self.a,
                                   ypoints*self.a, ypoints/s)/2.)
        self.rs = np.sqrt(self.ys**2 + self.xs**2)
        temp = np.zeros((spacing, spacing))
        temp[0][0] = 1
        self.lattice = np.tile((
            np.tile(temp, (int(ypoints), int(xpoints)))), (3, 1, 1))

    def theta_DW(self, r):  # allow other functions
        R = 13  # make this more flixible and in units of lattice
        w = 4
        return 2*np.arctan(np.sinh(R/w)/np.sinh(r/w))

    def makeskyrmion(self, nu=1, gamma=90.*np.pi/180):
        phi = nu*np.arctan2(self.ys, self.xs)+gamma
        sins = np.sin(self.theta_DW(self.rs))
        coss = np.cos(self.theta_DW(self.rs))
        m = np.array([sins*np.cos(phi), sins*np.sin(phi), coss])

        self.M = m

    def FTskyrmion(self):
        self.ftM = (np.fft.fft2(self.M))

    def makelattice(self):
        conv = np.abs(np.fft.fft2(self.lattice))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(sk.xs, sk.ys, conv[0])
        plt.xlabel('x')
        plt.show()

        self.lat = np.abs(
            np.fft.ifft2(self.ftM*self.ftlattice))
        print np.shape(self.lat)
        quit()


# print np.transpose(m)
# print np.sum(np.transpose(np.absolute(m)), axis=2)


'''
x = np.arange(0, 60, 0.1*1)
y = np.arange(0, 60, 0.1*0.8660254037844386)

x1, y1 = np.meshgrid(x+0.05, y+0.05)
x2, y2 = np.meshgrid(x, y)
xs = np.concatenate((x1, x2))
ys = np.concatenate((y1, y2))


nx = np.arange(0, 60, 18*1)
ny = np.arange(0, 60, 18*0.8660254037844386)

nx1, ny1 = np.meshgrid(nx+9, ny+9)
nx2, ny2 = np.meshgrid(nx, ny)
nxs = np.concatenate((nx1, nx2))
nys = np.concatenate((ny1, ny2))

zs = makelattice(xs, ys, nxs, nys)
'''


'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xs, ys, zs, cmap=cm.coolwarm,
                rstride=5, cstride=5)
plt.show()
'''

if __name__ == "__main__":
    sk = skyrmion2D()
    sk.makehexspace(20, 20, 70, 2)
    sk.makeskyrmion()
    sk.FTskyrmion()
    sk.makelattice()

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(sk.xs, sk.ys, sk.ftlattice[0])
    # plt.xlabel('x')
    fig, ax = plt.subplots()
    ax.quiver(sk.xs, sk.ys, sk.lat[0], sk.lat[1], sk.lat[2])
    # scale=70)
    ax.set_aspect('equal')

    plt.show()
