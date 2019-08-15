from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pyvtk import *
import sys
sys.path = ['..']+sys.path


lam_h = 7.
chi = 0. * np.pi/180
p1 = np.sqrt(3) / 2


class spintexture2D():

    def __makebasis(self, vector):
        n1 = vector / np.linalg.norm(vector)
        n2 = np.array([-n1[1]+1e-16,  n1[0]+1e-16])
        n2 -= n2.dot(n1)*n1
        n2 = n2 / np.linalg.norm(n2)
        basis = np.array([n1, n2])
        return basis

    def __rotmatrix(self, vector, theta):
        t = theta*np.pi/180.
        rot = np.array([[np.cos(t), -np.sin(t)],
                        [np.sin(t), np.cos(t)]])
        return np.matmul(rot, vector)

    def makespace(self, xmax, ymax, spacing):
        xvals = np.arange(-xmax, xmax+spacing, spacing)
        yvals = np.arange(-ymax, ymax+spacing, spacing)
        y, x = np.meshgrid(yvals, xvals)
        r = np.array(list(x.flatten()) +
                     list(y.flatten()))
        R = r.reshape((len(x.flatten()), 2), order='F')
        self.stepsize = spacing
        self.xmax = xmax
        self.ymax = ymax
        self.xvals = xvals
        self.yvals = yvals
        self.rvals = R
        self.y = y
        self.x = x

    def makehelix(self, lam_h, chi, q_h=[1, 0], k=0):
        chi = chi * np.pi/180.
        basis = self.__makebasis(q_h)
        rT = np.transpose(self.rvals)
        q_h = 2*np.pi/lam_h * basis[0]
        qr = np.zeros(np.shape(rT))
        m = np.zeros((3, np.shape(rT)[1]))
        M = np.zeros(np.shape(m))
        qr[0] = q_h[0]*rT[0]
        qr[1] = q_h[1]*rT[1]
        qr = np.sum(np.transpose(qr), axis=1)
        m[0] = np.cos(chi)*np.sin(qr + k)
        m[1] = np.sin(chi)*np.sin(qr + k)
        m[2] = np.cos(qr + k)
        M[0] = m[0]*basis[0][0]+m[1]*basis[1][0]
        M[1] = m[0]*basis[0][1]+m[1]*basis[1][1]
        M[2] = m[2]
        M = np.transpose(M)
        self.m = M / np.linalg.norm(M, axis=1)[:, np.newaxis]

    def makecone(self, lam_h, chi, tilt, q_h=[1, 0], k=0):
        chi = chi * np.pi/180.
        tilt = float(tilt)*np.pi/180.
        basis = self.__makebasis(q_h)
        rT = np.transpose(self.rvals)
        q_h = 2*np.pi/lam_h * basis[0]
        qr = np.zeros(np.shape(rT))
        m = np.zeros((3, np.shape(rT)[1]))
        M = np.zeros(np.shape(m))
        qr[0] = q_h[0]*rT[0]
        qr[1] = q_h[1]*rT[1]
        qr = np.sum(np.transpose(qr), axis=1)
        m[0] = np.cos(chi)*np.sin(qr + k) + np.sin(tilt)
        m[1] = np.sin(chi)*np.sin(qr + k)
        m[2] = np.cos(qr + k)
        M[0] = m[0]*basis[0][0]+m[1]*basis[1][0]
        M[1] = m[0]*basis[0][1]+m[1]*basis[1][1]
        M[2] = m[2]
        M = np.transpose(M)
        self.m = M / np.linalg.norm(M, axis=1)[:, np.newaxis]

    def makeskyrmion(self, lam_h, chi, q_h=[1, 0], k=0):
        chi = chi * np.pi/180.
        rT = np.transpose(self.rvals)
        m = np.zeros((3, np.shape(rT)[1]))
        M = np.zeros(np.shape(m))

        basis = self.__makebasis(q_h)
        q_h = 2*np.pi/lam_h * basis[0]
        for i in range(3):
            qr = np.zeros(np.shape(rT))
            qr[0] = q_h[0]*rT[0]
            qr[1] = q_h[1]*rT[1]
            qr = np.sum(np.transpose(qr), axis=1)
            m[0] = np.cos(chi)*np.sin(qr + k)
            m[1] = np.sin(chi)*np.sin(qr + k)
            m[2] = np.cos(qr + k)
            M[0] += m[0]*basis[0][0]+m[1]*basis[1][0]
            M[1] += m[0]*basis[0][1]+m[1]*basis[1][1]
            M[2] += m[2]

            q_h = self.__rotmatrix(q_h, 60)
            basis = self.__makebasis(q_h)

        M = np.transpose(M)
        self.m = M / np.linalg.norm(M, axis=1)[:, np.newaxis]

    def plotbasic(self):
        m = np.transpose(self.m)
        r = np.transpose(self.rvals)
        colors = [plt.cm.viridis(x) for x in 0.5+m[2]/(2*np.amax(m[2]))]
        fig = plt.figure()
        ax = fig.add_axes([0, 0, 1, 1])

        ax.quiver(r[0], r[1], m[0], m[1],
                  color=colors, pivot='middle')
        ax.set_aspect('equal')
        plt.show()

    def plotshilei(self):  # More thought needed here.
        m = np.transpose(self.m)
        r = np.transpose(self.rvals)
        colors = [plt.cm.viridis(x) for x in 0.5+m[2]/(2*np.amax(m[2]))]
        fig = plt.figure()
        ax = fig.add_axes([0, 0, 1, 1])
        Z = np.transpose(np.transpose(colors)[: -1])
        im = Z.reshape((len(self.xvals), len(self.yvals), 3))
        ax.quiver((-r[0] + self.xmax)/self.stepsize, (r[1]+self.ymax) /
                  self.stepsize, m[0], m[1])
        ax.imshow(np.rot90(im), interpolation='bilinear',
                  alpha=0.9, origin='upper')
        ax.axis('off')
        plt.show()


class spintexture3D():

    def __makebasis(self, vector):
        # need to make my own basis here
        n1 = vector / np.linalg.norm(vector)
        n2 = np.array([-n1[1]+1e-16,  n1[0]+1e-16, 0])
        n2 -= n2.dot(n1)*n1
        n2 = n2 / np.linalg.norm(n2)
        n3 = np.cross(n1, n2)
        basis = np.array([n1, n2, n3])
        return basis

    def __rotmatrix(self, u, theta):
        t = theta*np.pi/180.
        ux = u[0]
        uy = u[1]
        uz = u[2]
        return np.array([
            [np.cos(t)+ux**2*(1-np.cos(t))], [ux *
                                              uy*(1-np.cos(t))-uz*np.sin(t)],
            [ux*uz*(1-np.cos(t))+uy*np.sin(t)],
            [uy*uz*(1-np.cos(t))+uz*np.sin(t)
             ], [np.cos(t)+uy**2*(1-np.cos(t))],
            [uy*uz*(1-np.cos(t))-ux*np.sin(t)],
            [uz*ux*(1-np.cos(t))-uy*np.sin(t)], [uz *
                                                 uy*(1-np.cos(t))+ux*np.sin(t)],
            [np.cos(t)+uz**2*(1-np.cos(t))]]).reshape(3, 3)

    def makespace(self, xmax, ymax, zmax, stepsize):
        xvals = np.arange(-xmax, xmax+stepsize, stepsize)
        yvals = np.arange(-ymax, ymax+stepsize, stepsize)
        zvals = np.arange(-zmax, zmax+stepsize, stepsize)
        y, x, z = np.meshgrid(yvals, xvals, zvals)
        r = np.array(list(x.flatten()) +
                     list(y.flatten()) + list(z.flatten()))
        R = r.reshape((len(x.flatten()), 3), order='F')
        self.xmax = xmax
        self.ymax = ymax
        self.zmax = zmax
        self.xvals = xvals
        self.yvals = yvals
        self.zvals = zvals
        self.rvals = R
        self.y = y
        self.x = x
        self.z = z
        self.shape = 'rect'

    def makesphere(self, thangle, nthangle,
                   r=np.sqrt(3.)/(2.*np.pi)):
        npoints = np.sum(nthangle)
        self.rvals = np.zeros((npoints, 3))
        self.y = np.zeros(npoints)
        self.x = np.zeros(npoints)
        self.z = np.zeros(npoints)
        thangle = np.pi*np.array(thangle)/180.
        x = 0
        for i in range(len(thangle)):
            phi = np.linspace(0, 2*np.pi, nthangle[i]+1)[:-1]
            for j in phi:
                self.rvals[x] = np.array([r*thangle[i]*np.cos(j),
                                          r*thangle[i]*np.sin(j), 0])
                self.x[x] = r*np.cos(j)*np.sin(thangle[i])
                self.y[x] = r*np.sin(j)*np.sin(thangle[i])
                self.z[x] = r*np.cos(thangle[i])
                x = x+1
        self.shape = 'circ'

    def makehelix(self, lam_h, chi, q_vector=np.array([1, 0, 0]), k=0):
        chi = chi * np.pi/180.
        basis = self.__makebasis(q_vector)
        rT = np.transpose(self.rvals)
        q_h = 2*np.pi/lam_h * basis[0]
        qr = np.zeros(np.shape(rT))
        m = np.zeros(np.shape(rT))
        M = np.zeros(np.shape(rT))
        qr[0] = q_h[0]*rT[0]
        qr[1] = q_h[1]*rT[1]
        qr[2] = q_h[2]*rT[2]
        qr = np.sum(np.transpose(qr), axis=1)
        m[0] = np.cos(chi)*np.sin(qr + k)
        m[1] = np.sin(chi)*np.sin(qr + k)
        m[2] = np.cos(qr + k)
        M[0] = m[0]*basis[0][0]+m[1]*basis[1][0]+m[2]*basis[2][0]
        M[1] = m[0]*basis[0][1]+m[1]*basis[1][1]+m[2]*basis[2][1]
        M[2] = m[0]*basis[0][2]+m[1]*basis[1][2]+m[2]*basis[2][2]
        M = np.transpose(M)
        self.m = M / np.linalg.norm(M, axis=1)[:, np.newaxis]

    def makeskyrmion(self, lam_h, chi, plane=np.array([0, 0, 1]),
                     q_vector=np.array([1, 0, 0]), k=0, N=1):
        print plane
        chi = chi * np.pi/180.
        R = self.__rotmatrix(plane, 60)
        rT = np.transpose(self.rvals)
        M = np.zeros(np.shape(rT))
        q_vectors = np.array([q_vector, np.matmul(R, q_vector),
                              np.matmul(R, np.matmul(R, q_vector))])
        for i in range(3):
            basis = self.__makebasis(q_vectors[i])
            q_h = 2*np.pi/lam_h * basis[0]
            qr = np.zeros(np.shape(rT))
            m = np.zeros(np.shape(rT))
            qr[0] = q_h[0]*rT[0]
            qr[1] = q_h[1]*rT[1]
            qr[2] = q_h[2]*rT[2]
            qr = np.sum(np.transpose(qr), axis=1)
            m[0] = np.cos(chi)*np.sin((qr + k))
            m[1] = np.sin(chi)*np.sin((qr + k))
            m[2] = np.cos(qr + k)
            M[0] += m[0]*basis[0][0]+m[1]*basis[1][0]+m[2]*basis[2][0]
            M[1] += m[0]*basis[0][1]+m[1]*basis[1][1]+m[2]*basis[2][1]
            M[2] += m[0]*basis[0][2]+m[1]*basis[1][2]+m[2]*basis[2][2]
        M = np.transpose(M)
        self.m = M / np.linalg.norm(M, axis=1)[:, np.newaxis]

    def makecone(self, lam_h, chi, gam, q_vector=np.array([1, 0, 0]), k=0):
        chi = chi * np.pi/180.
        gam = gam * np.pi/180.
        basis = self.__makebasis(q_vector)
        rT = np.transpose(self.rvals)
        q_h = 2*np.pi/lam_h * basis[0]
        qr = np.zeros(np.shape(rT))
        m = np.zeros(np.shape(rT))
        M = np.zeros(np.shape(rT))
        qr[0] = q_h[0]*rT[0]
        qr[1] = q_h[1]*rT[1]
        qr[2] = q_h[2]*rT[2]
        qr = np.sum(np.transpose(qr), axis=1)
        m[0] = np.cos(chi)*np.sin(qr + k)+np.sin(gam)
        m[1] = np.sin(chi)*np.sin(qr + k)
        m[2] = np.cos(qr + k)
        M[0] = m[0]*basis[0][0]+m[1]*basis[1][0]+m[2]*basis[2][0]
        M[1] = m[0]*basis[0][1]+m[1]*basis[1][1]+m[2]*basis[2][1]
        M[2] = m[0]*basis[0][2]+m[1]*basis[1][2]+m[2]*basis[2][2]
        M = np.transpose(M)
        self.m = M / np.linalg.norm(M, axis=1)[:, np.newaxis]

    def makevecvtk(self, fname='vectors'):
        x = int(np.cbrt(len(self.rvals)))
        pos = [0, x-1, x**2 - 1, x**2 - x, x**3 - x**2,  # This is broken
               x**3 - x**2 + x - 1,  x**3 - 1, x**3 - x]
        if self.shape == 'circ':
            points = np.zeros(np.shape(np.transpose(self.rvals)))
            points[0] = self.x
            points[1] = self.y
            points[2] = self.z
            points = np.transpose(points)
        else:
            points = self.rvals
        print self.m
        print np.shape(self.m)
        vtk = VtkData(UnstructuredGrid(points,),
                      PointData(Vectors(self.m),
                                Scalars(np.sum(np.absolute(self.m), axis=1))),
                      # Scalars(range(len(self.rvals)))),
                      'test'
                      )
        vtk.tofile(fname)

    def makefouvtk(self, fname='fourier'):
        x = int(np.cbrt(len(self.rvals)))
        pos = [0, x-1, x**2 - 1, x**2 - x, x**3 - x**2,
               x**3 - x**2 + x - 1,  x**3 - 1, x**3 - x]
        # Need to fix the hexadron
        vtk = VtkData(UnstructuredGrid(self.rvals,
                                       hexahedron=pos,),
                      PointData(Vectors(np.absolute(self.f)),
                                Scalars(range(len(self.rvals)))),
                      'test'
                      )
        vtk.tofile('fourier')

    def plot_quiver(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        m = np.transpose(self.m)
        colors = [plt.cm.viridis(x)
                  for x in 1-(0.5+np.tensordot(self.m, np.array([0, 0, 1]),
                                               axes=1)/2.)]
        plt.xlabel('x')
        plt.ylabel('y')
        ax.quiver((self.x).ravel(), (self.y).ravel(), (self.z).ravel(),
                  m[0].ravel(), m[1].ravel(), m[2].ravel(), color=colors,
                  pivot='middle', arrow_length_ratio=0.5, length=0.05)
        plt.show()

    def fourier(self):  # Need to think of a better way of doing this - this is bad
        self.f = np.fft.fftshift(np.fft.fftn(self.m))

    def smooth_fourier(self, cutoff, size):
        print 'i do nothing'
        # cutoff = float(cutoff)
        # F = np.absolute(self.f)
        # x = np.argmax(F.ravel())
        # print np.argwhere((F > x/cutoff))
        #
        # rvals = np.zeros()


if __name__ == "__main__":
    q = spintexture3D()
    # , r=1/np.pi)
    # q.makesphere([0, 30, 60, 90, 120, 150, 180],
    #             [1, 12, 18, 24, 18, 12, 1])
    q.makespace(120, 120, 0, 12)
    q.makeskyrmion(1, 90)
    # q.makecone(70., 90, 0, q_vector=[1, 0, 0])  # make helix *2/rt(3)
    q.makevecvtk()
