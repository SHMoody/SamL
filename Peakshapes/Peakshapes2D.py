import numpy as np


def lorentzian2D((x, y), x0, y0, fwhmx, fwhmy, A, c):
    l1x = 0.5*fwhmx
    l1y = 0.5*fwhmy
    lx = (x-x0)**2
    ly = (y-y0)**2
    return np.array((A/np.pi**2)*((l1x*l1y) /
                                  (l1x**2*l1y**2 + ly*l1x**2 + lx*l1y**2 + ly*lx))+c).ravel()


def adv_lorentzian2D((x, y), x0, y0, fwhmx, fwhmy, th, A, c):
    x = x-x0
    y = y - y0
    Wx = fwhmx / (2*(2**(2./3.)-1))
    Wy = fwhmy / (2*(2**(2./3.)-1))
    v1 = Wx**2 * np.cos(th)**2 + Wy**2 * np.sin(th)**2
    v2 = Wx**2 * np.cos(th)*np.sin(th) - Wy**2*np.cos(th)*np.sin(th)
    v3 = Wx**2 * np.sin(th)**2 + Wy**2 * np.cos(th)**2
    V = np.array([[v1, v2], [v2, v3]])
    V_inv = np.linalg.inv(V)
    return (A*(1./(2*np.pi*np.linalg.det(V)**0.5) * (1 / (1 + V_inv[0][0]*x*x + 2*V_inv[0][1]*x*y + V_inv[1][1]*y*y)**1.5)) + c).ravel()


def gaussian2D((x, y), x0, y0, fwhmx, fwhmy, A, c):
    x0 = float(x0)
    y0 = float(y0)
    g1 = 2.*np.sqrt(np.log(2.) / np.pi)
    g2 = 4.*np.log(2.)
    return (A*(g1/fwhmx)*(g1/fwhmy)*np.exp(
        -g2*(x-x0)**2/(fwhmx**2) - g2*(y-y0)**2/(fwhmy**2)) + c).ravel()


def adv_gaussian2D((x, y), x0, y0, fwhmx, fwhmy, th, A, c):
    x = x-x0
    y = y - y0
    Wx = fwhmx / (2*np.sqrt(2*np.log(2)))
    Wy = fwhmy / (2*np.sqrt(2*np.log(2)))
    v1 = Wx**2 * np.cos(th)**2 + Wy**2 * np.sin(th)**2
    v2 = Wx**2 * np.cos(th)*np.sin(th) - Wy**2*np.cos(th)*np.sin(th)
    v3 = Wx**2 * np.sin(th)**2 + Wy**2 * np.cos(th)**2
    V = np.array([[v1, v2], [v2, v3]])
    V_inv = np.linalg.inv(V)
    expo = -0.5*(V_inv[0][0]*x*x + 2*V_inv[0][1]*x*y + V_inv[1][1]*y*y)
    return (A*((1. / np.sqrt(4*(np.pi**2)*np.linalg.det(V))) * np.exp(expo))+c).ravel()


def pseudovoight2D((x, y), x0, y0, fwhmx, fwhmy, n, A, c):
    return A*(n*lorentzian2D((x, y), x0, y0, fwhmx, fwhmy, 1, 0)
              + (1-n)*gaussian2D((x, y), x0, y0, fwhmx, fwhmy, 1, 0)) + c


def adv_pseudovoight2D((x, y), x0, y0, fwhmx, fwhmy, n, th, A, c):
    return A*(n*adv_lorentzian2D((x, y), x0, y0, fwhmx, fwhmy, th, 1, 0)
              + (1-n)*adv_gaussian2D((x, y), x0, y0, fwhmx, fwhmy, th, 1, 0)) + c


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import time
    from scipy.integrate import simps
    from scipy.integrate import dblquad

    xvals = np.linspace(-3, 3, 1000)
    yvals = np.linspace(-3, 3, 1000)
    x, y = np.meshgrid(xvals, yvals)
    #z = np.array(adv_lorenzian2D((x, y), 0., 0., 1. ,1., 0.,  1., 0.))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    z = np.array(adv_lorentzian2D((x, y), 0., 0., 0.2, 2., 10.,  1., 0.))
    z2 = np.array(adv_gaussian2D((x, y), 0., 0., 0.2, 2., 10., 1., 0.))
    z3 = np.array(adv_pseudovoight2D(
        (x, y), 0., 0., 0.2, 2., 0.5, 10., 1., 0.))

    #z2 = np.array(adv_lorentzian2D((x, y), 0., 0., 1. ,1., 0., 1., 0.))
    #x0 = adv_gaussian2D((0, 0), 0., 0., 1. ,1., 10., 1., 0.)
    #x1 =  adv_pseudovoight2D((0, 0.5), 0., 0., 1. ,1., 10., 1., 0.)
    # print x1/x0 , x0, x1

   # x0 = adv_lorentzian2D((0, 0), 0., 0., 1. ,1., 10., 1., 0.)
   # x1 =  adv_lorentzian2D((0, 0.5), 0., 0., 1. ,1., 10., 1., 0.)
   # print x1/x0 , x0, x1

    # print adv_lorentzian2D((0, 0.5), 0., 0., 1., 1., 0., 1., 0.)
    # print adv_lorentzian2D((0., 0.), 0., 0., 1., 1., 0., 1., 0.)
    ax.plot_surface(x, y, z.reshape(len(yvals), len(xvals)),
                    cmap=cm.viridis, alpha=0.6)
    ax.plot_surface(x, y, z2.reshape(len(yvals), len(xvals)),
                    cmap=cm.plasma, alpha=0.6)
    # ax.plot_surface(x, y, z3.reshape(len(yvals),len(xvals)),
    #                cmap = cm.winter, alpha = 0.5)
    # print dblquad(adv_lorenzian2D_int, -3, 3, lambda x:-3, lambda x:3,
    # args = (0., 0., 1. ,1., 0,  2., 0.))
    print simps(simps(z.reshape(len(yvals), len(xvals)), yvals), xvals)
    print simps(simps(z2.reshape(len(yvals), len(xvals)), yvals), xvals)
    print simps(simps(z3.reshape(len(yvals), len(xvals)), yvals), xvals)

    plt.show()
