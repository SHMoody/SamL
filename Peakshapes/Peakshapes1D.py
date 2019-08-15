import numpy as np


def lorentzian(x, x0, fwhm, A, c):
    l1 = fwhm / 2.
    l2 = (x-x0)**2
    return A*(1/np.pi)*(l1/(l1**2+l2)) + c


def gaussian(x, x0, fwhm, A, c):
    g1 = 2.*np.sqrt(np.log(2.) / np.pi)
    g2 = 4.*np.log(2.)
    return A*(g1/fwhm)*np.exp(-g2*(x-x0)**2/(fwhm**2))+c


def pseudovoight(x, x0, fwhm, n, A, c):
    return A*(n*lorentzian(x, x0, fwhm, 1, 0)
              + (1-n)*gaussian(x, x0, fwhm, 1, 0)) + c


def pearsonVII(x, x0, fwhm, m, A, c):  # Might be broken - not unit area
    l1 = fwhm / 2.
    l2 = (x-x0)**2
    return A*((1/np.pi)*(l1/(l1**2+l2)))**m + c


def asymgaussian(x, x0, sigL, sigR, A, c):  # not unit area
    mp = np.argmax(x > x0)
    xleft = x[:mp-1]
    xright = x[mp-1:]
    return A*np.concatenate((
        np.exp(-(xleft-x0)**2/sigL**2), np.exp(-(xright-x0)**2/sigR**2)), axis=0)+c


def asymlorentzian(x, x0, yL, yR, A, c):  # not unit area
    mp = np.argmax(x > x0)
    xleft = x[:mp-1]
    xright = x[mp-1:]
    return A*np.concatenate((
        (yL**2) / ((xleft-x0)**2+yL**2), (yR**2) / ((xright-x0)**2+yR**2)), axis=0)+c


def sinc(x, x0, q, A, c):
    return A*(np.sin(q*((x-x0)))/(q*(x-x0)))**2+c


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    xvals = np.linspace(-5., 5., 1000.)

    plt.figure()
    for i in range(1, 10):
        pvals = sinc(xvals, 1., 1/float(i), 1., 0.)**2
        plt.plot(xvals, pvals, label=str(i))
        plt.legend()
    plt.show()
