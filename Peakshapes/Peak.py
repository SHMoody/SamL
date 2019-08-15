import matplotlib.pyplot as plt
import numpy as np
from __builtin__ import any as b_any
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import h5py as h5
import scipy
import matplotlib.patches as patches
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.integrate import trapz
from scipy.integrate import simps
from Peakshapes1D import *
from Peakshapes2D import *


class Peak1D:

    def __init__(self, xdata, ydata, yerrs=0):

        if type(yerrs) == int:
            self.yerrs = np.sqrt(ydata)
        else:
            self.yerrs = yerrs
        self.xdata = xdata
        self.ydata = ydata
        self.center = self.xdata[list(self.ydata).index(max(self.ydata))]
        # finds the xvalue at which y is greatest
        self.maxht = self.ydata.max()*10
        self.minht = self.ydata[0]  # assumes first value is backgroud
        self.gfwhm = self.xdata[7] - self.xdata[0]
        # Assumes ~7 points across peak, xdata evenly spaced
        self.shapes = ['Sinc', 'Gaussian', 'Lorentzian', 'PearsonVII',
                       'Pseudovoight', 'Asymgauss', 'Asymloren',
                       'Multiple']
        self.funcs = [sinc, gaussian, lorentzian, pearsonVII,
                      pseudovoight, asymgaussian, asymlorentzian,
                      self.multipeak]

    def update_init(self, params):

        # changes initial parameters based on previous fits
        self.gx0 = params[0]
        self.gfwx = params[2]
        self.gmh = params[-2]
        self.gbg = params[-1]

    def __argsareeasy(self, *args):

        args = list(args)
        if len(args) == 1:
            args = args[0]
        fakelist = []
        nargs = []
        x = 0
        for i in self.functions:
            x += 1
            if i == 'Gaussian' or i == 'Lorentzian' or i == 'Sinc':
                fakelist.append(3)
            elif i == 'Pseudovoight' or i == 'Asymgauss' or i == 'Asymloren':
                fakelist.append(4)
            else:
                print i, 'does not feature in',  self.shapes[:-1]

        x = 0
        for i in fakelist:
            nargs.append(args[x:x+i]+[0])

            x = x + i
        return nargs + [args[-1]]

    def boundscheck(self):
        if hasattr(self, 'bounds'):
            pass
        else:
            self.makebounds()

    def makebounds(self, bounds=0):
        if bounds == 0:
            if self.shape == 'Gaussian' or self.shape == 'Lorentzian'\
                    or self.shape == 'Sinc':
                self.bounds = ((-np.inf, 0, 0, -np.inf),
                               (np.inf, np.inf, np.inf, np.inf))
            elif self.shape == 'Pseudovoight':
                self.bounds = ((-np.inf, 0, 0, 0, -np.inf),
                               (np.inf, np.inf, 1,  np.inf, np.inf))
            elif self.shape == 'Asymgauss' or self.shape == 'Asymloren':
                self.bounds = ((-np.inf, 0, 0, 0, 0),
                               (np.inf, np.inf, np.inf, np.inf, np.inf))
            elif self.shape == 'Multiple':
                # se np.tile(np.array([-np.inf, np.inf]),
                #              len(self.args))
                self.bounds = tuple(
                    map(tuple, np.tile(np.array([-np.inf, np.inf]),
                                       len(self.args)).reshape(2,
                                                               len(self.args),
                                                               order='F')))
            elif self.shape == 'Sinc':
                self.bounds = ((-np.inf, 0, 0, 0, 0),
                               (np.inf, np.inf, np.inf, np.inf, np.inf))
        elif type(bounds) == list:
            for i in range(len(bounds)):
                if bounds[i] == None:
                    'print yeah'
                    bounds[i] = -np.inf, np.inf

            x = np.transpose(np.array(bounds))
            self.bounds = tuple(map(tuple, x))

        else:
            print 'Bounds are: ', bounds
            print 'Bounds entered must be in list format, e.g for gaussian: '
            print 'Parameters: x0, fwhm, A, c'
            print 'Input: [(x0,x0),(fw,fw), (A,A), None'
            quit()

    def fit_sinc(self):
        self.shape = 'Sinc'
        self.boundscheck()
        self.popt, self.pcov = curve_fit(sinc, self.xdata, self.ydata,
                                         [self.center, 1, self.maxht,
                                          self.minht], sigma=self.yerrs,
                                         bounds=self.bounds)

    def fit_gaussian(self):
        self.shape = 'Gaussian'
        self.boundscheck()
        # print self.bounds
        # fits a gaussian function to the data and updates the peak shape
        self.popt, self.pcov = curve_fit(gaussian, self.xdata, self.ydata,
                                         [self.center, self.gfwhm,
                                          self.maxht, self.minht],
                                         sigma=self.yerrs, bounds=self.bounds)

    def fit_lorentzian(self):  # MIGHT NOT BE NORMALISED
        self.shape = 'Lorentzian'
        self.boundscheck()
        self.popt, self.pcov = curve_fit(lorentzian, self.xdata, self.ydata,
                                         [self.center, self.gfwhm,
                                          self.maxht, self.minht], sigma=self.yerrs,
                                         bounds=self.bounds)

    def fit_pearsonVII(self):  # Broken
        self.shape = 'PearsonVII'
        self.boundscheck()
        self.popt, self.pcov = curve_fit(pearsonVII, self.xdata, self.ydata,
                                         [self.center, self.gfwhm, 1,
                                          self.maxht, self.minht], sigma=self.yerrs)

    def fit_pseudovoight(self):
        self.shape = 'Pseudovoight'
        self.boundscheck()
        self.popt, self.pcov = curve_fit(pseudovoight, self.xdata, self.ydata,
                                         [self.center, self.gfwhm, 0.5,
                                          self.maxht, self.minht], sigma=self.yerrs,
                                         bounds=self.bounds)

    def fit_asymgauss(self):
        self.shape = 'Asymgauss'
        self.boundscheck()
        self.popt, self.pcov = curve_fit(asymgaussian, self.xdata, self.ydata,
                                         [self.center, self.gfwhm, self.gfwhm,
                                          self.maxht, self.minht], sigma=self.yerrs,
                                         bounds=self.bounds)

    def fit_asymloren(self):
        self.shape = 'Asymloren'
        self.boundscheck()
        self.popt, self.pcov = curve_fit(asymlorentzian, self.xdata, self.ydata,
                                         [self.center, self.gfwhm, self.gfwhm,
                                          self.maxht, self.minht], sigma=self.yerrs,
                                         bounds=self.bounds)

    def multipeak(self, xdata, *args):
        out = np.zeros(np.shape(xdata))  # creating output arra
        args = self.__argsareeasy(*args)
        for i in range(len(self.functions)):
            for j in range(len(self.shapes)):
                if self.shapes[j] == self.functions[i]:
                    function = self.funcs[self.shapes.index(self.functions[i])]
            if i > 0:
                args[i][-1] = 0
            out += function(xdata, *args[i])
        return out + args[-1]

    def fit_multipeak(self, functions, args, plot=False):
        self.shape = 'Multiple'
        self.functions = functions
        self.args = args
        self.boundscheck()
        if plot:
            plt.errorbar(self.xdata, self.ydata, yerr=(self.ydata)**0.5,
                         ls='--', marker='x')
            plt.plot(self.xdata, self.multipeak(self.xdata, args))
            plt.title('ititial parameters')
            plt.show()
        self.popt, self.pcov = curve_fit(self.multipeak,
                                         self.xdata, self.ydata,
                                         p0=args, sigma=self.yerrs, maxfev=5000,
                                         bounds=self.bounds)

    def funcfinder(self):
            # allows the correct peak shape to be used in other functions
        for i in self.shapes:
            if self.shape == i:
                return self.funcs[self.shapes.index(i)]

    def analyse_fit(self):

        f = self.funcfinder()
        chi = np.sum(
            ((self.ydata - f(self.xdata, *self.popt)) / (self.yerrs))**2)
        self.fit = chi / (len(self.xdata)-len(self.popt))
        self.errs = np.sqrt(np.diag(self.pcov))

    def plot_init(self):
        f = self.funcfinder()
        plt.figure()
        plt.errorbar(self.xdata, self.ydata, yerr=(self.ydata)**0.5,
                     ls='none', marker='x')
        newx = np.linspace(self.xdata.min(), self.xdata.max(), 100)
        plt.plot(newx, f(newx, self.center, self.gfwhm, 0.5,
                         self.maxht, self.minht), ls='-')
        plt.show()

    def plot(self, detail=False, hold=False):
        f = self.funcfinder()
        plt.figure()
        plt.errorbar(self.xdata, self.ydata, yerr=(self.ydata)**0.5,
                     ls='none', marker='x')
        newx = np.linspace(self.xdata.min(), self.xdata.max(), 100)
        plt.plot(newx, f(newx, *self.popt), ls='-')
        if self.shape == 'Multiple' and detail:
            x = 0
            params = self.__argsareeasy(self.popt)

            for i in self.functions:
                f = self.funcs[self.shapes.index(i)]
                nprms = np.array(list(params[x])+[params[-1]])
                plt.plot(newx, f(newx, *nprms), ls='--', label=i)
                plt.legend()
                x = x + 1
        if hold != True:
            plt.grid()
            plt.show()

    def area(self, bg=0):

        if b_any(self.shape in x for x in self.shapes[: 4]):
            self.area_anal = self.popt[-2]
        elif self.shape == 'Multiple':
            self.area_anal = 0.
            for i in self.__argsareeasy(self.popt)[:-1]:
                self.area_anal += i[-1]

        if bg == 0:
            bg = (self.xdata.max() - self.xdata.min())*self.popt[-1]
        elif bg == 1:
            bg = (self.ydata[0] + self.ydata[-1]) * \
                (self.xdata.max() - self.xdata.min())
        f = self.funcfinder()
        if (b_any(self.shape in x for x in self.shapes[: 4])
                or self.shape == 'Multiple'):
            total = quad(f, self.xdata.min(), self.xdata.max(),
                         args=tuple(self.popt))
            self.area_quad = total[0] - bg
            self.area_quade = total[1]
        xvals = np.linspace(self.xdata.min(), self.xdata.max(), 10000)
        self.area_trapz = trapz(f(xvals, *self.popt), xvals) - bg

        self.area_raw = trapz(self.ydata, x=self.xdata) - bg


class Peak2d:

    def __init__(self, xdata, ydata, zdata, zerrs=0):

        if type(zerrs) == int:
            self.zerrs = np.sqrt(zdata).ravel()
        else:
            self.zerrs = zerrs.ravel()
        self.xdata = xdata
        self.ydata = ydata
        self.zdata = zdata
        self.gy0, self.gx0 = np.unravel_index(self.zdata.argmax(),
                                              self.zdata.shape)
        self.gbg = zdata[0][0]
        self.gmh = self.zdata.ravel().max()
        self.gfwx = self.xdata[0][8] - self.xdata[0][0]
        self.gfwy = self.xdata[0][8] - self.xdata[0][0]

        self.shapes = ['Gaussian', 'Lorentzian', 'Pseudovoight',
                       'Adv_Gaussian', 'Adv_Lorentzian', 'Adv_Pseudovoight',
                       'Multiple']
        self.funcs = [gaussian2D, lorentzian2D, pseudovoight2D,
                      adv_gaussian2D, adv_lorentzian2D, adv_pseudovoight2D,
                      self.adv_2dmultipeak]
        self.bounds = ((0, 0, 0, 0, 0, 0, 0, 0),  # make a boundermaker func
                       (np.inf,  np.inf,  np.inf,  np.inf,
                        np.inf, np.inf,  np.inf, np.inf))

    def __argsareeasy(self, args):
        args = list(args)
        if len(args) == 1:
            args = args[0]
        fakelist = []
        nargs = []
        x = 0
        for i in self.functions:
            x += 1
            if i == 'Gaussian' or i == 'Lorentzian':
                fakelist.append(5)
            elif (i == 'Pseudovoight' or i == 'Adv_Gaussian'
                  or i == 'Adv_Lorentzian'):
                fakelist.append(6)
            elif i == 'Adv_Pseudovoight':
                fakelist.append(7)
            else:
                print i, 'does not feature in',  self.shapes[:-1]

        x = 0
        for i in fakelist:
            nargs.append(args[x:x+i]+[0])

            x = x + i

        return nargs + [args[-1]]

    def update_init(self, params):

        # changes initial parameters based on previous fits
        self.gx0 = params[0]
        self.gy0 = params[1]
        self.gfwx = params[2]
        self.gfwy = params[3]
        self.gmh = params[-2]
        self.gbg = params[-1]

    def fit_2dgaussian(self):

        self.popt, self.pcov = curve_fit(gaussian2D, (self.xdata, self.ydata),
                                         self.zdata.ravel(),
                                         p0=[self.gx0, self.gy0, self.gfwx,
                                             self.gfwy, self.gmh, self.gbg],
                                         sigma=self.zerrs)

        self.shape = 'Gaussian'

    def fit_adv_2dgaussian(self):

        self.popt, self.pcov = curve_fit(adv_gaussian2D, (self.xdata, self.ydata),
                                         self.zdata.ravel(),
                                         p0=[self.gx0, self.gy0, self.gfwx,
                                             self.gfwy, 0,  self.gmh, self.gbg],
                                         sigma=self.zerrs)

        self.shape = 'Adv_Gaussian'

    def fit_2dlorentzian(self):

        self.popt, self.pcov = curve_fit(lorentzian2D, (self.xdata, self.ydata),
                                         self.zdata.ravel(),
                                         p0=[self.gx0, self.gy0, self.gfwx,
                                             self.gfwy, 10000*self.gmh, self.gbg],
                                         sigma=self.zerrs)

        self.shape = 'Lorentzian'

    def fit_adv_2dgaussian(self):

        self.popt, self.pcov = curve_fit(adv_gaussian2D, (self.xdata, self.ydata),
                                         self.zdata.ravel(),
                                         p0=[self.gx0, self.gy0, self.gfwx,
                                             self.gfwy, 0,  self.gmh, self.gbg],
                                         sigma=self.zerrs)

        self.shape = 'Adv_Gaussian'

    def fit_2dpseudovoight(self):

        self.popt, self.pcov = curve_fit(pseudovoight2D, (self.xdata, self.ydata),
                                         self.zdata.ravel(),
                                         p0=[self.gx0, self.gy0, self.gfwx,
                                             self.gfwy, 0.5,  self.gmh, self.gbg],
                                         sigma=self.zerrs)

        self.shape = 'Pseudovoight'

    def fit_adv_2dpseudovoight(self):

        self.popt, self.pcov = curve_fit(adv_pseudovoight2D, (self.xdata, self.ydata),
                                         self.zdata.ravel(),
                                         p0=[self.gx0, self.gy0, self.gfwx,
                                             self.gfwy, 0.5, 0.,   self.gmh, self.gbg],
                                         sigma=self.zerrs, bounds=self.bounds)

        self.shape = 'Adv_Pseudovoight'

    def adv_2dmultipeak(self, (xdata, ydata), *args):
        out = np.zeros(np.shape(xdata.ravel()))  # creating output arra
        args = self.__argsareeasy(args)
        for i in range(len(self.functions)):
            for j in range(len(self.shapes)):
                if self.shapes[j] == self.functions[i]:
                    function = self.funcs[self.shapes.index(self.functions[i])]
            # fig = plt.figure()
            # ax = fig.add_subplot(111, projection='3d')
            # ax.plot_surface(xdata, ydata,
            #                function((xdata, ydata),
            #                         *args[i]).reshape(np.shape(xdata)))
            # plt.show()
            out += function((xdata, ydata), *args[i])
        return out + args[-1]

    def fit_2dmultipeak(self, functions, args, plot=False):
        self.functions = functions
        if plot:

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_surface(self.xdata, self.ydata, self.zdata,
                            alpha=0.4,  cmap=cm.coolwarm)

            ax.plot_surface(self.xdata, self.ydata,
                            self.adv_2dmultipeak((self.xdata, self.ydata),
                                                 *args).reshape(len(self.zdata),
                                                                len(self.zdata[0])),
                            cmap=cm.winter)

            plt.title('ititial parameters')
            plt.show()
        self.popt, self.pcov = curve_fit(self.adv_2dmultipeak,
                                         (self.xdata, self.ydata),
                                         self.zdata.ravel(), bounds=self.bounds,

                                         p0=args, sigma=self.zerrs,
                                         maxfev=3000)
        self.shape = 'Multiple'

    def funcfinder(self):

        for i in self.shapes:
            if self.shape == i:
                return self.funcs[self.shapes.index(i)]

    def analyse_fit(self):
        f = self.funcfinder()
        chi = np.sum((
            ((self.zdata).ravel() - f((self.xdata, self.ydata), * self.popt))
            / (self.zerrs))**2)
        self.fit = chi / (len(self.xdata)-len(self.popt))
        self.errs = np.sqrt(np.diag(self.pcov))

    def numarea(self):  # Needs fixing for 2D too.

        f = self.funcfinder()
        z = f((self.xdata, self.ydata), *self.popt)
        if len(self.xdata) == len(self.ydata):
            print simps(simps(z.reshape(len(self.ydata[0]), len(self.xdata[0])), self.ydata[0]), self.xdata[0])
            print self.popt[-2]
        else:
            print 'x and y dimensions must be identical - use .crop on image'
        print 'to be implemented'

    def area(self, bg=0):

        if b_any(self.shape in x for x in self.shapes[:6]):
            self.area_anal = pk.popt[-2]
        elif self.shape == 'Multiple':
            self.area_anal = 0.
            for i in self.__argsareeasy(self.popt)[:-1]:
                self.area_anal += i[-1]
        else:
            self.area_anal = 'No Analytical Area'

        if bg == 0:
            bg = (self.xdata.max() - self.xdata.min())*self.popt[-1]
        elif bg == 1:
            bg = (((self.xdata.max() - self.xdata.min()) *
                   (self.ydata.max() - self.ydata.min())) *
                  (self.zdata[0][0] + self.zdata[0][-1] +
                   self.zdata[-1][0] + self.zdata[-1][-1])/4)

    def plot(self):  # and the composing things

        f = self.funcfinder()
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(self.xdata, self.ydata, self.zdata,
                        alpha=0.4,  cmap=cm.coolwarm)
        ax.plot_surface(self.xdata, self.ydata,
                        f((self.xdata, self.ydata),
                          *self.popt).reshape(len(self.zdata),
                                              len(self.zdata[0])),
                        cmap=cm.winter)

        plt.xlabel('xaxis')
        plt.ylabel('yaxis')
        plt.show()

    def saveplot(self, fname):

        f = self.funcfinder()
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(self.xdata, self.ydata, self.zdata,
                        alpha=0.4,  cmap=cm.coolwarm)
        ax.plot_surface(self.xdata, self.ydata,
                        f((self.xdata, self.ydata),
                          *self.popt).reshape(len(self.zdata),
                                              len(self.zdata[0])),
                        cmap=cm.winter)

        plt.xlabel('xaxis')
        plt.ylabel('yaxis')
        plt.savefig(fname+'.png')


class lattice():

    def __init__(self, im):
        self.im = im
        self.X, self.Y = np.meshgrid(np.arange(len(im)),
                                     np.arange(len(im[0])))
        self.boxes = []

    def __argsareeasy(self, args):
        p = args[:2+2*len(self.functions)]
        args = args[2+2*len(self.functions):]
        nargs = []
        self.flatfunctions = self.functions * self.order * 2
        l = 0
        m = 0
        n = 1
        for i in range(self.order):
            i = i + 1
            for j in range(len(self.functions)):
                for k in range(2):
                    nargs.append(p[-2]+n*i*p[0+2*j]*np.cos(
                        np.pi*p[1+2*j]/180.))
                    nargs.append(p[-1]-n*i*p[0+2*j]*np.sin(
                        np.pi*p[1+2*j]/180.))
                    nargs.append(args[l+0])
                    nargs.append(args[l+1])
                    if self.functions[j] == 'Adv_Pseudovoight':
                        nargs.append(args[l+2])
                        l += 1
                    nargs.append(args[l+2]+(n+1)*90)
                    nargs.append(args[l+3])
                    n = - n
                if self.functions[j] == 'Adv_Pseudovoight':
                    l += 5
                else:
                    l += 4

        nargs.append(args[-1])
        # print nargs
        return tuple(nargs)

    def plotim(self):
        fig = plt.figure()
        plt.imshow(self.im)
        plt.show()

    def plotlogim(self):
        fig = plt.figure()
        plt.imshow(self.im)
        plt.show()

    def trybox(self, pixels, location):
        if (location[1]+pixels[1] > len(self.im) or
                location[0]+pixels[0] > len(self.im[0])):
            print 'box out of bounds'

        fig, ax = plt.subplots(1)
        rect = patches.Rectangle((location[0], location[1]),
                                 pixels[0], pixels[1],
                                 linewidth=1, edgecolor='r', facecolor='none')
        ax.imshow(self.im, aspect='auto')
        ax.add_patch(rect)
        plt.show()

    def makebox(self, pixels, location):
        self.boxes.append([pixels[0], pixels[1], location[0], location[1]])

    def showboxes(self):
        fig, ax = plt.subplots(1)
        ax.imshow(self.im, aspect='auto')
        for i in self.boxes:
            rect = patches.Rectangle((i[2], i[3]),
                                     i[0], i[1],
                                     linewidth=1, edgecolor='r',
                                     facecolor='none')
            ax.add_patch(rect)

        plt.show()

    def makelattice(self):
        for i in range(len(self.boxes)):
            z = self.im[self.boxes[i][3]: self.boxes[i][3]+self.boxes[i][1],
                        self.boxes[i][2]: self.boxes[i][2]+self.boxes[i][0]]
            x, y = np.meshgrid(
                range(self.boxes[i][2], self.boxes[i][2]+len(z[0])),
                range(self.boxes[i][3], self.boxes[i][3]+len(z)))

            if i == 0:
                self.xs = x
                self.ys = y
                self.zs = z
            else:
                self.xs = np.append(self.xs, x, axis=1)
                self.ys = np.append(self.ys, y, axis=1)
                self.zs = np.append(self.zs, z, axis=1)
        self.zerrs = np.sqrt(self.zs)

    def autolattice(self, dim, geometry):
        self.boxes = []
        for i in range(self.order):
            i = i + 1
            for j in range(geometry):
                self.popt[1] = self.popt[1] + j*60.
                self.makebox(
                    (dim[0], dim[1]),
                    (int(self.popt[2] +
                         i*self.popt[0]*np.cos(np.pi*self.popt[1]/180))-dim[0]/2.,
                     int(self.popt[2] -
                         i*self.popt[0]*np.sin(np.pi*self.popt[1]/180))-dim[1]/2.))

                self.makebox(
                    (dim[0], dim[1]),
                    (int(self.popt[2] -
                         i*self.popt[0]*np.cos(np.pi*self.popt[1]/180))-dim[0]/2.,
                     int(self.popt[2] +
                         i*self.popt[0]*np.sin(np.pi*self.popt[1]/180))-dim[1]/2.))

    def plot_raw_a(self):

        fig = plt.figure()
        ax = fig.add_subplot(121, projection='3d')
        ax.plot_surface(self.xs[:, :50], self.ys[:, :50], self.zs[:, :50]-1000,
                        alpha=0.4,  cmap=cm.coolwarm)

        plt.xlabel('x')
        plt.ylabel('y')

        ax = fig.add_subplot(122, projection='3d')
        ax.plot_surface(self.xs[:, 50:], self.ys[:, 50:], self.zs[:, 50:]-1000,
                        alpha=0.4,  cmap=cm.coolwarm)

        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def pair_peaks(self, (xdata, ydata), *args):  # needs rethink
        nargs = self.__argsareeasy(args)
        pk = Peak2d(self.xs, self.ys, self.zs)
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_surface(self.xs, self.ys, self.zs)
        # plt.show()
        pk.functions = self.flatfunctions
        # print nargs
        # print pk.functions
        # quit()
        return pk.adv_2dmultipeak((pk.xdata, pk.ydata),
                                  nargs)

    def fit_pairedpeaks(self, args, plot=False):
        if plot:
            fit = self.pair_peaks((self.xs, self.ys), *args)

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_surface(self.xs, self.ys, self.zs,
                            alpha=0.4,  cmap=cm.coolwarm)

            ax.plot_surface(self.xs, self.ys,
                            fit.reshape(len(self.zs), len(self.zs[0])),
                            cmap=cm.winter)

            plt.title('initial parameters')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.show()
        if hasattr(self, 'bounds'):
            self.popt, self.pcov = curve_fit(self.pair_peaks,
                                             ((self.xs).ravel(),
                                              (self.ys).ravel()),
                                             (self.zs).ravel(),
                                             p0=args, sigma=(
                                                 self.zerrs).ravel(),
                                             bounds=self.bounds,
                                             maxfev=50000)
            self.shape = 'Pairedpeaks'
        else:
            self.popt, self.pcov = curve_fit(self.pair_peaks,
                                             ((self.xs).ravel(),
                                              (self.ys).ravel()),
                                             (self.zs).ravel(),
                                             p0=args, sigma=(
                                                 self.zerrs).ravel(),
                                             maxfev=50000)
            self.shape = 'PairedPeaks'


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    for j in range(5000, 40001, 5000):
        data = np.transpose(np.loadtxt(
            '/home/sam/SamL/SpinTextures/'+str(j)+'_datafromatoms.txt',
            skiprows=140))
        print trapz(data[1], x=data[0])
    # pk.area()
    # print pk.area_anal, pk.area_quad, pk.area_trapz, pk.area_raw
    # print
    # pk.fit_lorentzian()
    # pk.area()
    # print pk.area_anal, pk.area_quad, pk.area_trapz, pk.area_raw
    # print
    # pk.fit_pseudovoight()
    # pk.area()
    # print pk.area_anal, pk.area_quad, pk.area_trapz, pk.area_raw
    # print
    # pk.area()
    # pk.analyse_fit()
    # print pk.popt
    # print pk.errs
    # print
    # print pk.area_anal, pk.area_quad, pk.area_trapz, pk.area_raw
    # pk.fit_asymloren() ; pk.analyse_fit() ; print pk.fit
    # pk.fit_gaussian() ; pk.analyse_fit()  ; print pk.fit
    # pk.fit_pseudovoight(); pk.analyse_fit();print pk.fit
    # pk.fit_lorentzian() ; pk.analyse_fit(); print pk.fit
    # pk.plot()
'''
    def h5toarray():
        im = h5.File('test_images/skyrmion.nxs', 'r')
        im = scipy.squeeze(im[im.keys()[0]]["scan_data"]["data_14"])
        return im

    im = lattice(np.log(h5toarray()))
    im.makebox((100, 100), (550, 1730))
    im.makebox((100, 100), (1050,  580))
    im.makelattice()
    im.fit_pairedpeaks(['Adv_Pseudovoight', 'Adv_Pseudovoight'],
                       args=(630, 66, 847, 1210, 20, 20, 0, 0.5, 50, 10.7),
                       plot=False)

    zdata = im.paired_peaks((im.xs, im.ys), *tuple(map(tuple, [im.popt]))[0])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(im.xs, im.ys, zdata.reshape(
        np.shape(im.xs)), cstride=5,  rstride=5, cmap=cm.coolwarm)
    ax.plot_surface(im.xs+100, im.ys, im.zs, cstride=5,  rstride=5,
                    cmap=cm.viridis)
    print im.popt
    plt.show()

    data = np.transpose(np.loadtxt(
        '/home/sam/Data/i16_Diamond/GdPd2Si3/730435.dat',
        skiprows=140))

    from PIL import Image

    im = Image.open(
        '/home/sam/Data/i16_Diamond/GdPd2Si3/731599-pilatus3_100k-files/00008.tif').crop([190, 50, 300, 160])
    zdata = np.array(im)
    # plt.imshow(zdata)
    # plt.show()
    x, y = np.meshgrid(range(len(zdata[0])), range(len(zdata)))

    pk = Peak2d(x, y, zdata)
    pk.fit_adv_2dgaussian()
    pk.plot()


    def h5toarray():
        im = h5.File('test_images/skyrmion.nxs', 'r')
        im = scipy.squeeze(im[im.keys()[0]]["scan_data"]["data_14"])
        return im

    im = lattice(np.log(h5toarray()))
    im.makebox((100, 100), (550, 1730))
    im.makebox((100, 100), (1050,  580))
    im.fit_pairedpeaks(['Adv_Pseudovoight', 'Adv_Pseudovoight'], [[
                       630, 66, [847, 1210], [20, 20, 0, 0.5, 50, 10.7]]])


    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_surface(x, y, np.log(a), cmap=cm.plasma, cstride=5, rstride=5


    import matplotlib.pyplot as plt

    data = np.transpose(np.loadtxt(
        '/home/sam/Data/i16_Diamond/GdPd2Si3/730435.dat',
        skiprows=140))

       from PIL import Image

       im = Image.open(
           '/home/sam/Data/i16_Diamond/GdPd2Si3/731599-pilatus3_100k-files/00008.tif').crop([190, 50, 300, 160])
       zdata = np.array(im)
       # plt.imshow(zdata)
       # plt.show()
       x, y = np.meshgrid(range(len(zdata[0])), range(len(zdata)))

       pk = Peak2d(x, y, zdata)
       pk.fit_adv_2dgaussian()

       # pk.fit_2dmultipeak(['Adv_Pseudovoight', 'Adv_Pseudovoight'],
       #                   [54, 46, 8, 8, 0, 0.5, 21030,
       #                    55, 43, 7, 7, 0, 0.5, 20000, 1279], plot=False)

       pk.area()

       newparams = pk.popt
       for i in range(1, 18):
           im = Image.open('/home/sam/Data/i16_Diamond/GdPd2Si3/731599-pilatus3_100k-files/' +
                           str(i).zfill(5)+'.tif').crop([150, 50, 300, 150])
           zdata = np.array(im)
           x, y = np.meshgrid(range(len(zdata[0])), range(len(zdata)))
           pk = Peak2d(x, y, zdata)
           pk.update_init(newparams)
           pk.fit_2dgaussian()
           pk.plot()
'''
