import numpy as np 
from scipy.special import gamma

def gaussian(X, Mu, Sigma): # make it a FWHM matrix
    '''
    Stuff about function

    Note Sigma is the covariance matrix not FWHM
    '''
    k = np.shape(X[0])[0] 
    return (np.exp(-0.5*(np.matmul(np.matmul((X-Mu), np.linalg.inv(Sigma)), np.transpose(X-Mu))))
           / np.sqrt((2*np.pi)**k * np.linalg.det(Sigma))).diagonal() # Need to make more efficient it's taking too long :( - throwing away many diagonal elements here
           )

def lorenzian(X, Mu, Sigma): # make it a fwhm matrix
    '''
    Stuff about function

    Note Sigma is the covariance matrix not FWHM
    '''
    k = np.shape(X[0])[0]  # Need to make more efficient it's taking too long :( - throwing away many diagonal elements here
    return (gamma((k+1.)/2.) / (gamma(0.5)*np.pi**(k/2.)*np.sqrt(np.linalg.det(Sigma))*
                              (1 + np.matmul(np.matmul((X-Mu), np.linalg.inv(Sigma)),
                               np.transpose(X-Mu)))**((1.+k)/2.))).diagonal()

def pseudovoight(X, Mu, Sigma, eta):
    '''
    Stuff about function

    Note Sigma is NOT the covariance matrix and is in fact the FWHM
    '''
    SigmaG = Sigma*2.*np.sqrt(2.*np.log(2.))
    SigmaL = Sigma/2.

    return eta*gaussian(X, Mu, SigmaG) + (1-eta)*lorenzian(X, Mu, SigmaL)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    import time
    xs, ys = np.meshgrid(np.linspace(0, 10, 100), np.linspace(0, 10, 100))
    X = np.reshape(np.append(xs,ys), (10000, 2), order = 'F')

    t0 =  time.time()
    zs = gaussian(X, np.array([5, 5]), np.array([2, 0, 0, 3]).reshape((2,2))).reshape(100,100)
    print time.time() - t0
    ax.plot_surface(xs, ys, zs, cmap = plt.cm.coolwarm)
    plt.show()
