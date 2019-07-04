from Peak import *
from pyvtk import *
import scipy
import scipy.misc


class double_peak:


if __name__ == "__main__":
    from mpl_toolkits.mplot3d import Axes3D
    from PIL import Image
    from matplotlib import cm

    data = np.array(Image.open('test_images/cone.jpg'))
    x, y = np.meshgrid(np.arange(len(data)), np.arange(len(data[0])))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, data)
    plt.show()
