


import LFPyStim as LFPy
import numpy as np
import matplotlib
matplotlib.use('TkAgg')


import matplotlib.pyplot as plt
import pickle

import sys
import pdb

import util

from matplotlib import animation

def loadcell(cellstr, folder=''):

    f = open(folder + cellstr, 'rb')

    cell = pickle.load(f)

    f.close()

    return cell


def electrode_LFP(cell):

    res = 10

    xmin = min(cell.xend)

    ymin = int(min(cell.yend) / res) - 5
    ymax = int(max(cell.yend) / res) + 5

    zmin = int(min(cell.zend) / res) - 5
    zmax = int(max(cell.zend) / res) + 5

    X, Y, Z = np.mgrid[1:2, ymin:ymax:1, zmin:zmax:1] * res
    X = X - xmin
    # define parameters for extracellular recording electrode,
    # using optional method

    electrodeParameters = {
        'sigma': 0.3,  # extracellular conductivity
        'x': X.flatten(),  # x,y,z-coordinates of contacts
        'y': Y.flatten(),
        'z': Z.flatten(),
        'N': np.array([[1, 0, 0]] * X.size),  # surface normals
        'r': 6,  # contact site radius
        'n': 20,  # datapoints for averaging
    }

    print("creating electrode...")

    print("z-direction from : %d %d" % (zmin, zmax))
    print("y-direction from : %d %d" % (ymin, ymax))

    print("total electrodes: %d" % len(X.flatten()))
    # create extracellular electrode object for LFPs on grid
    electrode = LFPy.RecExtElectrode(cell, **electrodeParameters)

    # calcing electrode LFP
    print("calcing  electrode LFP...")
    electrode.calc_lfp()
    return electrode, X, Y, Z


def getspike_basic(lfp):
    pk = np.unravel_index(np.argmax(lfp), lfp.shape)

    spike = lfp[:, pk[1] - 70:pk[1] + 50]
    return spike


def debug_plots(cell, electrode, cellsoma=None, electrodeLFP=None):

    if cellsoma:
        plt.figure()
        plt.plot(cell.somav)

    if electrodeLFP:
        print("_*_*_*_*_ ELEC LFP HERE _*_*_*_*_*_*_*__")
        """
        plt.figure()
        plt.matshow(ele`ctrode.LFP)
        plt.colorbar()
        plt.axis('tight')
        """
    plt.show()


def update(frame, spike, shape, im, im2, ln):
    im.set_array(spike[:, frame].reshape(shape))
    # im2.set_array(util.csd(spike[:, frame].reshape(shape)))
    im2.set_array(util.csd(spike[:, frame].reshape(shape)))
    ln.set_data(np.linspace(0, frame - 1, frame), spike[200, :frame])
    return [im, im2, ln]


def animate_spike(spike, shape):
    print("animate spike")
    fig1 = plt.figure()
    plt.subplot(1,2,1)
    print("animate spike")
    ax1 = plt.gca()
    print("imshow!!!")
    im = ax1.imshow(spike[:, 0].reshape(shape))
    im.set_clim(np.min(spike.flatten()), np.max(spike.flatten()))
    fig1.colorbar(im, ax=ax1)
   
    plt.subplot(1,2,2)
    ax3 = plt.gca()
    print("imshow!!!")
    im3 = ax3.imshow(spike[:,0].reshape(shape))
    cmax = np.max([np.max(util.csd(x).flatten()) for x in spike])
    cmin = np.min([np.min(util.csd(x).flatten()) for x in spike])

    im3.set_clim(cmax, cmin)
    fig1.colorbar(im3, ax=ax3)

    plt.figure()
    ax2 = plt.axes(xlim=(0, 500), ylim=(np.min(spike.flatten()),
                                 np.max(spike.flatten())))

    ln, = ax2.plot([], [])

    update_func = lambda frame, spike=spike, shape=shape, im=im, im3=im3, ln=ln: \
                                               update(frame, spike, shape, im, im3, ln)
    
    print("trying animation...")
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Liam Long'), bitrate=1800)

    anim = animation.FuncAnimation(fig1, update_func, init_func = lambda : update_func(0),
                               frames=120, interval=5, blit=True)

    anim.save('basic_animation.mp4', writer=writer)

    plt.show()


def main():
    celltoload = "L4_ChC_cACint209_2.results"
    directory = "./results/cellvoltages_axons_"

    cell = loadcell(celltoload, directory)

    electrode, X, Y, Z  = electrode_LFP(cell)
    print("_________ELECTRODE INFO______________________\nY shape: ")
    print(Y.shape)
    print(electrode.LFP.shape)
    print("____________________________________________")

    shape = Y.shape[1:]

    spike = getspike_basic(electrode.LFP)
    noisy_spike = util.addNoise(spike)
    # util.wrap_kcsd(X, Y, Z, spike)
    animate_spike(noisy_spike, shape)

    return cell, electrode, spike


if __name__ == '__main__':

    electrodeLFP = '-elecplot' in sys.argv
    cellsoma = '-cellsoma' in sys.argv

    print(cellsoma)
    
    cell, electrode, spike = main()
    debug_plots(cell, electrode, cellsoma=cellsoma, electrodeLFP=electrodeLFP)

    
