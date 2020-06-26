import cv2 as cv
import scipy.ndimage as nd

import numpy as np

from kcsd import KCSD2D
import matplotlib.pyplot as plt
import pdb

def csd(lfp):
    return -nd.gaussian_laplace(lfp, 5)

def plotcell():
    pass


def wrap_kcsd(X, Y, Z, lfp, h=50, sigma=0.3):
    """
    reference: https://github.com/Neuroinflab/kCSD-python
    """
    ele_pos = np.vstack((Y.flatten(), Z.flatten())).T
    lfp_tmp = lfp[:,0].reshape(len(X.flatten()), 1)
    pdb.set_trace()
    k = KCSD2D(ele_pos, lfp)
    pdb.set_trace()
    

    return k


def addNoise(matrix):
    rms = calcNoiseRMS()
    matrix = np.add(matrix, np.random.normal(scale=rms, size=matrix.shape))
    return matrix


def calcNoiseRMS(l=6e-6, w=6e-6, rlead=30e3):
        """
        Noise model from:

        Can One Concurrently Record Electrical Spikes from Every Neuron in a Mammalian Brain?
        David Kleinfeld, Lan Luan, Partha P.Mitra, Jacob T.Robinson, Rahul Sarpeshkar, Kenneth Shepard, Chong Xie, Timothy D.Harris
        Sept. 25, 2019

        This is just a rough estimate
        All info int his func is taken from this paper and credit to the authors for compiling these numbers and equations
        """
        eps0 = 8.8541878182e-12 # free space permitivity
        eps_double = 2 # er for double layer
        kB = 1.38064852e-23 # boltzmann constant
        T = 310.05  # Kelvin - 36.9 celsius same as internal body

        Apad = l * w # area of electrode

        # tdouble - thickness of double layerfrom:
        # Y. Burak, D. Andelman Hydration interactions: aqueous solvent effects in
        # electric double layers Phys. Rev. E, 62 (2000), pp. 5295-5312
        t_double = 0.3e-9

        Cpad = eps0 * eps_double * Apad / t_double

        # coating can increase noise threshold by effective increase in surface area e.g. TIn, PEDOT
        # a conservative increase in capacitance (some say 100times cacitance)
        factor = 20
        delt_Vpad = np.sqrt(kB * T / (Cpad * factor))

        fAP = 30e3
        delt_Vlead = np.sqrt(4 * kB * T * rlead * fAP)

        # ignoring spreading resistance and degradation of coupling
        return delt_Vpad + delt_Vlead
