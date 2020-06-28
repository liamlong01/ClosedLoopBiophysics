import numpy as np
import matplotlib.pyplot as plt

import util

import pdb


class CLStim():
    sigma = 0.3
    current = 100e-6

    def __init__(self,  verbose = True, label="Unlabeled CLStim"):
        """
        network is a collection of NetworkCells
        """
        self.trigger = 0
      

        self.label = label
        self.clims = [-100e-6, 100e-6]


        self.stimming = 0
        self.currents = []
        self.current = 0

        self.verbose = verbose


    def doStim(self, time, electrode, LFPs, hocObj):
        
        if self.verbose:
            pdb.set_trace()
            print("CLSTIM received LFP with shape:", LFPs.shape)


        LFPs = self.addNoise(LFPs)
        self.trigger = 1
        if self.time > 20e-3:
            self.biphasicpulse()

        if np.max(LFPs) > 0.01 and not self.stimming:
            self.plot(LFPs)

        self.currents = np.append(self.currents, self.current)
        print(self.time, "s, current setting:",
              self.current, "maxlfp", np.max(LFPs))

        microStim(self.current, 0, 0, 0, 0.3, hocObj)

    def biphasicpulse(self):
        if not self.stimming:
            if self.trigger:
                self.stimming = True
                self.start = self.time

        else:
            if self.time - self.start < 100e-6:
                self.current = 100e-6

            elif self.time - self.start < 200e-6:
                self.current = -100e-6

            elif self.time - self.start < 1e-3:
                self.current = 0

            else:
                self.stimming = 0
                self.start = 0
                self.trigger = 0

    def plot(self, lfps):
        vec = lfps.flatten()
        if np.any(vec > self.clims[1]):
            self.clims[1] = np.max(vec)
        if np.any(vec < self.clims[0]):
            self.clims[0] = np.min(vec)

        plt.clf()
        plt.subplot(2, 1, 1)
        im = plt.imshow(lfps, aspect='auto')
        cbar = plt.colorbar()
        im.set_clim(self.clims[0], self.clims[1])
        plt.title("LFPS")
        # self.csd(lfps)

        plt.subplot(2, 1, 2)
        plt.plot(self.currents)
        plt.plot()

        plt.pause(0.01)

    def csd(self, lfps):
        csd = cv.Laplacian(lfps, cv.CV_64F, ksize=3)
        plt.subplot(1, 2, 2)
        plt.clf()
        im = plt.imshow(csd, aspect='auto')
        cbar = plt.colorbar()
        plt.title("CSD")

        plt.pause(0.05)

    def addNoise(self, matrix):

        rms = self.calcNoiseRMS()
        matrix = np.add(matrix, np.random.normal(scale=rms, size=matrix.shape))
        return matrix

    def calcNoiseRMS(self, l=6e-6, w=6e-6, rlead=30e3):
        """
        Noise model from:

        Can One Concurrently Record Electrical Spikes from Every Neuron in a Mammalian Brain?
        David Kleinfeld, Lan Luan, Partha P.Mitra, Jacob T.Robinson, Rahul Sarpeshkar, Kenneth Shepard, Chong Xie, Timothy D.Harris
        Sept. 25, 2019

        This is just a rough estimate
        All info is taken from this paper and credit to the authors for compiling these numbers and equations
        """
        eps0 = 8.8541878182e-12
        eps_double = 2
        kB = 1.38064852e-23
        T = 310.05  # Kelvin - 36.9 celsius same as internal body

        Apad = l * w

        # tdouble from:
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

    def print_info(self, info):
        print(self.label + '----' + info)


def microStim(current, x, y, z, sigma, hocObj):
    print("calling microstim....", current)
    count = 0
    hocObj._ref_stim_xtra[0] = 1

    for sec in hocObj.allsec():

        interval = 100
        if hocObj.ismembrane("xtra", sec=sec):
            for seg in sec:

                r = np.sqrt((seg.x_xtra - x)**2 + (seg.y_xtra - y)
                            ** 2 + (seg.z_xtra - z)**2)
                if r == 0:
                    r = 1e-9
                """
                if count%interval == 0:
                    print(seg, "({},{},{}) ".format(x,y,z), "({},{},{}) ".format(seg.x_xtra,seg.y_xtra,seg.z_xtra), "voltage is:", current*1e-3/(4*np.pi*sigma*r))
                    count = 0
                """
                seg._ref_es_xtra[0] = current * 1e-3 / (4 * np.pi * sigma * r)

                count += 1
                print("{} e_extracellular: {}, global stim: {}, es_xtra: {}".format(
                    seg, seg.e_extracellular, hocObj("stim_xtra"), seg.es_xtra))
