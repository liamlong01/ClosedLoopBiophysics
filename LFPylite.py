"""
A remake of LFPy library with a limited but more focused feature set

Dependencies:
 lfpcalc.py from the original LFPy package
"""

import neuron  # this currently has just one use (n3d -x3d etc)
# TODO remove neuron import at some point with slight refactoring

import sys
import pdb
import numpy as np
from aberraAxon import MyelinatedCell

import lfpcalc
import matplotlib.pyplot as plt
import pickle


class CellData():

    """
    Contains all relevant cell data in a picklable data format
    """

    def __init__(self, cell):
        """
        This is adapted directly from the
        run_simulation._collect_geometry_neuron
        function found in the LFPy library.

        Copyright (C) 2012 Computational Neuroscience Group, NMBU.
        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILdITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more detail
        """
        self.imem = None
        if hasattr(cell, 'totnsegs'):
            print("totnsegs")
            raise NotImplementedError

        else:
            self.totseg = 0
            for sec in cell.sections:
                self.totseg += sec.nseg

        self.area = np.zeros(self.totseg)
        self.diam = np.zeros(self.totseg)
        self.length = np.zeros(self.totseg)

        self.xstart = np.zeros(self.totseg)
        self.xend = np.zeros(self.totseg)
        self.ystart = np.zeros(self.totseg)
        self.yend = np.zeros(self.totseg)
        self.zstart = np.zeros(self.totseg)
        self.zend = np.zeros(self.totseg)

        self.x = np.zeros(self.totseg)
        self.y = np.zeros(self.totseg)
        self.z = np.zeros(self.totseg)

        counter = 0

        for sec in cell.sections:

            n3d = int(neuron.h.n3d(sec=sec))
            nseg = sec.nseg
            gsen2 = 1. / 2 / nseg
            if n3d:
                L = np.zeros(n3d)
                x = np.zeros(n3d)
                y = np.zeros(n3d)
                z = np.zeros(n3d)
                for i in range(n3d):
                    L[i] = neuron.h.arc3d(i)
                    x[i] = neuron.h.x3d(i)
                    y[i] = neuron.h.y3d(i)
                    z[i] = neuron.h.z3d(i)
                L /= sec.L

                segx = np.zeros(nseg)
                for i, seg in enumerate(sec):
                    segx[i] = seg.x

                segx0 = (segx - gsen2).round(decimals=6)
                segx1 = (segx + gsen2).round(decimals=6)

                # temporary store position of segment midpoints
                segx = np.zeros(nseg)
                for i, seg in enumerate(sec):
                    segx[i] = seg.x

                # can't be >0 which may happen due to NEURON->Python float transfer:
                segx0 = (segx - gsen2).round(decimals=6)
                segx1 = (segx + gsen2).round(decimals=6)

                # fill tors with interpolated coordinates of start and end points
                self.xstart[counter:counter + nseg] = np.interp(segx0, L, x)
                self.xend[counter:counter + nseg] = np.interp(segx1, L, x)

                self.ystart[counter:counter + nseg] = np.interp(segx0, L, y)
                self.yend[counter:counter + nseg] = np.interp(segx1, L, y)

                self.zstart[counter:counter + nseg] = np.interp(segx0, L, z)
                self.zend[counter:counter + nseg] = np.interp(segx1, L, z)

                # fill in values area, diam, length
                for i, seg in enumerate(sec):
                    self.area[counter] = neuron.h.area(seg.x)
                    self.diam[counter] = seg.diam
                    self.length[counter] = sec.L / nseg

                    counter += 1


class basicsim():
    # to get lfp mapping we need:
    # xstart, mid, end
    # same for y
    # same for z

    # xzy start/END IN CELL._REAL_POSITIONS - NEEDS REL_START REL_END and somapos

    # x,y,z coordinate
    # sigma
    # r_limit

    def __init__(self, celldata):
        self.celldata = celldata

    def sim(self, mcell):
        cell = mcell.cell
        # stimuli = self.create_stimuli(cell, 3)

        recordings = {}

        recordings['time'] = neuron.h.Vector()
        recordings['soma(0.5)'] = neuron.h.Vector()
        # recordings['isoma(0.5)'] = neuron.h.Vector()

        recordings['time'].record(neuron.h._ref_t, 0.1)
        # recordings['soma(0.5)'].record(init.CellLoader.cell.soma[0](0.5)._ref_v, 0.1)
        recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, 0.1)
        # recordings['isoma(0.5)'].record(cell.soma[0](0.5)._ref_i_membrane, 0.1)

        neuron.h.tstop = 2000
        # create_stimuli(init.CellLoader.cell, 3)

        # assert(sys.getrefcount(stimuli) > 1)
        memireclist = neuron.h.List()
        for sec in mcell.sections:
            for seg in sec:

                memirec = neuron.h.Vector()
                memirec.record(seg._ref_i_membrane, .1)
                memireclist.append(memirec)

        print('Disabling variable timestep integration')
        neuron.h.cvode_active(0)

        neuron.h.finitialize(-65)
        print("finitialize")
        neuron.h.fcurrent()
        print('Running for %f ms' % neuron.h.tstop)

        counter = 0

        neuron.h.t = 0
        interval = 1000
        while neuron.h.t < neuron.h.tstop:

            neuron.h.fadvance()
            if counter == interval:

                print(neuron.h.t, cell.soma[0](
                    0.5)._ref_v[0], cell.soma[0](0.5)._ref_i_membrane[0])
                counter = 0

                lfps = np.matmul(self.mapping.T, np.array(memireclist))
               
                plt.clf()
                plt.subplot(2, 1, 1)
                plt.imshow(lfps[:, -1].reshape([10, 100]))
                plt.subplot(2, 1, 2)
                plt.plot(recordings['time'], recordings['soma(0.5)'])
                plt.pause(0.05)


            counter += 1

        self.celldata.imem = np.array(memireclist)

        time = np.array(recordings['time'])
        soma_voltage = np.array(recordings['soma(0.5)'])

        soma_voltage_filename = 'soma_voltage_step.dat'

        np.savetxt(
            soma_voltage_filename,
            np.transpose(
                np.vstack((
                    time,
                    soma_voltage))))

        print('Soma voltage for step saved to: %s'
              % soma_voltage_filename)

        # import matplotlib.pyplot as plt
        plt.plot(recordings['time'], recordings['soma(0.5)'])
        plt.show()

    def create_stimuli(self, cell, step_number):
        """Create the stimuli"""

        print('Attaching stimulus electrodes')

        stimuli = []
        step_amp = [0] * 3

        with open('current_amps.dat', 'r') as current_amps_file:
            first_line = current_amps_file.read().split('\n')[0].strip()
            hyp_amp, step_amp[0], step_amp[1], step_amp[2] = first_line.split(
                ' ')

        iclamp = neuron.h.IClamp(0.5, sec=cell.soma[0])
        iclamp.delay = 700
        iclamp.dur = 2000
        iclamp.amp = float(step_amp[step_number - 1])
        print('Setting up step current clamp: '
              'amp=%f nA, delay=%f ms, duration=%f ms' %
              (iclamp.amp, iclamp.delay, iclamp.dur))

        stimuli.append(iclamp)

        hyp_iclamp = neuron.h.IClamp(0.5, sec=cell.soma[0])
        hyp_iclamp.delay = 0
        hyp_iclamp.dur = 3000
        hyp_iclamp.amp = float(hyp_amp)
        print('Setting up hypamp current clamp: '
              'amp=%f nA, delay=%f ms, duration=%f ms' %
              (hyp_iclamp.amp, hyp_iclamp.delay, hyp_iclamp.dur))

        stimuli.append(hyp_iclamp)

        return stimuli


if __name__ == '__main__':
    mcell = MyelinatedCell()
    mcell.loadcell(16, myelinate_ax=True)
    celld = CellData(mcell)
    sim = basicsim(celld)
    res = 10
    xmin = min(celld.xend)

    ymin = int(min(celld.yend) / res) - 5
    ymax = int(max(celld.yend) / res) + 5

    zmin = int(min(celld.zend) / res) - 5
    zmax = int(max(celld.zend) / res) + 5

    X, Y, Z = np.mgrid[1:2, -50:50:1, -5:5:1] * res
    X = X - xmin
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()

    sim.mapping = np.zeros([celld.totseg, X.shape[0]])

    for i, (x, y, z) in enumerate(zip(X, Y, Z)):
        # print(i,x,y,z)
        sim.mapping[:, i] = lfpcalc.calc_lfp_linesource(
            celld, x, y, z, 0.3, celld.diam)
        # print(x,y,z, mapping.shape)

    sim.sim(mcell)


    with open('testpickle.celldata', 'wb') as pfile:
        pickle.dump(celld, pfile)
