"""
A remake of LFPy library with a limited but more focused feature set

Dependencies:
 lfpcalc.py from the original LFPy package
"""

import neuron  # this currently has just one use (n3d -x3d etc)
# TODO remove neuron import at some point with slight refactoring

from  multiprocessing import Pipe,Process

import sys
import pdb
import numpy as np
from aberraAxon import MyelinatedCell

import lfpcalc
import matplotlib.pyplot as plt
import pickle
import plot3d



class CellData():

    """
    Contains all relevant cell data in a picklable data format
    """

    def shift(self, x, y,z ):
        self.data['xstart'] += x 
        self.data['xend'] += x 
        self.data['ystart'] += y
        self.data['yend'] += y 
        self.data['zstart'] += z 
        self.data['zend'] += z 


    def __iter__():
        self.idx = 0
        self.high = len(self.xstart)
        return self

    def __next__():
        output = self.data[self.idx]
        self.idx+=1
        return output

         

    def __init__(self, cell, hocObj = None):
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
        if hocObj is None:
            hocObj = neuron.h


        
        if hasattr(cell, 'totnsegs'):
            print("totnsegs")
            raise NotImplementedError

        else:
            self.totseg = 0
            for sec in cell.sections:
                self.totseg += sec.nseg

        # using a structured array
        dtype = {
        'names': ('xstart','xend','ystart', 'yend', 'zstart','zend', 'area', 'diam', 'length'),
        'formats': (float, float,  float,   float,   float,   float,  float,  float,  float)
        }

        self.data = np.zeros(self.totseg, dtype=dtype)
        counter = 0

        for sec in cell.sections:

            n3d = int(hocObj.n3d(sec=sec))
            nseg = sec.nseg
            gsen2 = 1. / 2 / nseg
            if n3d:
                L = np.zeros(n3d)
                x = np.zeros(n3d)
                y = np.zeros(n3d)
                z = np.zeros(n3d)
                for i in range(n3d):
                    L[i] = hocObj.arc3d(i)
                    x[i] = hocObj.x3d(i)
                    y[i] = hocObj.y3d(i)
                    z[i] = hocObj.z3d(i)
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
                self.data['xstart'][counter:counter + nseg] = np.interp(segx0, L, x)
                self.data['xend'][counter:counter + nseg] = np.interp(segx1, L, x)

                self.data['ystart'][counter:counter + nseg] = np.interp(segx0, L, y)
                self.data['yend'][counter:counter + nseg] = np.interp(segx1, L, y)

                self.data['zstart'][counter:counter + nseg] = np.interp(segx0, L, z)
                self.data['zend'][counter:counter + nseg] = np.interp(segx1, L, z)

                # fill in values area, diam, length
                for i, seg in enumerate(sec):
                    self.data['area'][counter] = hocObj.area(seg.x)
                    self.data['diam'][counter] = seg.diam
                    self.data['length'][counter] = sec.L / nseg

                    counter += 1

        # we need these references to have exactly thes names to be compatible with LFP libraries
        # but is also convenient to have the structured array format
        self.xstart = self.data['xstart']
        self.xend = self.data['xend']

        self.ystart = self.data['ystart']
        self.yend = self.data['ystart']

        self.zstart = self.data['zstart']
        self.zend = self.data['zstart']

        self.area = self.data['area']
        self.diam = self.data['diam']
        self.length = self.data['length']

        self.imem = np.zeros(self.xstart.shape)



class basicsim():
    # to get lfp mapping we need:
    # xstart, mid, end
    # same for y
    # same for z

    # xzy start/END IN CELL._REAL_POSITIONS - NEEDS REL_START REL_END and somapos

    # x,y,z coordinate
    # sigma
    # r_limit
    doDebug = False

    def __init__(self, celldata, hocObj, renderPipe = None, title='A cell soma plot'):
        self.celldata = celldata
        self. hocObj = hocObj

        self.title = title
        self.mapping = None

        self.renderPipe = renderPipe

    def  debug(self):
        def dumFunc():
            return 0
        if self.doDebug:
            return pdb.set_trace
        else:
            return dumFunc
     

    def sim(self, mcell, shape, stepFunction = None):

        self.debug()()      
        
        cell = mcell.cell
        # stimuli = self.create_stimuli(cell, 3)

        recordings = {}

        recordings['time'] = self.hocObj.Vector()
        recordings['soma(0.5)'] = self.hocObj.Vector()
        # recordings['isoma(0.5)'] = self.hocObj.Vector()

        recordings['time'].record(self.hocObj._ref_t, 0.025)
        # recordings['soma(0.5)'].record(init.CellLoader.cell.soma[0](0.5)._ref_v, 0.1)
        recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, 0.025)
        # recordings['isoma(0.5)'].record(cell.soma[0](0.5)._ref_i_membrane, 0.1)

        self.hocObj.tstop = 3000 #TODO: make parameter
        self.hocObj.dt = 0.025 # 0.025 -> 25 us = 40kHz sampling
        # create_stimuli(init.CellLoader.cell, 3)

        # assert(sys.getrefcount(stimuli) > 1)
        memireclist = self.hocObj.List()
        for sec in mcell.sections:
            for seg in sec:

                memirec = self.hocObj.Vector()
                memirec.record(seg._ref_i_membrane, .1)
                memireclist.append(memirec)

        for sec in self.hocObj.allsec():
            if self.hocObj.ismembrane("xtra", sec=sec):
                for x in sec:
                    x.es_xtra = 1


        self.hocObj._ref_stim_xtra[0] = 1

        if cell.synapses:
            print("synapses are working!!!!")
            # self.microStim(0,0,0,0, 0.3)
            for i in range(54): #TODO: paramaterize these synapse indices to make sure we know what we initialize
                cell.synapses.active_pre_mtypes.x[i] = 1

            cell.synapses.update_synapses(self.hocObj.Shape())

        print('Disabling variable timestep integration')
        self.hocObj.cvode_active(0)

        self.hocObj.finitialize(-65)
        self.hocObj.fcurrent()
        print('Running for %f ms' % self.hocObj.tstop)

        counter = 0

        self.hocObj.t = 0
        interval = 1


        while self.hocObj.t < self.hocObj.tstop:

            self.hocObj.fadvance()
            #if self.hocObj.t > 10:
            #    self.hocObj._ref_stim_xtra[0] = -50


            if counter == interval:
                print(self.hocObj.t, cell.soma[0](
                    0.5)._ref_v[0], cell.soma[0](0.5)._ref_i_membrane[0])
                counter = 0
                imem = np.array(memireclist)[:,-1]
                
                if self.renderPipe:
                 
                  self.renderPipe.send(imem)
                
                if self.mapping is not None:
                    lfp = np.matmul(self.mapping, imem).reshape(shape)
                    if stepFunction:
                        cmd = stepFunction(lfp)

                        self.handleCmd(cmd)
                """ 
                plt.clf()
                plt.title(self.title)
                plt.plot(recordings['time'], recordings['soma(0.5)'])
                plt.pause(0.05)
                """

            counter += 1

        self.celldata.imem= np.array(memireclist)

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
    def handleCmd(self, cmd):
        print("basicsim -- " + str(cmd))
        cmd.DO(self)

    def constStim(self, theta, phi):
        theta = theta * np.pi / 180
        phi = phi * np.pi / 180
        Ex = np.sin(theta) * np.cos(phi)
        Ey = np.sin(theta) * np.sin(phi)
        Ez = np.cos(theta)

        for sec in self.hocObj.allsec():
            if self.hocObj.ismembrane("xtra", sec=sec):
                for seg in sec:
                    seg.es_xtra = -50 * \
                        (-Ex * seg.x_xtra + Ey * seg.y_xtra + Ez * seg.z_xtra) * 1e-3


    def create_stimuli(self, cell, step_number):
        """Create the stimuli"""

        print('Attaching stimulus electrodes')

        stimuli = []
        step_amp = [0] * 3

        with open('current_amps.dat', 'r') as current_amps_file:
            first_line = current_amps_file.read().split('\n')[0].strip()
            hyp_amp, step_amp[0], step_amp[1], step_amp[2] = first_line.split(
                ' ')

        iclamp = self.hocObj.IClamp(0.5, sec=cell.soma[0])
        iclamp.delay = 700
        iclamp.dur = 2000
        iclamp.amp = float(step_amp[step_number - 1])
        print('Setting up step current clamp: '
              'amp=%f nA, delay=%f ms, duration=%f ms' %
              (iclamp.amp, iclamp.delay, iclamp.dur))

        stimuli.append(iclamp)

        hyp_iclamp = self.hocObj.IClamp(0.5, sec=cell.soma[0])
        hyp_iclamp.delay = 0
        hyp_iclamp.dur = 3000
        hyp_iclamp.amp = float(hyp_amp)
        print('Setting up hypamp current clamp: '
              'amp=%f nA, delay=%f ms, duration=%f ms' %
              (hyp_iclamp.amp, hyp_iclamp.delay, hyp_iclamp.dur))

        stimuli.append(hyp_iclamp)

        return stimuli


def startRenderer(celldata, elec):
    parent_conn, conn = Pipe()
    print("rendering")
    p = Process(target=plot3d.start, args=(celldata, conn, elec))
    p.start()
    print("unrendering")
    return parent_conn


if __name__ == '__main__':
    mcell = MyelinatedCell()

    for arg in sys.argv[1:]:
        assert(arg in ["-noaxon", "-nosynapse"])

    myelinate_ax = "-noaxon" not in sys.argv
    synapse = '-nosynapse' not in sys.argv

    mcell.loadcell('L5_TTPC2_cADpyr232_1', myelinate_ax=myelinate_ax, synapses=synapse)
    neuron.h.finitialize(-65)
    celld = CellData(mcell, neuron.h)
    sim = basicsim(celld, neuron.h)


    # sim.renderPipe = startRenderer(celld) #TODO?: include dummy elec param here?

    res = 10
    xmin = min(celld.xend)

    ymin = int(min(celld.yend) / res) - 5
    ymax = int(max(celld.yend) / res) + 5

    zmin = int(min(celld.zend) / res) - 5
    zmax = int(max(celld.zend) / res) + 5

    X, Y, Z = np.mgrid[1:2, ymin - 50:ymax + 50:1, zmin - 5:zmin + 5:1] * res
    X = X - xmin - 10

    sim.mapping = np.zeros([len(X.flatten()), celld.totseg])
    sim.doDebug = True

    for i, (x, y, z) in enumerate(zip(X.flatten(), Y.flatten(), Z.flatten())):
        sim.mapping[i, :] = lfpcalc.calc_lfp_linesource(
            celld, x, y, z, 0.3, celld.diam)
        # print(x,y,z, mapping.shape)

    sim.sim(mcell, X.shape[1:])

    with open('testpickle.celldata', 'wb') as pfile:
        pickle.dump(celld, pfile)
