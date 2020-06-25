from multiprocessing import Process, Pipe

from aberraAxon import MyelinatedCell
import LFPylite
from ElectrodeArray import ElectrodeArray

import os
import pdb
import numpy as np
import traceback
import time
import lfpcalc

import cv2 as cv

import matplotlib.pyplot as plt


class Message():
    def Do(self, singlecellNEURON):
        pass

    def __str__(self):
        return "Message"


class TestMessage(Message):
    def DO(self, singlecellNEURON):
        print("we are running test message")


class CreateMessage(Message):
    def __init__(self, msg):
        self.msg = msg


class GeometryMessage(Message):

    def __init__(self, projection):
        self.projection = projection

    def DO(self, cell):
        print("we are running geometry message")

        cell.conn.send("hello mom")

        # singlecellNEURON.get_pt3d_polygons(projection=self.projection)


class StopMessage(Message):
    def DO(self, singlecellNEURON):
        print("we are running stop message")
        singlecellNEURON.__stopped = True


class cell_init_msg(Message):

    def __init__(self, hoccommand):
        self.hoccommand = hoccommand

    def DO(self, ncell):
        ncell.hoc(self.hoccommand)
        print("reloaded again!!!!")

class eval_msg(Message):

    def __init__(self, command):
        self.command = command

    def DO(self, ncell):
        eval(self.command)

class generic_func_msg(Message):

    def __init__(self, func, cellarg=True):
        self.func = func
        self.cellarg = cellarg

    def DO(self, ncell):
        if self.cellarg:
            self.func(ncell)
        else:
            self.func()

class sim_msg(Message):

    def DO(self, ncell, render=False):
        ncell.simulator.sim(ncell.mcell, ncell.shape)
        if render: 
            ncell.simulator.renderPipe = LFPylite.startRenderer(ncell.celld, ncell.elec)

class sim_msg(Message):

    def DO(self, ncell):
        ncell.simulator.sim(ncell.mcell, ncell.shape, ncell.stepFunction)


class add_elec_msg(Message):
    def __init__(self, X, Y, Z):
        self.X = X
        self.Y = Y
        self.Z = Z

    def DO(self, ncell):
        X, Y, Z = self.X, self.Y, self.Z
        ncell.shiftcell(X, Y, Z)

        ncell.simulator.mapping = np.zeros([len(X.flatten()), ncell.celld.totseg])

        for i, (x, y, z) in enumerate(zip(X.flatten(), Y.flatten(), Z.flatten())):

            ncell.simulator.mapping[i, :] = lfpcalc.calc_lfp_linesource(
                ncell.celld, x, y, z, 0.3, ncell.celld.diam)
            # print(x,y,z, mapping.shape)
            ncell.shape = X.shape[1:]

        ncell.elec = ElectrodeArray(X,Y,Z)

       
        print("initialized the render")

class Terminate(Exception):
    pass

class terminate_msg(Message):

    def DO(self, ncell):
        raise Terminate()

class SetStim(Message):
    def __init__(self, current, x,y,z,sigma):
        self.x = x
        self.y = y
        self.z = z
        self.sigma = sigma
        self.current = current
  
    def DO(self, simcontrol):
        simcontrol.microStim(self. current, self.x, self.y, self.z, self.sigma)


class CLStim():
    sigma = 0.3
    current = 100e-6

    def __init__(self, network, shape, dt, label="Ulabeled CLStim"):
        """
        network is a collection of NetworkCells
        """
        self.trigger = 0
        self.shape = shape
        self.pipes = []  # DO NOT send to pipes
        for conn in network:
            self.pipes.append(conn)

        self.parent_conn, self.conn = Pipe()
        self.label = label
        self.clims = [-100e-6, 100e-6]

        self.time = 0
        self.dt = dt

        self.stimming = 0 
        self.currents = []
        self.current = 0

        p = Process(target=self.Core, args=(self.pipes, self.conn))
        p.start()

    def Core(self, pipes, conn):
       
        try:

            while True:
                LFPs = np.zeros(self.shape)
                for pipe in pipes:

                    lfp = pipe.recv()
                  
                    LFPs = np.add(LFPs, lfp.reshape(self.shape))
                    LFPs = self.addNoise(LFPs)

                self.trigger = 1

                if self.time > 1e-3: 
                    self.plot(LFPs)
                    
                    
                self.biphasicpulse()
                self.currents = np.append(self.currents, self.current)
                print(self.time, "s, current setting:", self.current, "maxlfp", np.max(LFPs))
                self.current = 0
                for pipe in pipes:
                   pipe.send(SetStim(self.current, x=0,y=0,z=0,sigma=0.3)) # send stim

                 

                self.time+=self.dt

        except Exception as e:
            self.print_info("Unhandled exception:" +
                            str(e.__class__) + ' - ' + str(e) + "- terminating process")
            traceback.print_exc()
            return e

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


            elif self.time - self.start < 10e-3:
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
        plt.subplot(2,1,1)
        im = plt.imshow(lfps, aspect='auto')
        cbar = plt.colorbar()
        im.set_clim(self.clims[0], self.clims[1])
        plt.title("LFPS")
        # self.csd(lfps)
        
        plt.subplot(2,1,2)
        plt.plot(self.currents)
        plt.plot()
        
        plt.pause(0.05)
        

    def csd(self, lfps):
        csd = cv.Laplacian(lfps, cv.CV_64F, ksize=3)
        plt.subplot(1,2,2)
        plt.clf()
        im = plt.imshow(csd, aspect='auto')
        cbar = plt.colorbar()
        plt.title("CSD")
       
        plt.pause(0.05)


    def addNoise(self, matrix):


        rms = self.calcNoiseRMS()
        matrix = np.add(matrix, np.random.normal(scale = rms, size=matrix.shape))
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
        factor = 20 # a conservative increase in capacitance (some say 100times cacitance)
        delt_Vpad = np.sqrt(kB * T / (Cpad * factor))


        fAP = 30e3
        delt_Vlead = np.sqrt(4*kB*T*rlead*fAP)

        return delt_Vpad+delt_Vlead #ignoring spreading resistance and degradation of coupling

    def print_info(self, info):
        print(self.label + '----' + info)


class NetworkCell():

    def __init__(self, cell, label="Unlabeled"):
        self.cellstr = cell
        self.parent_conn, self.conn = Pipe()

        self.label = label + cell
        try:
            p = Process(target=self.Core, args=(self.conn,))
            p.start()
        except BrokenPipeError:  # TODO (liam): Figure out why this occurs
            self.print_info("Broken Pipe Error occured")
        self.__stopped = False

    def print_info(self, info):
        print(self.label + '----' + info)

    def Core(self, conn):
        try:
            # initialize cell here
            # we need to delay these imports until we are in a new process
            # this will keep hocObj's unique
            # ODO move these imports? doing mechanism loading in bash now

            import LFPy
            import neuron

            import sys
            sys.stdout = open("logging/" + self.label + ".txt", "w")
            self.print_info("loading nrn")
            self.hoc = neuron.h
            self.conn = conn
            self.print_info("sending msg")
            try:
                self.conn.send(CreateMessage(self.label))
            except BrokenPipeError:
                return

            # cwd = os.getcwdnce            # os.chdir(neuronfolder)
            self.print_info("Instatiating cell")
            self.mcell = MyelinatedCell(hocObj=self.hoc)
            self.mcell.loadcell(self.cellstr, myelinate_ax=True, synapses=False)

            self.celld = LFPylite.CellData(self.mcell, self.hoc)

            self.simulator = LFPylite.basicsim(
                self.celld, self.hoc, title=self.label)

            self.__stopped = False

            while not self.__stopped:
                Msg = self.conn.recv()
                try:
                    self.handleMsg(Msg)
                except Terminate as e:
                    self.print_info(str(e) + ' occurred in ' + str(Msg))
                    if e.__class__ == Terminate:
                        return

        except Exception as e:
            self.print_info("Unhandled exception:" +
                            str(e.__class__) + ' - ' + str(e) + "- terminating process")
            traceback.print_exc()
            return e

    def handleMsg(self, Msg):

        self.print_info("running msg.do")
        Msg.DO(self)

        self.print_info(str(Msg))

    def stepFunction(self, lfp):
        if self.conn.poll():
            self.handleMsg(self.conn.recv())
        self.conn.send(lfp)

        msg = self.conn.recv()
        self.print_info(str(msg))
        return msg

    def shiftcell(self, x, y, z):
        xshift = - np.min(self.celld.xstart) +np.mean(x)
        self.print_info("shifting cell by {} units".format(xshift))
        self.celld.shift(xshift, 0, 0)
