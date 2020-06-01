from multiprocessing import Process, Pipe
from aberraAxon import MyelinatedCell
import LFPylite
import os
import pdb
import numpy as np

import lfpcalc


class Message():
    def Do(self, singlecellNEURON):
        pass


class TestMessage(Message):
    def DO(self, singlecellNEURON):
        print("we are running test message")


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

    def DO(self, ncell):
        ncell.sim.sim(ncell.mcell)


class Terminate(Exception):
    pass

class terminate_msg(Message):

    def DO(self, ncell):
        raise Terminate()


class NetworkCell():

    def __init__(self, cell, label="UnlabeledNetworkCell"):
        self.cellstr = cell
        self.parent_conn, self.conn = Pipe()

        self.label = label
        p = Process(target=self.Core, args=(self.conn,))
        p.start()

        self.__stopped = False

    def print_info(self, info):
        print(self.label + '----' + info)

    def Core(self, conn):
        try:
            # initialize cell here
            # we need to delay these imports until we get the correct mechanism file
            # ODO move these imports? doing mechanism loading in bash now

            import LFPy
            import neuron
            self.hoc = neuron.h
            self.conn = conn

            # cwd = os.getcwd()
            # os.chdir(neuronfolder)
            self.print_info("Instatiating cell")
            self.mcell = MyelinatedCell(hocObj=self.hoc)
            self.mcell.loadcell(16, myelinate_ax=True, synapses=True)
            self.celld = LFPylite.CellData(self.mcell)
            self.sim = LFPylite.basicsim(self.celld)
            res = 10

            xmin = min(self.celld.xend)

            ymin = int(min(self.celld.yend) / res) - 5
            ymax = int(max(self.celld.yend) / res) + 5

            zmin = int(min(self.celld.zend) / res) - 5
            zmax = int(max(self.celld.zend) / res) + 5

            X, Y, Z = np.mgrid[1:2, -50:50:1, -5:5:1] * res
            X = X - xmin
            X = X.flatten()
            Y = Y.flatten()
            Z = Z.flatten()

            self.sim.mapping = np.zeros([self.celld.totseg, X.shape[0]])

            for i, (x, y, z) in enumerate(zip(X, Y, Z)):
                # print(i,x,y,z)
                self.sim.mapping[:, i] = lfpcalc.calc_lfp_linesource(
                    self.celld, x, y, z, 0.3, self.celld.diam)
                # print(x,y,z, mapping.shape)

            self.__stopped = False

            while not self.__stopped:
                Msg = conn.recv()
                try:
                    self.print_info("running msg.do")
                    Msg.DO(self)

                    self.print_info(str(Msg))
                except Terminate as e:
                    self.print_info(str(e) + ' occurred in '+ Msg)
                    if e.__class__ == Terminate:
                        return

        except Exception as e:
            self.print_info("Unhandled exception:" +
                            str(e.__class__) +' - ' + str(e) + "- terminating process")
            return e
