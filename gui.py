
import sys
import os

import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtGui as QtGui
from PyQt5.QtCore import QThread

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io as io
import pdb


from importlib import reload

import singlecellNEURON
import LFPylite


class cellcontrol():

    def __init__(self, cell_directory):
        self.cell_directory = cell_directory
        self.cells = [f.name for f in os.scandir(cell_directory) if f.is_dir()]

        self._layers()

    def gendirectory(self):

        return self.layerstr + '_' + self.morphostr + '_' + self.electstr + '_' + self.idxstr

    def _layers(self):
        layer_dict = {}
        for cell in self.cells:
            layer_key = cell.split('_')[0]
            if layer_key in layer_dict.keys():
                layer_dict[layer_key].append(cell)
            else:
                layer_dict[layer_key] = [cell]

        self.layers = layer_dict

    def _morpho(self, layerstr):
        morpho_dict = {}
        self.layerstr = layerstr

        for cell in self.cells:

            if cell.split('_')[0] == layerstr:

                morpho_key = cell.split('_')[1]
                if morpho_key in morpho_dict.keys():
                    morpho_dict[morpho_key].append(cell)
                else:
                    morpho_dict[morpho_key] = [cell]

        return morpho_dict

    def _electtype(self, layerstr, morphostr):
        elect_dict = {}
        self.morphostr = morphostr
        for cell in self.cells:

            if cell.split('_')[0] == layerstr \
                    and cell.split('_')[1] == morphostr:

                elect_key = cell.split('_')[2]
                if elect_key in elect_dict.keys():
                    elect_dict[elect_key].append(cell)
                else:
                    elect_dict[elect_key] = [cell]

        return elect_dict

    def _cellidx(self, layerstr, morphostr, electstr):
        morpho_dict = {}
        self.electstr = electstr
        for cell in self.cells:
            if cell.split('_')[0] == layerstr \
                and cell.split('_')[1] == morphostr \
                    and cell.split('_')[2] == electstr:

                morpho_key = cell.split('_')[3]
                if morpho_key in morpho_dict.keys():
                    morpho_dict[morpho_key].append(cell)
                else:
                    morpho_dict[morpho_key] = [cell]

        return morpho_dict

    def setidx(self, toWhat):
        self.idxstr = toWhat


class simGUI(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()
        self.cellcontrol = cellcontrol('neurons/')
        self.conns = []
        self.initUI()

    def initUI(self):
        layout = QtWidgets.QHBoxLayout()

        self.setLayout(layout)

        layout.addWidget(self.initCellSelector('Select Cells'))

        layout.addWidget(self.cellBuilder('CellNetwork'))
        
        layout.addWidget(self.miscControls("Misc. Controls"))

        self.setGeometry(300, 300, 250, 150)
        self.setWindowTitle('Neuron Simulations')

        self.show()

    def cellBuilder(self, title):
        groupbox = QtWidgets.QGroupBox(title, parent=self)

        layout = QtWidgets.QVBoxLayout()
        groupbox.setLayout(layout) 

        self.cellListView = QtWidgets.QListWidget(parent=groupbox)
        layout.addWidget(self.cellListView)

        addbtn = QtWidgets.QPushButton('Add cell', groupbox)
        addbtn.clicked.connect(lambda: self.addCell(self.cellcontrol.gendirectory()))
        layout.addWidget(addbtn)

        return groupbox



    def miscControls(self, title):
        groupbox = QtWidgets.QGroupBox(title, parent=self)

        layout = QtWidgets.QVBoxLayout()
        groupbox.setLayout(layout)

        self.renderCbox = QtWidgets.QCheckBox("Render", parent=self)
        layout.addWidget(self.renderCbox)
        
        btn3 = QtWidgets.QPushButton('Try LFP sim', self)
        btn3.clicked.connect(self.simcells)
        btn3.resize(btn3.sizeHint())
        layout.addWidget(btn3)

        btn = QtWidgets.QPushButton('Kill All', self)
        btn.clicked.connect(lambda: self.terminateAll())
        btn.resize(btn.sizeHint())
        layout.addWidget(btn)

        return groupbox

    def terminateAll(self):
        for conn in self.conns:
            conn.send(singlecellNEURON.terminate_msg())
        self.conns = []
        self.cellListView.clear()
        reload(singlecellNEURON)

    def addCell(self, cell_str):
        newConn, newLabel = self.createNetworkCell(cell_str)

        self.conns.append(newConn)

        QtWidgets.QListWidgetItem(newLabel, parent= self.cellListView)

    def simcells(self):
        res = 10
        X, Y, Z = np.mgrid[1:2, -100:100:1, -20:+20:1] * res

        self.cls = singlecellNEURON.CLStim(
            self.conns, X.shape[1:], dt=25e-6)

        for conn in self.conns:
            conn.send(singlecellNEURON.add_elec_msg(X,Y,Z))
            conn.send(singlecellNEURON.sim_msg(self.renderCbox.checkState()))

    def createNetworkCell(self, cell_str):
        """
        wrapper for network cell initializer that accounts
        for neuron crashing on import for subprocesses
        """
        loaded = False
        while not loaded:
            cell = singlecellNEURON.NetworkCell(cell_str)
            conn = cell.parent_conn
            print("waiting for neuron import ....")
            try:
                if conn.poll(5):  # check if msg received after initializing
                    print("got pipe ....")
                    loaded = True
            except BrokenPipeError:  # process crashed
                pass

        conn.recv()  # clear once something is received

        return conn, cell.label

    def initNEURON(self, cell_str):
        #TODO: remove - currently here in short term for reference
        print("Parent - initializing child...")
        
        res = 10
        self.X, self.Y, self.Z = np.mgrid[1:2, -50:50:1, -5:+5:1] * res

        self.pipe,l = self.createNetworkCell(cell_str)
        self.pipe.send(singlecellNEURON.add_elec_msg(X,Y,Z))

        self.pipe2,l = self.createNetworkCell("L4_ChC_cACint209_1")
        self.pipe2.send(singlecellNEURON.add_elec_msg(X,Y,Z))

        self.cls = singlecellNEURON.CLStim(
            [self.pipe, self.pipe2], X.shape[1:])

        print("Parent - sending...")
        self.pipe.send(singlecellNEURON.cell_init_msg("create axon"))
        self.pipe.send(singlecellNEURON.eval_msg("ncell.hoc.axon"))

        print("Parent - receiving...")

    def initCellSelector(self, title):
        groupbox = QtWidgets.QGroupBox(title, parent=self)
        dv = QtGui.QDoubleValidator(parent=groupbox)
        layout = QtWidgets.QFormLayout(parent=groupbox)

        nbox = QtWidgets.QLineEdit(str('hello'), parent=self)
        nbox.setValidator(dv)
        layout.addRow('noise level: ', nbox)

        Lcombobox = QtWidgets.QComboBox(parent=self)
        Lcombobox.addItems(self.cellcontrol.layers.keys())
        Lcombobox.setCurrentText("L1")
        layout.addRow('layer: ', Lcombobox)

        self.morphcombobox = QtWidgets.QComboBox(parent=self)
        self.morphcombobox.addItems(self.cellcontrol._morpho("L1"))
        layout.addRow('morphology: ', self.morphcombobox)

        Lcombobox.currentTextChanged.connect(self.updatemorpho)

        self.electcombobox = QtWidgets.QComboBox(parent=self)
        self.electcombobox.addItems(self.cellcontrol._electtype("L1", "DAC"))
        layout.addRow('electrictype: ', self.electcombobox)

        self.morphcombobox.currentTextChanged.connect(
            lambda x: self.updateelect(Lcombobox.currentText(), x))

        self.idxcombobox = QtWidgets.QComboBox(parent=self)
        self.idxcombobox.addItems(
            self.cellcontrol._cellidx("L1", "DAC", "bNAC219"))
        layout.addRow('idx: ', self.idxcombobox)

        self.electcombobox.currentTextChanged.connect(
            lambda x: self.updateidx(Lcombobox.currentText(),
                                     self.morphcombobox.currentText(),
                                     x))

        self.idxcombobox.currentTextChanged.connect(self.cellcontrol.setidx)
        self.cellcontrol.setidx('1')
        btn = QtWidgets.QPushButton('Plot Potential', self)
        btn.clicked.connect(self.debug)
        btn.resize(btn.sizeHint())
        layout.addRow('Debug func', btn)

        return groupbox

    def updatemorpho(self, layer_str):
        while self.morphcombobox.count() > 0:
            self.morphcombobox.removeItem(0)
        self.morphcombobox.addItems(self.cellcontrol._morpho(layer_str).keys())

    def updateelect(self, layer_str, morpho_str):
        while self.electcombobox.count() > 0:
            self.electcombobox.removeItem(0)
        self.electcombobox.addItems(
            self.cellcontrol._electtype(layer_str, morpho_str).keys())

    def updateidx(self, layer_str, morpho_str, elect_str):
        while self.idxcombobox.count() > 0:
            self.idxcombobox.removeItem(0)
        self.idxcombobox.addItems(self.cellcontrol._cellidx(
            layer_str, morpho_str, elect_str).keys())

    def debug(self):
        pdb.set_trace()


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ex = simGUI()
    sys.exit(app.exec_())
