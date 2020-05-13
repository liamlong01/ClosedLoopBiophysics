
import sys
import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtGui as QtGui
from PyQt5.QtCore import QThread
from FitzNagNeuron import FitzNagNeuron
from CoupledOscillator import CoupledOscillator, noisyCO
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io as io


class simGUI(QtWidgets.QWidget):


    def __init__(self):
        super().__init__()

        self.initUI()


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ex = simGUI()
    sys.exit(app.exec_()