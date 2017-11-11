## @file miniGUI.py
#  Defines mini GUI class for drawing molecules from the command line
#
#  Written by Terry Gani for HJK Group
#
#  Dpt of Chemical Engineering, MIT

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from molSimplify.Classes.mWidgets import *
from molSimplify.Classes.globalvars import *
import sys, os, random, shutil, unicodedata, inspect, glob, time, tempfile
from pkg_resources import resource_filename, Requirement
import xml.etree.ElementTree as ET

## Mini GUI class for drawing molecules from the command line
class miniGUI(QApplication):
    def __init__(self, args):
        QApplication.__init__(self,args)
        p = QPalette()
        p.setColor(QPalette.Background,Qt.white)
        self.window = mQMainWindow()
        self.lwindow = QWidget()
        self.lwindow.setPalette(p)
        self.lgrid = QGridLayout()
        self.lwindow.setLayout(self.lgrid)
        self.lwindow.setWindowTitle('SVG Drawing')

    #def qcloseligs(self):
        self.lwindow.hide()

    def addsvg(self,filename):
        self.svgwidget = mSvgWidget(filename)
        self.lgrid.addWidget(self.svgwidget,0,1)
        #self.lwclose = QPushButton('Close')
        #self.lwclose.clicked.connect(self.qcloseligs)
        #self.lgrid.addWidget(self.lwclose,1,0)

    def show(self):
        self.lwindow.show()
        sys.exit(self.exec_())
