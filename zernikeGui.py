# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 37:49:19 2014

@author: Nathan
"""
import os
import platform
import sys
import time
import ast

import Zernike as zernike
import Dither as dither
import numpy as np
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.Qwt5 import *

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs 

from Dither import *

import numpy
from numpy import arange, sin, pi

class Zplot(QMainWindow):
    #Main window containing all widgits

    def __init__(self,parent=None):
        super(Zplot,self).__init__(parent)
        self.t = 0        
        self.calcTime = 0
        self.zBoxes = 37*[0]
        self.zBoxLabels = 37*[0]
        self.cofs = 37*[0]
        self.wavelength = 633
        self.voxelSize = 20
        self.diameter = 50
        self.xRes = self.yRes = int(self.diameter*1000/self.voxelSize)
        self.height = 13
        self.voxelsize = 20000
        self.n1 = 1.5-0.021
        self.n2 = 1.5+0.021
        self.is3d = False
        self.isCalculated = False
        self.isDithered = False
        self.isLayersGenerated = False
                
        
        #Status Bar Messages        
        self.calcMessage = "CALCULATING..."
        self.readyMessage = "READY"
        
        
        
        #Draw Everything
        self.make_widgets()  
        self.set_geometry()
        
    
    
    def save_system(self):
        dialog = QFileDialog(self)
        filename = dialog.getSaveFileName()
        f = open(filename,'w')
        
        self.params = [self.cofs, self.wavelength, self.height, self.n1, 
                       self.n2, self.xRes, self.yRes,self.voxelsize]
        
        
        for i in self.params:
            f.write(str(i)+"\n")        
        
        """if self.isCalculated:
            np.savetxt(f,self.z.opd)
            
        if self.isDithered:
            self.params.append(self.d.dithered)
        
        if self.isLayersGenerated:
            self.params.append(self.d.data)"""
            
        
        
                      
        
    def load_system(self):
        dialog = QFileDialog(self)
        filename = dialog.getOpenFileName()
        f = open(filename,'r')
        self.cofs = ast.literal_eval(f.readline())
        self.wavelength = float(f.readline())
        self.height = int(f.readline())
        self.n1 = float(f.readline())
        self.n2 = float(f.readline())
        self.xRes = int(f.readline())
        self.yRes = int(f.readline())
        self.voxelsize = int(f.readline())
        self.update_boxes()
        print f.readlines()

    
    def update_values(self):
        """Queries all text boxes and updates the system parameters"""
        self.voxelSize = float(self.voxelSizeBox.text())
        self.voxelsize=self.voxelSize*1000
        self.diameter = float(self.diameterBox.text())
        self.xRes = self.yRes = int(self.diameter*1000/self.voxelSize)
        self.wavelength = float(self.wavelengthBox.text())
        self.height = int(self.heightBox.text())
        self.n1 = float(self.n1Box.text())
        self.n2 = float(self.n2Box.text())
        
    def update_boxes(self):
        for i in range(len(self.zBoxes)):
            self.zBoxes[i].setText(str(self.cofs[i]))
        self.wavelengthBox.setText(str(self.wavelength))
        self.heightBox.setText(str(self.height))
        self.n1Box.setText(str(self.n1))
        self.n2Box.setText(str(self.n2))
        self.voxelSizeBox.setText(str(self.voxelSize))
        self.diameterBox.setText(str(self.diameter))
        
    
    def toggle(self, state):
        """Assigns bool 'state' to the value of the setEnabled() parameter 
        for each button"""
        self.buttonCalc.setEnabled(state)
        self.buttonDither.setEnabled(state)
        #self.button2D.setEnabled(state)
        if self.isDithered:
            self.buttonLayers.setEnabled(state)
            self.buttonLayersFast.setEnabled(state)
            self.layersMenu.setEnabled(state)
        if self.isLayersGenerated:
            self.buttonExportAll.setEnabled(state)
        #self.buttonSlice.setEnabled(state)
        if not state:
            self.t = time.time()
            self.status.showMessage(self.calcMessage)
        else:
            self.calcTime = round(time.time()-self.t,5)
            self.readyMessage ="READY       Calculation Time: "+str(self.calcTime)+"s"
            self.status.showMessage(self.readyMessage)
        #self.button3D.setEnabled(state)
        
        
    def export_all(self):
        """Opens a dialog and asks the user to select a folder to save bitmaps in"""
        self.layerFolder = str(QFileDialog.getExistingDirectory())
        self.d.savelayers(self.layerFolder)
        return 1
    

    def __generate_zernike_axes(self):
        """Clears the previous figure and creates new axes for the ideal wavefront 
        plot defined by the zernike coefficients"""
        self.zernikeFig.clear()
        self.zernikeGridSpec = gs.GridSpec(1,2,width_ratios = [1,.1])
        self.zernikeAxes1 = self.zernikeFig.add_subplot(self.zernikeGridSpec[0])
        self.zernikeAxes2 = self.zernikeFig.add_subplot(self.zernikeGridSpec[1])
        self.zernikeAxes1.set_title(str(self.diameter)+" mm GRIN Lens")
        self.zernikeAxes1.set_xlabel("GRIN Lens Width (mm)")
        self.zernikeAxes2.set_ylabel("Grin Lens Height (mm)")
        self.zernikeAxes2.set_xlabel("OPD ($\mu$m)")
        
    def __generate_dither_axes(self):
        """Same as above function, but for the dither plot"""
        self.ditherFig.clear()
        self.ditherGridSpec = gs.GridSpec(1,2,width_ratios = [1,.1])
        self.ditherAxes1 = self.ditherFig.add_subplot(self.ditherGridSpec[0])
        self.ditherAxes2 = self.ditherFig.add_subplot(self.ditherGridSpec[1])
        self.ditherAxes1.set_title(str(self.diameter)+" mm GRIN Lens Dithered")
        self.ditherAxes1.set_xlabel('GRIN Lens Width (mm)')
        self.ditherAxes2.set_xlabel('Effective n')
        
    def __generate_3d_axes(self):
        self.zernikeFig.clear()
        self.zernikeAxes3d = self.zernikeFig.add_subplot(111,projection = '3d')
        self.zernikeAxes3d.set_xlabel("GRIN Lens Width (mm)")
        self.zernikeAxes3d.set_zlabel("OPD ($\mu$m)")
        
    """def __generate_slice_axes(self):
        #Same as above function, but for the dither plot
        self.sliceFig.clear()
        self.sliceGridSpec = gs.GridSpec(1,2,width_ratios = [1,.1])
        self.sliceAxes1 = self.sliceFig.add_subplot(self.sliceGridSpec[0])
        self.sliceAxes2 = self.sliceFig.add_subplot(self.sliceGridSpec[1])
        self.sliceAxes2.set_title('Effective n')"""
    
    def generate_layers(self):
        """calls the 'calc3d()' function from 'dither.py' to generate dithered layers
        that stack up to an effective index of refraction. Also creates the menu
        that allows the user to select a layer to view"""
        self.update_values()
        self.isLayersGenerated = True        
        self.toggle(0)
        self.d.calc3d()
        self.layersMenu.setEnabled(True)
        self.layersMenu.clear()
        self.layersMenu.addItem("Composite")
        for i in range(self.height):
            self.layersMenu.addItem("Layer %g" %(i+1))
            
        self.isLayersGenerated = True
        self.toggle(1)
    
    def generate_layers_fast(self):
        """Same as 'generate_layers', but with the 'calc3dfast()' function from 'dither.py'"""
        
        self.update_values()
        self.toggle(0)
        self.d.calc3dfast()
        self.layersMenu.setEnabled(True)
        self.layersMenu.clear()
        self.layersMenu.addItem("Composite")
        for i in range(self.height):
            self.layersMenu.addItem("Layer %g" %(i+1))
            
        self.isLayersGenerated = True
        self.toggle(1)
        
        
        
        
        
    def show_dithered_layer(self, layer):
        try:
            if layer[-1] == 'e':
                self.__generate_dither_axes()
                image = self.ditherAxes1.imshow(np.float32(self.d.dithered))
                self.ditherFig.colorbar(image, self.ditherAxes2)
                self.ditherCanvas.draw()
            else:
                if layer[-2] == ' ':
                    index = int(layer[-1])-1
                else:
                    index = 10*int(layer[-2])+int(layer[-1])-1
                
                self.ditherFig.clear()
                self.ditherAxes = self.ditherFig.add_subplot(111)
                self.ditherAxes.imshow(self.d.data[:,:,index])
                self.ditherCanvas.draw()
            
        except:
            return 1
    
        
    def dither(self):
        self.toggle(0)
        self.update_values()
        self.isDithered = True
        self.status.showMessage(self.calcMessage)
        
        #try:
        self.d = dither.Dither(self.z.opd, self.height,self.voxelsize,self.n1,self.n2)                 
        
        """except:        
            self.d = dither.Dither(self.z.opd,
                               self.height,self.voxelsize,
                               self.n1,self.n2)"""
        
        print self.d.height
        self.status.showMessage(self.readyMessage)
        
        self.__generate_dither_axes()
        image = self.ditherAxes1.imshow(np.float32(self.d.dithered),extent = (0,self.diameter,0,self.diameter))
        self.ditherFig.colorbar(image, self.ditherAxes2)
        self.ditherCanvas.draw()
        self.toggle(True)
        self.buttonLayers.setEnabled(True)
        self.buttonLayersFast.setEnabled(True)
   
    """def slice_layers(self):
        self.toggle(0)
        self.update_values()   
        self.s = sl.Slice(self.z.opd,self.height,self.voxelsize,self.n1,self.n2)
        self.__generate_slice_axes()
        image = self.sliceAxes1.imshow(self.s.total)
        self.sliceFig.colorbar(image,self.sliceAxes2)
        self.sliceCanvas.draw()
        self.toggle(1)"""
    
    def minheight(self,opd,pvheight,n1,n2):
        '''
        Given the desired optical path difference (opd), the height of each voxel (vheight),
        and the two indices of refraction, this function with tell you the minimum needed layers
        that are needed to acheive the desired opd.
        '''
        
        if n1>n2:
            raise ValueError('n2 must be greater than n1')
        minheight=opd/(pvheight*(n2-n1))
        return minheight    
    
    
    def enable_3d(self, state):
        if state == 0:
            self.dimensionMenu.setEnabled(False)
        else:
            self.dimensionMenu.setEnabled(True)
        
     
    def plot_zernike(self, dimension):
        if self.isCalculated:
            if not dimension:
                self.__generate_zernike_axes()
                image = self.zernikeAxes1.imshow(self.z.opd/1000,cmap='jet',extent=[0, self.diameter, 0, self.diameter])
                self.zernikeFig.colorbar(image, self.zernikeAxes2)        
                self.zernikeFig.hold(False)
                self.zernikeCanvas.draw()
            else:
                self.__generate_3d_axes()
                self.zernikeAxes3d.plot_surface(self.z.x,self.z.y,self.z.opd,cstride=self.z.opd.shape[0]//50,rstride=self.z.opd.shape[1]//50,linewidth=0)
                self.zernikeFig.hold(False)
                self.zernikeCanvas.draw()
        
        
        
        
        
        
    def calculate(self):
        self.toggle(0)
        self.update_values()
        
        for i in range(37):
            try:
                self.cofs[i] = float(self.zBoxes[i].text())
            except:
                self.cofs[i] = 0;        
        self.z = zernike.Zernike(self.cofs, self.wavelength,self.xRes)
        
        idealLayers = self.minheight(np.nanmax(self.z.opd)-np.nanmin(self.z.opd),self.voxelsize,self.n1,self.n2)
        self.idealLayersLabel.setText("Optimum Layer #: "+str(int(np.ceil(idealLayers))))
        self.isCalculated = True
        self.plot_zernike(self.dimensionMenu.currentIndex())
        
        print idealLayers
        self.toggle(1)
        
   
        
        
    def make_widgets(self):
        
        #Menu Action Declarations
        
        self.saveSystemAction = QAction('&Save System',self)
        self.saveSystemAction.setShortcut('Ctrl+S')
        self.saveSystemAction.triggered.connect(self.save_system)
        
        self.loadSystemAction = QAction('&Load System',self)
        self.loadSystemAction.setShortcut('Ctrl+L')
        self.loadSystemAction.triggered.connect(self.load_system)
        
        
        #Menus and StatusBar Declarations / Implementations
        self.menus = self.menuBar()
        self.fileMenu = self.menus.addMenu('&File')
        self.fileMenu.addAction(self.saveSystemAction)
        self.fileMenu.addAction(self.loadSystemAction)
        #self.fileMenu.addAction(self.exitAction)
        
        #self.settingsMenu = self.menus.addMenu('&Settings')
        #self.settingsMenu.addAction(self.changeWavelengthAction)
        #self.settingsMenu.addAction(self.changeHeightAction)
        
        self.status = self.statusBar()
        self.status.showMessage(self.readyMessage)
        
        
        #Text boxes and labels for zernike coefficients    
        for i in range(37):
            self.zBoxes[i] = QLineEdit('0')
            self.zBoxLabels[i] = QLabel("Z%s" %str(i+1))
        
        #Dither 3D Stuff
        self.layersMenu = QComboBox()
        self.layersMenu.addItem("Composite")
        self.layersMenu.setEnabled(False)
        self.layersMenu.currentIndexChanged[str].connect(self.show_dithered_layer)

        #Zernike 3D Stuff
        
        self.enable3dButton = QCheckBox("Enable 3D (slow)")
        self.enable3dButton.setToolTip("Recommended for a diameter of less than 10mm")
        self.enable3dButton.stateChanged[int].connect(self.enable_3d)
        
        self.dimensionMenu = QComboBox()
        self.dimensionMenu.addItem("2d")
        self.dimensionMenu.addItem("3d")
        self.dimensionMenu.setEnabled(False)
        self.dimensionMenu.currentIndexChanged[int].connect(self.plot_zernike)
        
        #Other Labels
        self.idealLayersLabel = QLabel("Minimum Ideal Layers")
        self.ditherLabel = QLabel("Dithered Image")
        
        self.wavelengthBox = QLineEdit(str(self.wavelength))
        self.wavelengthLabel = QLabel("Wavelength (nm):")
        
        self.heightBox = QLineEdit(str(self.height))
        self.heightLabel = QLabel("Layers:")
        
        self.n1Box = QLineEdit(str(self.n1))
        self.n1Label = QLabel("n1:")
        
        self.n2Box = QLineEdit(str(self.n2))
        self.n2Label = QLabel("n2")
        
        self.voxelSizeBox = QLineEdit(str(self.voxelSize))
        self.voxelSizeLabel = QLabel(u"Voxel Size (μm):")
        
        self.diameterBox = QLineEdit(str(self.diameter))
        self.diameterLabel = QLabel("Diameter (mm):")
        
        #Button Declarations and Statuses
        """self.buttonSlice = QPushButton('Slice')
        self.buttonSlice.setEnabled(False)
        self.buttonSlice.pressed.connect(self.slice_layers)"""       
        
        self.buttonExportAll = QPushButton('Export All')
        self.buttonExportAll.setEnabled(False)
        self.buttonExportAll.pressed.connect(self.export_all)
        
        
        
        self.buttonDither = QPushButton('Dither')
        self.buttonDither.setEnabled(False)
        self.buttonDither.pressed.connect(self.dither)
        
        
        self.buttonLayersFast = QPushButton('Generate Layers (Packed)')
        self.buttonLayersFast.setEnabled(False)
        self.buttonLayersFast.pressed.connect(self.generate_layers_fast)        
        
        
        self.buttonLayers = QPushButton('Generate Layers (sparse)')
        self.buttonLayers.setEnabled(False)
        self.buttonLayers.pressed.connect(self.generate_layers)
        
        self.buttonCalc = QPushButton('Calculate')
        self.buttonCalc.pressed.connect(self.calculate)
        
        """self.button2D = QPushButton('2D')
        self.button2D.setCheckable(True)
        self.button2D.setChecked(True)
        
        self.button3D = QPushButton('3D')
        self.button3D.setCheckable(True)
        self.button3D.setEnabled(False)
        
        self.zButtons = QButtonGroup()
        self.zButtons.addButton(self.button2D)
        self.zButtons.addButton(self.button3D)"""
        
        #Graphs
        self.zernikeFig = plt.figure()
        self.zernikeCanvas = FigureCanvas(self.zernikeFig)
        plt.close(self.zernikeFig)
        
        self.ditherFig = plt.figure()
        self.ditherCanvas = FigureCanvas(self.ditherFig)
        plt.close(self.ditherFig)
        
        self.sliceFig = plt.figure()
        self.sliceCanvas = FigureCanvas(self.sliceFig)
        plt.close(self.sliceFig)
        
        
        
        
        
        
        
    def set_geometry(self):
        
        #Section Declarations
        self.top = QWidget()
        self.left = QWidget()
        self.middle = QWidget()
        self.right = QWidget()
        self.bottom = QWidget()
        self.mainWindow = QWidget()
        
        
        #Layout Declarations
        self.topLayout = QGridLayout()
        self.leftLayout = QGridLayout()        
        self.middleLayout = QGridLayout()        
        self.rightLayout = QGridLayout()
        self.bottomLayout = QGridLayout()        
        self.mainLayout = QGridLayout()  
        
        
        #Adding Widgets to Layouts
        
        self.topLayout.addWidget(self.wavelengthLabel,0,0)
        self.topLayout.addWidget(self.wavelengthBox,0,1)
        
        self.topLayout.addWidget(self.heightLabel,1,0)
        self.topLayout.addWidget(self.heightBox,1,1)        
        
        self.topLayout.addWidget(self.n1Label,0,2)
        self.topLayout.addWidget(self.n1Box,0,3)
        
        self.topLayout.addWidget(self.n2Label,1,2)
        self.topLayout.addWidget(self.n2Box,1,3)
        
        self.topLayout.addWidget(self.voxelSizeLabel,0,4)
        self.topLayout.addWidget(self.voxelSizeBox,0,5)
        
        self.topLayout.addWidget(self.diameterLabel,1,4)
        self.topLayout.addWidget(self.diameterBox,1,5)
        
        self.topLayout.addWidget(self.idealLayersLabel,2,0,1,2)
                
        
        #self.leftLayout.addWidget(self.buttonSlice,1,2,1,2)
        self.leftLayout.addWidget(self.buttonDither,1,2,1,2)
        self.leftLayout.addWidget(self.buttonCalc,1,0,1,2)
        for i in range(37):
            if i <= 18:
                self.leftLayout.addWidget(self.zBoxLabels[i],i+2,0)
                self.leftLayout.addWidget(self.zBoxes[i],i+2,1)
            else:
                self.leftLayout.addWidget(self.zBoxLabels[i],i-17,2)
                self.leftLayout.addWidget(self.zBoxes[i],i-17,3)
        
        self.middleLayout.addWidget(self.enable3dButton)
        self.middleLayout.addWidget(self.dimensionMenu,0,1)
        #self.middleLayout.addWidget(self.button3D,0,1)
        self.middleLayout.addWidget(self.zernikeCanvas,1,0,5,5)
        
        self.rightLayout.addWidget(self.buttonLayers,0,0)
        self.rightLayout.addWidget(self.buttonLayersFast,0,1)
        self.rightLayout.addWidget(self.layersMenu,0,2)
        self.rightLayout.addWidget(self.ditherCanvas,1,0,5,5)
        self.rightLayout.addWidget(self.buttonExportAll,0,3)
        
        self.bottomLayout.addWidget(self.sliceCanvas,0,0)
        
        #Layout Implementations
        self.top.setLayout(self.topLayout)
        self.left.setLayout(self.leftLayout)
        self.middle.setLayout(self.middleLayout)
        self.right.setLayout(self.rightLayout)
        self.bottom.setLayout(self.bottomLayout)
        
        
    
        #Main Window Stuff
        self.mainLayout.addWidget(self.top,0,0,1,2)
        self.mainLayout.addWidget(self.left,1,0)
        self.mainLayout.addWidget(self.middle,1,1)
        self.mainLayout.addWidget(self.right,1,2)
        #self.mainLayout.addWidget(self.bottom,2,2)
        self.mainWindow.setLayout(self.mainLayout)
        self.setCentralWidget(self.mainWindow)
        self.setGeometry(100, 100, 1500, 600)
        self.setWindowTitle('Grin Phase Plate Design')
        self.show()
            
        
    
app = QApplication(sys.argv)
form = Zplot()
app.exec_()   