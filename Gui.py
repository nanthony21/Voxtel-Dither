# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 37:49:19 2014

@author: Nathan Neibauer and Nick Anthony
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
        self.wavelength = 0.633
        self.voxelSize = 20
        self.xRes = 1000
        self.height = 13
        self.printedradius=20
        self.normratio=1
        self.n1 = 1.5064
        self.n2 = 1.5438
        self.is3d = False
        self.isCalculated = False
        self.isDithered = False
        self.isLayersGenerated = False
        self.invert=False
        self.fringe=False
        self.params = [self.cofs, self.wavelength, self.height, self.n1, 
                       self.n2, self.xRes, self.voxelSize,self.normratio,self.invert,self.fringe]
        
        sys.stdout=EmittingStream()        
        self.connect(sys.stdout,SIGNAL('textWritten(QString)'),self.normalOutputWritten)        
        sys.stderr=EmittingStream()        
        self.connect(sys.stderr,SIGNAL('textWritten(QString)'),self.normalOutputWritten)        
                #Status Bar Messages        
        self.calcMessage = "CALCULATING..."
        self.readyMessage = "READY"
        
        
        
        
        
        #Draw Everything
        self.make_widgets()  
        self.set_geometry()
        
        self.update_boxes()
        
        
        
        
        
     
    def save_system(self):
        #Saves all variables to a text file.
        self.update_values()
        dialog = QFileDialog(self)
        filename = dialog.getSaveFileName()
        f = open(filename,'w')
        
        for i in self.params:
            f.write(str(i)+"\n")
            
        
        
                      
        
    def load_system(self):
        #Loads all variables from a text file
        dialog = QFileDialog(self)
        filename = dialog.getOpenFileName()
        f = open(filename,'r')
        
        def str_bool_check(a):
            if a=="False\n":
                return False
            elif a=='True\n':
                return True
            else:
                raise ValueError('Expected boolean value but instead got %s.'%a)
        
        self.cofs = ast.literal_eval(f.readline())
        self.wavelength = float(f.readline())
        self.height = int(f.readline())
        self.n1 = float(f.readline())
        self.n2 = float(f.readline())
        self.xRes = int(f.readline())
        self.voxelSize = float(f.readline())
        self.normratio=float(f.readline())
        self.invert=str_bool_check(f.readline())        
        self.fringe=str_bool_check(f.readline())        
        print self.invert, self.fringe
        self.update_boxes()
    
    def helptext(self):
        ''' adds the "help" option to the file menu and brings up the help dialog'''
        dialog=QDialog(parent=None)
        dialog.main=QWidget()
        dialog.layout=QGridLayout()
        dialog.text=QLabel('This Software is capable of generating print '
        'files for Standard Zernike coefficients 1-37 and Fringe Zernike '
        'coefficients 1-25 as defined by Zemax OpticStudio.\n'
        'Note that for both cases the first coefficient, Piston, '
        'will not affect your results.\n\n'
        'The coefficients are in units of waves, so if you have a set of coefficients '
        'that are in units of mm you can achieve this by setting the "Wavelength" box to 1 mm.\n\n'
        'Normalization Ratio: This is the ratio of the diameter of the unit circle to which the Zernike '
        'coefficients are normalized to the diameter of the printed optic.\nThis value must be '
        'greater than one because the Zernike polynomials are not orthogonal outside of the unit circle.\n\n'
        'Invert: By default the input zernike coefficients describe the optical path '
        'length (OPL) of the printed optic. If you want the Zernike coefficients to describe\nthe wavefront'
        'that is created by the optic then you will want to have the "Invert" box checked, this will effectively multiply all of the coefficients by -1.')
        dialog.layout.addWidget(dialog.text,1,1)
        dialog.setLayout(dialog.layout)
        dialog.setWindowTitle('Help')
        dialog.show()
        dialog.exec_()
        

    
    def update_values(self):
        """Queries all text boxes and updates the system parameters"""
        self.voxelSize = float(self.voxelSizeBox.text())
        self.xRes = int(self.resolutionBox.text())
        self.wavelength = float(self.wavelengthBox.text())
        self.normratio=float(self.normratioBox.text())
        self.n1 = float(self.n1Box.text())
        self.n2 = float(self.n2Box.text())
        
        if self.invertCheck.isChecked():
            self.invert=True
        else:
            self.invert=False
        for i in range(37):
            try:
                self.cofs[i] = float(self.zBoxes[i].text())
            except:
                self.cofs[i] = 0; 
       
        self.params = [self.cofs, self.wavelength, self.height, self.n1, 
                       self.n2, self.xRes, self.voxelSize,self.normratio,self.invert,self.fringe]
        
        
    def update_boxes(self):
        for i in range(len(self.zBoxes)):
            self.zBoxes[i].setText(str(self.cofs[i]))
        self.wavelengthBox.setText(str(self.wavelength))
        self.heightBox.setText(str(self.height))
        self.n1Box.setText(str(self.n1))
        self.n2Box.setText(str(self.n2))
        self.normratioBox.setText(str(self.normratio))
        self.voxelSizeBox.setText(str(self.voxelSize))
        self.resolutionBox.setText(str(self.xRes))
        self.invertCheck.setChecked(self.invert)
        
        if self.fringe==True:
            self.fringeButton.setText('Use Standard Coefficients')
            for i in range(25,37):
                self.zBoxes[i].setEnabled(False)
                self.zBoxes[i].setText('0')
                self.cofs[i]=0
        else:
            self.fringeButton.setText('Use Fringe Coefficients')
            for i in range(25,37):
                self.zBoxes[i].setEnabled(True)
                
    
    def toggle(self, state):
        """Assigns bool 'state' to the value of the setEnabled() parameter 
        for each button"""
        self.buttonCalc.setEnabled(state)
        self.buttonDither.setEnabled(state)
        #self.button2D.setEnabled(state)
        if self.isDithered:
            self.buttonLayers.setEnabled(state)
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
        self.zernikeAxes1.set_title('Max OPD: %.4f ($\mu$m)'%(np.nanmax(self.z.opd)-np.nanmin(self.z.opd)))
        self.zernikeAxes2.set_xlabel("OPD ($\mu$m)")
        
    def __generate_dither_axes(self):
        """Same as above function, but for the dither plot"""
        self.ditherFig.clear()
        self.ditherGridSpec = gs.GridSpec(1,2,width_ratios = [1,.1])
        self.ditherAxes1 = self.ditherFig.add_subplot(self.ditherGridSpec[0])
        self.ditherAxes2 = self.ditherFig.add_subplot(self.ditherGridSpec[1])
        self.ditherAxes2.set_xlabel('Effective n')
        
        
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
       
        self.toggle(0)
        try:        
            self.d.calc3d()
            self.layersMenu.setEnabled(True)
            self.layersMenu.clear()
            self.layersMenu.addItem("Composite")
            for i in range(self.height):
                self.layersMenu.addItem("Layer %g" %(i+1))
                
            self.isLayersGenerated = True
            
        except Exception, err:
            sys.stderr.write('Error: %s\n' % str(err))
        self.toggle(1)
        
        
        
        
        
    def show_dithered_layer(self, layer):
        try:
            if layer[-1] == 'e':
                self.__generate_dither_axes()
                image = self.ditherAxes1.imshow(np.float32(self.d.dithered),interpolation='none',cmap='jet')
                self.ditherFig.colorbar(image, self.ditherAxes2)
                self.ditherCanvas.draw()
            else:
                if layer[-2] == ' ':
                    index = int(layer[-1])-1
                else:
                    index = 10*int(layer[-2])+int(layer[-1])-1
                
                self.ditherFig.clear()
                self.ditherAxes = self.ditherFig.add_subplot(111)
                self.ditherAxes.imshow(self.d.data[:,:,index],interpolation='none',cmap='jet')
                self.ditherAxes.set_title('Red=N1 Blue=N2')                
                self.ditherCanvas.draw()
            
        except:
            return 1
    
        
    def dither(self):
        self.toggle(0)
        self.update_values()
        try:        
            self.isDithered = True
            self.status.showMessage(self.calcMessage)
            
    
            self.d = dither.Dither(self.z.opd, self.height,self.voxelSize,self.n1,self.n2)                 
            
            self.status.showMessage(self.readyMessage)
            
            self.__generate_dither_axes()
            image = self.ditherAxes1.imshow(np.float32(self.d.dithered),interpolation='none',cmap='jet')
            self.ditherFig.colorbar(image, self.ditherAxes2)
            self.ditherCanvas.draw()
        except Exception,err:
            sys.stderr.write('Error: %s\n'%str(err))
        self.toggle(True)
 
   
    
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
    
        
     
    def plot_zernike(self):
        if self.isCalculated:
            self.__generate_zernike_axes()
            image = self.zernikeAxes1.imshow(self.z.opd,cmap='jet',interpolation='none')
            
            self.zernikeFig.colorbar(image, self.zernikeAxes2)        
            self.zernikeFig.hold(False)
            self.zernikeCanvas.draw()

        
        
        
        
    def calculate(self):
        self.toggle(0)
        self.update_values()
        global zopd
        
       
        try:
            self.z = zernike.Zernike(self.cofs, self.wavelength,self.xRes,self.normratio,correcting=self.invert,fringe=self.fringe)
            
            zopd=self.z.opd
            self.idealLayers = int(np.ceil(self.minheight(np.nanmax(self.z.opd)-np.nanmin(self.z.opd),self.voxelSize,self.n1,self.n2)))
            self.idealLayersLabel.setText("Optimum Layer #: "+str(self.idealLayers))
            self.isCalculated = True
            self.isDithered=False
            self.plot_zernike()       
            self.height = self.idealLayers
            self.heightBox.setText(str(self.height))
            print "Done Calculating"
        except Exception, err:
            sys.stderr.write('Error: %s\n' % str(err))

            
        self.toggle(1)
        
       
    def fringeswitch(self):
        if self.fringe==False:
            self.fringe=True
            self.fringeButton.setText('Use Standard Coefficients')
            for i in range(25,37):
                self.zBoxes[i].setEnabled(False)
                self.zBoxes[i].setText('0')
                self.cofs[i]=0
        else:
            self.fringe=False
            self.fringeButton.setText('Use Fringe Coefficients')
            for i in range(25,37):
                self.zBoxes[i].setEnabled(True)
        
        
    def make_widgets(self):
        
        #Menu Action Declarations
        
        self.saveSystemAction = QAction('&Save System',self)
        self.saveSystemAction.setShortcut('Ctrl+S')
        self.saveSystemAction.triggered.connect(self.save_system)
        
        self.loadSystemAction = QAction('&Load System',self)
        self.loadSystemAction.setShortcut('Ctrl+L')
        self.loadSystemAction.triggered.connect(self.load_system)
        
        self.help=QAction('Help',self)
        self.help.triggered.connect(self.helptext)
        
        #Menus and StatusBar Declarations / Implementations
        self.menus = self.menuBar()
        self.fileMenu = self.menus.addMenu('&File')
        self.fileMenu.addAction(self.saveSystemAction)
        self.fileMenu.addAction(self.loadSystemAction)
        self.fileMenu.addAction(self.help)        
        
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

        
        #Other Labels
        self.idealLayersLabel = QLabel("Optimum Layer #:")
        self.ditherLabel = QLabel("Dithered Image")
        
        self.wavelengthBox = QLineEdit(str(self.wavelength))
        self.wavelengthLabel = QLabel(u"Wavelength (μm):")
        
        self.heightBox = QLineEdit(str(self.height))
        self.heightLabel = QLabel("Layers:")
        
        self.n1Box = QLineEdit(str(self.n1))
        self.n1Label = QLabel("n1:")
        
        self.n2Box = QLineEdit(str(self.n2))
        self.n2Label = QLabel("n2:")
        
        self.voxelSizeBox = QLineEdit(str(self.voxelSize))
        self.voxelSizeLabel = QLabel(u"Voxel Height (μm):")
        
        self.resolutionBox = QLineEdit(str(self.xRes))
        self.resolutionLabel = QLabel("Diameter (pixels):")
        
        self.normratioBox=QLineEdit(str(self.normratio))
        self.normratioLabel=QLabel('Normalization Ratio:')
        
        #Button Declarations and Statuses
       
        self.invertCheck=QCheckBox()
        self.invertLabel=QLabel('Invert:')
        
        self.fringeButton=QPushButton('Use Fringe Coefficients')
        self.fringeButton.setEnabled(True)
        self.fringeButton.pressed.connect(self.fringeswitch)
        
        
        self.buttonExportAll = QPushButton('Export All')
        self.buttonExportAll.setEnabled(False)
        self.buttonExportAll.pressed.connect(self.export_all)
        
        
        
        self.buttonDither = QPushButton('Dither')
        self.buttonDither.setEnabled(False)
        self.buttonDither.pressed.connect(self.dither)
        
                
        
        
        self.buttonLayers = QPushButton('Generate Layers')
        self.buttonLayers.setEnabled(False)
        self.buttonLayers.pressed.connect(self.generate_layers)
        
        self.buttonCalc = QPushButton('Calculate')
        self.buttonCalc.pressed.connect(self.calculate)
        
        self.text=QTextEdit()
        self.text.setReadOnly(True)
        
        
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
        
        self.topLayout.addWidget(self.wavelengthLabel,0,0,1,1)
        self.topLayout.addWidget(self.wavelengthBox,0,1)
        
        self.topLayout.addWidget(self.heightLabel,1,0,alignment=Qt.AlignRight)
        self.topLayout.addWidget(self.heightBox,1,1)        
        
        self.topLayout.addWidget(self.n1Label,0,2)
        self.topLayout.addWidget(self.n1Box,0,3)
        
        self.topLayout.addWidget(self.n2Label,1,2)
        self.topLayout.addWidget(self.n2Box,1,3)
        
        self.topLayout.addWidget(self.voxelSizeLabel,0,4,alignment=Qt.AlignRight)
        self.topLayout.addWidget(self.voxelSizeBox,0,5)
        
        self.topLayout.addWidget(self.resolutionLabel,1,4)
        self.topLayout.addWidget(self.resolutionBox,1,5)
        
        self.topLayout.addWidget(self.normratioLabel,0,6)
        self.topLayout.addWidget(self.normratioBox,0,7)
        
        self.topLayout.addWidget(self.idealLayersLabel,2,0,1,2)
        
        self.topLayout.addWidget(self.invertLabel,1,6)
        self.topLayout.addWidget(self.invertCheck,1,6,1,1,alignment=Qt.AlignCenter)
        
        self.topLayout.addWidget(self.fringeButton,1,7)
        
        self.topLayout.addWidget(self.text,0,12,2,1)
                
        
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
        
        self.middleLayout.addWidget(self.zernikeCanvas,1,0,5,5)
        
        self.rightLayout.addWidget(self.buttonLayers,0,0)
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
          
    def normalOutputWritten(self, text):
        """Append text to the QTextEdit."""
        # Maybe QTextEdit.append() works as well, but this is how I do it:
        cursor = self.text.textCursor()
        cursor.movePosition(QTextCursor.End)
        cursor.insertText(text)
        self.text.setTextCursor(cursor)
        self.text.ensureCursorVisible()

class EmittingStream(QObject):

    textWritten = pyqtSignal(str)

    def write(self, text):
        self.textWritten.emit(str(text))
       

app = QApplication(sys.argv)
form = Zplot()
app.exec_()   