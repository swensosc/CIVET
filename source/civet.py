#! /usr/bin/env python
import sys
import argparse
import string
import subprocess
import numpy as np
import netCDF4 as netcdf4
import matplotlib
import matplotlib.pyplot as plt
import libcivet

# append path before importing pyqt libraries
libdir = './pyqtgraph/pyqtgraph-develop'
libdir = '/project/tss/swensosc/pylibs/pyqtgraph/pyqtgraph-develop'
sys.path.append(libdir)

# import pyqtgraph libraries  --------------------------------------
import pyqtgraph as pg
import pyqtgraph.exporters
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
from pyqtgraph.graphicsItems.ROI import *

#----------  Definitions  -----------------------
#class CivetGUI(QtGui.QWidget):
class CivetGUI(pg.GraphicsWindow):
    
    def __init__(self):
        #QtGui.QWidget.__init__(self)
        pg.GraphicsWindow.__init__(self)

        xs=1000
        ys=800
        #self.resize(xs,ys)

        self.setupGUI()
        self.createLUT()

        self.addAttributes()

        self.addParameterTree()

        ctitle = 'CIVET: The CTSM Interactive Visualization and Exploration Tool'
        self.setWindowTitle(ctitle)
        #self.setBackground((255,255,255))
        #self.setBackground((249,243,222))
        #self.setBackground((227,238,247))
        self.setBackground((127,165,122))
        
    def Civet(self,infile):
        self.readData(infile,decimalTime=self.decimalTime)
        self.create_colorbar(cbx=10,cby=256)
        self.plotMap()

        if self.param_tree.param('Map Annotation', 'Coastlines').value():
            self.plot_coastlines()
        self.addGlobalRange()
        self.plotTimeSeries()

    def addAttributes(self):
        self.tsGlobalYMin = np.float(0)
        self.tsGlobalYMax = np.float(1)
        self.tsGlobalXMin = np.float(0)
        self.tsGlobalXMax = np.float(1)

        self.tsLocalYMin = self.tsGlobalYMin
        self.tsLocalYMax = self.tsGlobalYMax

        self.tsCurrentYMin = self.tsGlobalYMin
        self.tsCurrentYMax = self.tsGlobalYMax
        self.tsCurrentXMin = self.tsGlobalXMin
        self.tsCurrentXMax = self.tsGlobalXMax

        self.decimalTime = True
        
        self.tsMean  = np.float(-999)
        self.tsTrend = np.float(-999)
        self.tsAmp   = np.float(-999)
        
        self.lonMin = np.float(0)
        self.lonMax = np.float(360)
        self.latMin = np.float(-90)
        self.latMax = np.float(90)

        self.westLon  = self.lonMin
        self.eastLon  = self.lonMax
        self.southLat = self.latMin
        self.northLat = self.latMax

        self.mapTime = 0
        self.xTimeSeries = 20
        self.yTimeSeries = 90
        
        self.mapCurrentMin  = 0
        self.mapCurrentMax  = 1
        self.mapGlobalMin  = 0
        self.mapGlobalMax  = 1

        self.gaussWidth = 0
        self.maxGaussWidth = 10
        self.bilinResMultiplier = 1
        self.maxBRMultiplier = 4
        self.logTransform = False

        self.doSelectPoint = False
        self.doSelectRegion = False

        self.doAutoScaleUpdateXY = False
        self.doTSDecomposition = False

        self.tsPenWidth = 2
        self.coastPenWidth = 1.5
        self.riverPenWidth = 1.5
        self.countryPenWidth = 1.0
        self.statePenWidth = 1.0

        self.tsPenColor = (0,0,0)
        self.coastPenColor = (0,0,0)
        self.riverPenColor = (0,0,0)
        self.countryPenColor = (0,0,0)
        self.statePenColor = (255,0,0)

    def WhatIsACivet(self):
        # a window created in the function closes after function call
        #self.wiac_window = pg.GraphicsWindow()
        xs,ys = (600,400)
        self.wiac_window.resize(xs,ys)
        self.wiac_window.setWindowTitle('This is a Civet')
        self.wiac_window.setBackground((255,255,255))
        self.wiac_window.setVisible(True)

        logo_gv = pg.GraphicsView(parent=self.wiac_window)
        logo_gv.resize(0.95*xs,0.95*ys)

        vb_logo  = pg.ViewBox()
        logo_gv.addItem(vb_logo)
        img_logo = pg.ImageItem()
        vb_logo.addItem(img_logo)
        vb_logo.resize(0.95*xs,0.95*ys)
        vb_logo.setBackgroundColor('w')
        self.addLogo(img_logo)

        #self.wiac_window.show()
        logo_gv.show()
        
    def setupGUI(self):

        # create main layout  ---------------------------------
        self.main_layoutwidget=pg.LayoutWidget()
        self.setLayout(self.main_layoutwidget.layout)
        xs=self.width()
        ys=self.height()
        self.main_layoutwidget.resize(xs,ys)

        # create sub-layouts  ---------------------------------
        self.data_layoutwidget = pg.LayoutWidget()
        self.main_layoutwidget.addWidget(self.data_layoutwidget,col=0,row=0)
        self.param_layoutwidget = pg.LayoutWidget()
        self.main_layoutwidget.addWidget(self.param_layoutwidget,col=1,row=0)

        # set main layout attributes
        self.dfrac = 0.67
        data_xs = self.dfrac * xs
        param_xs = (1.0 - self.dfrac) * xs
        self.main_layoutwidget.layout.setColumnMinimumWidth(0,data_xs)
        self.main_layoutwidget.layout.setColumnMinimumWidth(1,param_xs)
        
        # set parameter layout attributes
        self.ptfrac = 0.95
        param_top_height    = ys*self.ptfrac
        param_bottom_height = ys*(1.0-self.ptfrac)
        self.param_layoutwidget.layout.setRowMinimumHeight(0, param_top_height)

        self.param_layoutwidget.nextRow()
        self.ptwidget = ParameterTree()
        self.param_layoutwidget.addWidget(self.ptwidget,  col=0, row=0, colspan=1, rowspan=1)

        wiac_button = QtGui.QPushButton("What is a Civet?")
        wiac_button.clicked.connect(self.WhatIsACivet)
        self.wiac_window = pg.GraphicsWindow()
        self.wiac_window.setVisible(False)
        self.param_layoutwidget.addWidget(wiac_button,  col=0, row=1, colspan=1, rowspan=1)
        # add glw to data_layoutwidget to hold viewboxes
        self.data_glw = pg.GraphicsLayoutWidget()
        self.data_layoutwidget.addWidget(self.data_glw,col=0)
        self.data_glw.setBackground((255,255,255))

        # from QGraphicsGridLayout
        qldata=self.data_layoutwidget.layout
        qlparam=self.param_layoutwidget.layout
    
        dcolminwidth=np.array((0.95,)) * data_xs
        for i in range(dcolminwidth.size):
            qldata.setColumnMinimumWidth(i, dcolminwidth[i])
            
        pcolminwidth=np.array((0.95,)) * param_xs
        for i in range(pcolminwidth.size):
            qlparam.setColumnMinimumWidth(i, pcolminwidth[i])

        # from QGraphicsGridLayout
        self.data_layout_widget = self.data_glw.addLayout(colspan=1)

        # top layout has 4 columns (2 vb, 2 axes)
        # and 3 rows (label, images, bottom axis)
        self.top_layout = self.data_layout_widget.addLayout(colspan=3)
        self.qltop=self.top_layout.layout
        top_ys = ys*0.6
        self.qltop.setRowMinimumHeight(1, top_ys)
        self.qltop.setRowMaximumHeight(1, top_ys)

        self.tcol_width = np.array((0.075,0.05,0.8,0.075))
        tcolminwidth = self.tcol_width * data_xs * 0.9
        tcolmaxwidth = self.tcol_width * data_xs * 0.9
        for i in range(tcolminwidth.size):
            self.qltop.setColumnMinimumWidth(i, tcolminwidth[i])
            self.qltop.setColumnMaximumWidth(i, tcolmaxwidth[i])

        # create labels in top layout
        self.top_data_label_left = pg.LabelItem(justify='left')
        self.top_layout.addItem(self.top_data_label_left,row=0,col=0,colspan=2)
        self.top_data_label_right = pg.LabelItem(justify='right')
        self.top_layout.addItem(self.top_data_label_right,row=0,col=2)
        self.top_data_label_left.setText('')
        self.top_data_label_right.setText('')
        
        # create viewboxes in top layout
        self.vb_cbar  = self.top_layout.addViewBox(row=1, col=1)
        self.vb_image = self.top_layout.addViewBox(row=1, col=2)
        self.map_image = pg.ImageItem()
        self.vb_image.addItem(self.map_image)
    
        self.colorbar = pg.ImageItem()
        self.vb_cbar.addItem(self.colorbar)

        # add axes to top layout
        self.colorbar_left_axis=pg.AxisItem('left')
        self.colorbar_left_axis.linkToView(self.vb_cbar)
        self.top_layout.addItem(self.colorbar_left_axis,row=1, col=0)

        map_right_axis=pg.AxisItem('right')
        map_right_axis.linkToView(self.vb_image)
        self.top_layout.addItem(map_right_axis,row=1, col=3)

        map_bottom_axis=pg.AxisItem('bottom')
        map_bottom_axis.linkToView(self.vb_image)
        self.top_layout.addItem(map_bottom_axis,row=2, col=2)

        self.map_overplot_rivers = pg.PlotDataItem()
        self.vb_image.addItem(self.map_overplot_rivers)
        self.map_overplot_rivers.setZValue(10)

        self.map_overplot_countries = pg.PlotDataItem()
        self.vb_image.addItem(self.map_overplot_countries)
        self.map_overplot_countries.setZValue(11)

        self.map_overplot_states = pg.PlotDataItem()
        self.vb_image.addItem(self.map_overplot_states)
        self.map_overplot_states.setZValue(12)

        self.map_overplot_coasts = pg.PlotDataItem()
        self.vb_image.addItem(self.map_overplot_coasts)
        self.map_overplot_coasts.setZValue(13)

       # shift down 1 row and add bottom layout
        self.data_layout_widget.nextRow()
        # bottom layout has 2 columns (1 plot, 1 axis)
        # and 4 rows (rangeplot, label, plot, axis)
        self.bottom_layout = self.data_layout_widget.addLayout(colspan=3)
        self.qlbottom=self.bottom_layout.layout
        bot_ys = (ys - top_ys)
        range_time_ys = 0.15 * bot_ys
        plot_ys  = bot_ys - range_time_ys - 0.05
        self.qlbottom.setRowMinimumHeight(1, range_time_ys)
        self.qlbottom.setRowMaximumHeight(1, range_time_ys)
        self.qlbottom.setRowMinimumHeight(2, plot_ys)
        self.qlbottom.setRowMaximumHeight(2, plot_ys)

        self.bcol_width = np.array((0.075,0.35,0.425,0.1))
        bcolminwidth = self.bcol_width * data_xs * 0.9
        bcolmaxwidth = self.bcol_width * data_xs * 0.9
        for i in range(bcolminwidth.size):
            self.qlbottom.setColumnMinimumWidth(i, bcolminwidth[i])
            self.qlbottom.setColumnMaximumWidth(i, bcolmaxwidth[i]) 

        # create labels for bottom plot
        self.bottom_data_label_left = pg.LabelItem(justify='left')
        self.bottom_data_label_left.setText('')
        self.bottom_layout.addItem(self.bottom_data_label_left,row=0,col=1)
        self.bottom_data_label_right = pg.LabelItem(justify='right')
        self.bottom_data_label_right.setText('')
        self.bottom_layout.addItem(self.bottom_data_label_right,row=0,col=2)

        # create plots to contain linearregionitem
        self.range_time_plot = pg.PlotItem()
        self.bottom_layout.addItem(self.range_time_plot,row=1,col=1,colspan=2)
        self.range_data_plot = pg.PlotItem()
        self.bottom_layout.addItem(self.range_data_plot,row=2,col=3)

        # create main plot in bottom layout
        self.time_series_plot = pg.PlotItem()
        self.time_series_plot.autoBtn.clicked.connect(self.autoBtnClicked)
        self.bottom_layout.addItem(self.time_series_plot,row=2,col=1,colspan=2)

        self.ts_bottom_axis=pg.AxisItem('bottom')
        self.ts_bottom_axis.linkToView(self.time_series_plot.vb)
        self.bottom_layout.addItem(self.ts_bottom_axis,row=3, col=1,colspan=2)
        
        ts_top_left_axis=pg.AxisItem('left')
        ts_top_left_axis.linkToView(self.time_series_plot.vb)
        self.bottom_layout.addItem(ts_top_left_axis,row=1, col=0)
        self.ts_bot_left_axis=pg.AxisItem('left')
        self.ts_bot_left_axis.linkToView(self.time_series_plot.vb)
        self.bottom_layout.addItem(self.ts_bot_left_axis,row=2, col=0)
        ts_top_left_axis.setVisible(False)
        #self.time_series_plot.getAxis('right').setStyle(showValues=False)
        
    def resizeGUI(self):
        xs=self.width()
        ys=self.height()
        self.main_layoutwidget.resize(xs,ys)
        # set main layout attributes
        data_xs = self.dfrac * xs
        param_xs = (1.0 - self.dfrac) * xs
        self.main_layoutwidget.layout.setColumnMinimumWidth(0,data_xs)
        self.main_layoutwidget.layout.setColumnMinimumWidth(1,param_xs)
        
        # set parameter layout attributes
        param_top_height    = ys*self.ptfrac
        param_bottom_height = ys*(1.0-self.ptfrac)
        self.param_layoutwidget.layout.setRowMinimumHeight(0, param_top_height)
        self.param_layoutwidget.layout.setRowMinimumHeight(1, param_bottom_height)

        #self.vb_logo.resize(0.95*param_xs,param_bottom_height)

        # from QGraphicsGridLayout
        qldata=self.data_layoutwidget.layout
        qlparam=self.param_layoutwidget.layout
    
        dcolminwidth=np.array((0.925,)) * data_xs
        for i in range(dcolminwidth.size):
            qldata.setColumnMinimumWidth(i, dcolminwidth[i])
            
        pcolminwidth=np.array((0.925,)) * param_xs
        for i in range(pcolminwidth.size):
            qlparam.setColumnMinimumWidth(i, pcolminwidth[i])

        # from QGraphicsGridLayout
        self.data_layout_widget.layout.setRowMinimumHeight(0, ys*0.5)
        self.data_layout_widget.layout.setRowMaximumHeight(0, ys*0.5)
        
        tcolminwidth = self.tcol_width * data_xs * 0.9
        tcolmaxwidth = self.tcol_width * data_xs * 0.9
        for i in range(tcolminwidth.size):
            self.qltop.setColumnMinimumWidth(i, tcolminwidth[i])
            self.qltop.setColumnMaximumWidth(i, tcolmaxwidth[i])

        bcolminwidth = self.bcol_width * data_xs * 0.9
        bcolmaxwidth = self.bcol_width * data_xs * 0.9
        for i in range(bcolminwidth.size):
            self.qlbottom.setColumnMinimumWidth(i, bcolminwidth[i])
            self.qlbottom.setColumnMaximumWidth(i, bcolmaxwidth[i]) 

    def addLogo(self,imageItem):
        from PIL import Image
        import os
        dir_path = os.path.dirname(os.path.realpath(__file__))
        imgfile = dir_path +'/../images/civet_woodcut.jpg'
        jpgfile = Image.open(imgfile)
        logo = np.array(jpgfile.getdata()).reshape(jpgfile.size[1],
                                                   jpgfile.size[0], 3)
        logo = np.swapaxes(logo,0,1)
        logo = logo[:,::-1,:]
        imageItem.setImage(logo)

    def addParameterTree(self):
        # Define parameter tree  --------------------------------------------
        self.params = [
            {'name': 'Mouse Behavior', 'type': 'group', 'children': [
                {'name': 'Select Point', 'type': 'action'},
                {'name': 'Select Region', 'type': 'action'},
            ]},
            {'name': 'Map View', 'type': 'group', 'children': [
                {'name': 'Update Map', 'type': 'action'},
                {'name': 'Map Scale', 'type': 'group', 'children': [
                    {'name': 'Max', 'type': 'float', 'value': self.mapCurrentMax, 'default': self.mapGlobalMax},
                    {'name': 'Min', 'type': 'float', 'value': self.mapCurrentMin, 'default': self.mapGlobalMin},
                ]},
                {'name': 'Map Bounds', 'type': 'group', 'children': [
                    {'name': 'Max Longitude', 'type': 'float', 'value': self.eastLon, 'default': self.lonMax},
                    {'name': 'Min Longitude', 'type': 'float', 'value': self.westLon, 'default': self.lonMin},
                    {'name': 'Max Latitude', 'type': 'float', 'value': self.northLat, 'default': self.latMax},
                    {'name': 'Min Latitude', 'type': 'float', 'value': self.southLat, 'default': self.latMin},
                ]},
                {'name': 'Restore Map Bounds', 'type': 'action'},
            ]},
            {'name': 'Map Annotation', 'type': 'group', 'expanded': False, 'children': [
                {'name': 'Coastlines', 'type': 'bool', 'value': True},
                {'name': 'Rivers', 'type': 'bool', 'value': False},
                {'name': 'Countries', 'type': 'bool', 'value': False},
                {'name': 'US States', 'type': 'bool', 'value': False},
                
            ]},
            {'name': 'Data Transformations', 'type': 'group', 'expanded': False, 'children': [
                {'name': 'Smoothing', 'type': 'group', 'expanded': False, 'children': [
                    {'name': 'Bilinear', 'type': 'int', 'value': self.bilinResMultiplier, 'default': 1, 'limits': (1,self.maxBRMultiplier)},
                    {'name': 'Gaussian', 'type': 'int', 'value': self.gaussWidth, 'default': 0, 'limits': (0,self.maxGaussWidth)},
                ]},
                {'name': 'Log', 'type': 'bool', 'value': False},
            ]},
            {'name': 'Time Series View', 'type': 'group', 'children': [
                {'name': 'Update Time Series', 'type': 'action'},
                {'name': 'Time Series Location', 'type': 'group', 'children': [
                    {'name': 'X', 'type': 'int', 'value': self.xTimeSeries},
                    {'name': 'Y', 'type': 'int', 'value': self.yTimeSeries},
                ]},
                {'name': 'Time Series Bounds', 'type': 'group', 'children': [
                    {'name': 'Max Y-Axis', 'type': 'float', 'value': self.tsCurrentYMax, 'default': self.tsGlobalYMax},
                    {'name': 'Min Y-Axis', 'type': 'float', 'value': self.tsCurrentYMin, 'default': self.tsGlobalYMin},
                    {'name': 'Max X-Axis', 'type': 'float', 'value': self.tsCurrentXMax, 'default': self.tsGlobalXMax},
                    {'name': 'Min X-Axis', 'type': 'float', 'value': self.tsCurrentXMin, 'default': self.tsGlobalXMin},
                ]},
                {'name': 'Restore Time Series Bounds', 'type': 'action'},
            ]},
            {'name': 'Time Series Analysis', 'type': 'group', 'children': [
                {'name': 'Calculate Decomposition', 'type': 'bool', 'value': False},
                {'name': 'Mean', 'type': 'float', 'value': self.tsMean, 'readonly':True},
                {'name': 'Trend','type': 'float', 'value': self.tsTrend, 'readonly':True},
                {'name': 'Amplitude Annual Cycle','type': 'float', 'value': self.tsAmp, 'readonly':True},
                {'name': 'Plot Mean', 'type': 'bool', 'value': False},
                {'name': 'Plot Trend', 'type': 'bool', 'value': False},
                {'name': 'Plot Annual Cycle', 'type': 'bool', 'value': False},
            ]},
            {'name': 'Export Image', 'type': 'group', 'expanded': False, 'children': [
                {'name': 'Map', 'type': 'bool', 'value':True},
                {'name': 'Time Series', 'type': 'bool', 'value': True},
                {'name': 'Export', 'type': 'action'},
            ]},
        ]

        # Create tree of Parameter objects
        self.param_tree = Parameter.create(name='params', type='group', children=self.params)
        self.ptwidget.setParameters(self.param_tree, showTop=False)

        # Define actions taken when action buttons selected
        self.param_tree.param('Mouse Behavior', 'Select Point').sigActivated.connect(self.selectPoint)
        self.param_tree.param('Mouse Behavior', 'Select Region').sigActivated.connect(self.selectRegion)
        
        self.param_tree.param('Map View', 'Update Map').sigActivated.connect(self.updateMap)

        self.param_tree.param('Map View', 'Restore Map Bounds').sigActivated.connect(self.restoreMapBounds)

        # Define actions taken when parameters changed
        self.param_tree.param('Map View', 'Map Scale','Min').sigValueChanged.connect(self.updateMapScaleMin)
        self.param_tree.param('Map View', 'Map Scale','Max').sigValueChanged.connect(self.updateMapScaleMax)

        self.param_tree.param('Map View', 'Map Bounds','Min Longitude').sigValueChanged.connect(self.updateMapWestLon)
        self.param_tree.param('Map View', 'Map Bounds','Max Longitude').sigValueChanged.connect(self.updateMapEastLon)
        self.param_tree.param('Map View', 'Map Bounds','Min Latitude').sigValueChanged.connect(self.updateMapSouthLat)
        self.param_tree.param('Map View', 'Map Bounds','Max Latitude').sigValueChanged.connect(self.updateMapNorthLat)

        self.param_tree.param('Map Annotation', 'Coastlines').sigValueChanged.connect(self.updateMapCoasts)
        self.param_tree.param('Map Annotation', 'Rivers').sigValueChanged.connect(self.updateMapRivers)
        self.param_tree.param('Map Annotation', 'Countries').sigValueChanged.connect(self.updateMapCountries)
        self.param_tree.param('Map Annotation', 'US States').sigValueChanged.connect(self.updateMapStates)

        self.param_tree.param('Data Transformations', 'Smoothing', 'Bilinear').sigValueChanged.connect(self.updateBilinearSmoothing)
        self.param_tree.param('Data Transformations', 'Smoothing', 'Gaussian').sigValueChanged.connect(self.updateGaussianSmoothing)
        self.param_tree.param('Data Transformations', 'Log').sigValueChanged.connect(self.updateLogTransform)

        self.param_tree.param('Time Series View', 'Update Time Series').sigActivated.connect(self.updateTimeSeries)
        self.param_tree.param('Time Series View', 'Restore Time Series Bounds').sigActivated.connect(self.restoreTimeSeriesBounds)
        self.param_tree.param('Time Series View', 'Time Series Location','X').sigValueChanged.connect(self.updateXTimeSeries)
        self.param_tree.param('Time Series View', 'Time Series Location','Y').sigValueChanged.connect(self.updateYTimeSeries)
        self.param_tree.param('Time Series View', 'Time Series Bounds','Min Y-Axis').sigValueChanged.connect(self.updateTsYMin)
        self.param_tree.param('Time Series View', 'Time Series Bounds','Max Y-Axis').sigValueChanged.connect(self.updateTsYMax)
        self.param_tree.param('Time Series View', 'Time Series Bounds','Min X-Axis').sigValueChanged.connect(self.updateTsXMin)
        self.param_tree.param('Time Series View', 'Time Series Bounds','Max X-Axis').sigValueChanged.connect(self.updateTsXMax)

        self.param_tree.param('Time Series Analysis', 'Calculate Decomposition').sigValueChanged.connect(self.updateTSDecomposition)

        self.param_tree.param('Export Image', 'Export').sigActivated.connect(self.exportImage)

        # Define functions of parameter tree items

    def getViewRange(self):
        x=self.vb_image.viewRange()
        lonmin = np.max((self.lonMin,np.min((x[0][0],self.lonMax))))
        lonmax = np.max((self.lonMin,np.min((x[0][1],self.lonMax))))
        latmin = np.max((self.latMin,np.min((x[1][0],self.latMax))))
        latmax = np.max((self.latMin,np.min((x[1][1],self.latMax))))
        return [(lonmin, lonmax),(latmin,latmax)]

    def dragUpdateRange(self):
        x=self.getViewRange()[0]
        y=self.getViewRange()[1]     
        self.updateMapBounds(x,y)
        
    def mouseMovedMap(self,evt):
        if self.doSelectPoint:
            pos = evt[0]  # using signal proxy turns original arguments into a tuple
            if self.vb_image.sceneBoundingRect().contains(pos):
                mousePoint = self.vb_image.mapSceneToView(pos)
                mx = np.int(mousePoint.x()/self.mapLonScale - self.mapLonTrans)
                my = np.int(mousePoint.y()/self.mapLatScale - self.mapLatTrans)
                if np.logical_and((mx > 0 and mx < self.im),
                                  (my > 0 and my < self.jm)):
                    self.top_data_label_right.setText("<span style='font-size: 12pt' style='color: black'>x=%0.1f,y=%0.1f,z=%0.1f</span>" % (self.mapLon[mx], self.mapLat[my],self.mapData[mx,my]))
                    self.vLine_top.setPos(mousePoint.x())
                    self.hLine_top.setPos(mousePoint.y())

        if self.doSelectRegion:
            self.mapROI.sigRegionChanged.connect(self.updateMapBoundsFromROI)

        if not self.doSelectPoint and not self.doSelectRegion:            
            self.vb_image.sigRangeChanged.connect(self.dragUpdateRange)

    def mouseMovedTimeSeries(self,evt):
        if self.doSelectPoint:
            pos = evt[0]  # using signal proxy turns original arguments into a tuple
            if self.time_series_plot.sceneBoundingRect().contains(pos):
                mousePoint = self.time_series_plot.vb.mapSceneToView(pos)
                index = int(mousePoint.x())
                if index > 0 and index < len(self.tsTimeAll):
                    self.bottom_data_label_right.setText("<span style='color: black' style='font-size: 12pt'>time=%0.1f,y=%0.1f</span>" % (mousePoint.x(), mousePoint.y()))
                    self.vLine_bot.setPos(mousePoint.x())
                    self.hLine_bot.setPos(mousePoint.y())

    # Update bottom plot when mouse clicked in top plot
    # i.e. choose a location in the map, then plot its time series
    def mouseClickedMap(self,evt):
        if self.doSelectPoint:
            pos = evt[0].scenePos()
            if self.vb_image.sceneBoundingRect().contains(pos):
                mousePoint = self.vb_image.mapSceneToView(pos)
                mx = np.int(mousePoint.x()/self.mapLonScale - self.mapLonTrans)
                my = np.int(mousePoint.y()/self.mapLatScale - self.mapLatTrans)
                if np.logical_and((mx > 0 and mx < self.im),(my > 0 and my < self.jm)):
                    self.xTimeSeries = mx
                    self.yTimeSeries = my
                    self.setTimeSeriesData()
                    # clear if left-clicked, overplot if right-clicked
                    if evt[0].button() > 1:
                        clear = False
                    else:
                        clear = True
                    self.plotTimeSeries(clear=clear)
                    self.bottom_data_label_left.setText("<span style='font-size: 12pt' style='color: black'>lon=%0.2f,  lat=%0.2f</span>" % (self.mapLon[mx],self.mapLat[my]))

    # Update top plot when mouse clicked in bottom plot
    # i.e. choose a time from the time series, and plot the map at that time
    def mouseClickedTimeSeries(self,evt):
        if self.doSelectPoint:
            pos = evt[0].scenePos()
            if self.time_series_plot.sceneBoundingRect().contains(pos):
                mousePoint = self.time_series_plot.vb.mapSceneToView(pos)
                mx = np.argmin(np.abs(mousePoint.x() - self.tsTimeAll))
                if (mx > 0 and mx < len(self.tsTimeAll)):
                    ptime = self.tsTimeAll[mx]
                    self.mapTime = mx
                    self.updateMap()
                    self.top_data_label_left.setText("<span style='font-size: 12pt' style='color: black'>time=%0.2f</span>" % (ptime))

    def selectPoint(self):
        if self.doSelectRegion:
            self.deselectRegion()

        if not self.doSelectPoint:
            # create crosshairs for each plot
            self.vLine_top = pg.InfiniteLine(angle=90, movable=False)
            self.hLine_top = pg.InfiniteLine(angle=0, movable=False)
            self.vb_image.addItem(self.vLine_top, ignoreBounds=True)
            self.vb_image.addItem(self.hLine_top, ignoreBounds=True)
            self.vLine_bot = pg.InfiniteLine(angle=90, movable=False)
            self.hLine_bot = pg.InfiniteLine(angle=0, movable=False)
            self.time_series_plot.addItem(self.vLine_bot, ignoreBounds=True)
            self.time_series_plot.addItem(self.hLine_bot, ignoreBounds=True)

            self.doSelectPoint = True
        else:
            self.deselectPoint()
                    
    def deselectPoint(self):
        self.doSelectPoint = False
        self.vb_image.removeItem(self.vLine_top)
        self.vb_image.removeItem(self.hLine_top)
        self.time_series_plot.removeItem(self.vLine_bot)
        self.time_series_plot.removeItem(self.hLine_bot)

    def selectRegion(self):

        if self.doSelectPoint:
            self.deselectPoint()

        if not self.doSelectRegion:
            # create region of interest
            xmean = np.mean((self.westLon,self.eastLon))
            ymean = np.mean((self.southLat,self.northLat))
            xsize = 0.25*(self.eastLon - self.westLon)
            ysize = 0.25*(self.northLat - self.southLat)
            xpos, ypos = xmean - 0.5*xsize, ymean - 0.5*ysize
            # get bounds
            bnds = self.vb_image.itemBoundingRect(self.map_image)
            self.mapROI = pg.ROI([xpos,ypos],[xsize,ysize],
                                 maxBounds=bnds,pen=(0,0,0))
            self.mapROI.setZValue(15)
            self.mapROI.addScaleHandle([0,1],[0.5,0.5])
            self.mapROI.addScaleHandle([1,0],[0.5,0.5])
            self.vb_image.addItem(self.mapROI)
            self.doSelectRegion = True
        else:
            self.deselectRegion()

    def deselectRegion(self):
        self.doSelectRegion = False
        self.vb_image.removeItem(self.mapROI)

    def updateMap(self):
        # Update map image and annotations
        if self.doSelectRegion:
            self.deselectRegion()

        self.plotMap(rescale=True)

        self.map_overplot_coasts.clear()
        self.updateMapCoasts()

        self.map_overplot_rivers.clear()
        self.updateMapRivers()

        self.map_overplot_countries.clear()
        self.updateMapCountries()

        self.map_overplot_states.clear()
        self.updateMapStates()

    def updateMapBoundsFromROI(self):
        x=self.mapROI.getAffineSliceParams(self.mapData,self.map_image)
        # x[0] is shape, x[2] is origin
        x0 = np.int(x[2][0])
        y0 = np.int(x[2][1])
        x1 = x0 + np.int(x[0][0])
        y1 = y0 + np.int(x[0][1])
        self.updateMapBounds((self.mapLon[x0],self.mapLon[x1]),
                             (self.mapLat[y0],self.mapLat[y1]))

    def updateMapBounds(self,x,y):
        #update parameter tree values
        self.param_tree.param('Map View', 'Map Bounds','Min Longitude').setValue(x[0])
        self.param_tree.param('Map View', 'Map Bounds','Max Longitude').setValue(x[1])
        self.param_tree.param('Map View', 'Map Bounds','Min Latitude').setValue(y[0])
        self.param_tree.param('Map View', 'Map Bounds','Max Latitude').setValue(y[1])
        self.updateMapWestLon
        self.updateMapEastLon
        self.updateMapSouthLat
        self.updateMapNorthLat

    def restoreMapBounds(self):
        # Restore bounds to global values
        self.param_tree.param('Map View', 'Map Bounds','Min Longitude').setValue(self.lonMin)
        self.param_tree.param('Map View', 'Map Bounds','Max Longitude').setValue(self.lonMax)
        self.param_tree.param('Map View', 'Map Bounds','Min Latitude').setValue(self.latMin)
        self.param_tree.param('Map View', 'Map Bounds','Max Latitude').setValue(self.latMax)

    def updateMapScaleMin(self):
        newVal = self.param_tree.param('Map View', 'Map Scale','Min').value()
        self.mapCurrentMin = np.min((self.mapGlobalMax,np.max((self.mapGlobalMin,newVal))))

    def updateMapScaleMax(self):
        newVal = self.param_tree.param('Map View', 'Map Scale','Max').value()
        self.mapCurrentMax = np.min((self.mapGlobalMax,np.max((self.mapGlobalMin,newVal))))

    def updateMapWestLon(self):
        newVal = self.param_tree.param('Map View', 'Map Bounds','Min Longitude').value()
        self.westLon = np.min((self.lonMax,np.max((self.lonMin,newVal))))

    def updateMapEastLon(self):
        newVal = self.param_tree.param('Map View', 'Map Bounds','Max Longitude').value()
        self.eastLon = np.min((self.lonMax,np.max((self.lonMin,newVal))))

    def updateMapSouthLat(self):
        newVal = self.param_tree.param('Map View', 'Map Bounds','Min Latitude').value()
        self.southLat = np.min((self.latMax,np.max((self.latMin,newVal))))

    def updateMapNorthLat(self):
        newVal = self.param_tree.param('Map View', 'Map Bounds','Max Latitude').value()
        self.northLat = np.min((self.latMax,np.max((self.latMin,newVal))))

    def updateMapRivers(self):
        if self.param_tree.param('Map Annotation', 'Rivers').value():
            self.plot_rivers()
        else:
            self.map_overplot_rivers.clear()
            
    def updateMapCountries(self):
        if self.param_tree.param('Map Annotation', 'Countries').value():
            self.plot_countries()
        else:
            self.map_overplot_countries.clear()

    def updateMapStates(self):
        if self.param_tree.param('Map Annotation', 'US States').value():
            self.plot_states()
        else:
            self.map_overplot_states.clear()

    def updateMapCoasts(self):
        if self.param_tree.param('Map Annotation', 'Coastlines').value():
            self.plot_coastlines()
        else:
            self.map_overplot_coasts.clear()
            
    def plot_coastlines(self):
        lonbounds=np.asarray((self.westLon,self.eastLon),dtype=np.float)
        latbounds=np.asarray((self.southLat,self.northLat),dtype=np.float)
        libcivet.plot_coastlines_nc(self.map_overplot_coasts,
                                center='Dateline',
                                pen=pg.mkPen(self.coastPenColor,
                                width=self.coastPenWidth),
                                xbounds=lonbounds,ybounds=latbounds)

    def plot_rivers(self):
        lonbounds=np.asarray((self.westLon,self.eastLon),dtype=np.float)
        latbounds=np.asarray((self.southLat,self.northLat),dtype=np.float)
        libcivet.plot_rivers_nc(self.map_overplot_rivers,
                                center='Dateline',
                                pen=pg.mkPen(self.riverPenColor,
                                width=self.riverPenWidth),
                                xbounds=lonbounds,ybounds=latbounds)

    def plot_countries(self):
        lonbounds=np.asarray((self.westLon,self.eastLon),dtype=np.float)
        latbounds=np.asarray((self.southLat,self.northLat),dtype=np.float)
        libcivet.plot_countries_nc(self.map_overplot_countries,
                                center='Dateline',
                                pen=pg.mkPen(self.countryPenColor,
                                width=self.countryPenWidth),
                                xbounds=lonbounds,ybounds=latbounds)

    def plot_states(self):
        lonbounds=np.asarray((self.westLon,self.eastLon),dtype=np.float)
        latbounds=np.asarray((self.southLat,self.northLat),dtype=np.float)
        libcivet.plot_states_nc(self.map_overplot_states,
                                center='Dateline',
                                pen=pg.mkPen(self.statePenColor,
                                width=self.statePenWidth),
                                xbounds=lonbounds,ybounds=latbounds)

    def updateBilinearSmoothing(self):
        newVal = self.param_tree.param('Data Transformations', 'Smoothing', 'Bilinear').value()
        self.bilinResMultiplier = np.min((self.maxBRMultiplier,np.max((1,newVal))))

    def updateGaussianSmoothing(self):
        newVal = self.param_tree.param('Data Transformations', 'Smoothing', 'Gaussian').value()
        self.gaussWidth = np.min((self.maxGaussWidth,np.max((0,newVal))))

    def updateLogTransform(self):
        self.logTransform = self.param_tree.param('Data Transformations', 'Log').value()

    def updateTimeSeries(self):
        # Update time series
        self.plotTimeSeries()

    def restoreTimeSeriesBounds(self):
        # Restore bounds to global values
        self.param_tree.param('Time Series View', 'Time Series Bounds','Min X-Axis').setValue(self.tsGlobalXMin)
        self.param_tree.param('Time Series View', 'Time Series Bounds','Max X-Axis').setValue(self.tsGlobalXMax)
        self.param_tree.param('Time Series View', 'Time Series Bounds','Min Y-Axis').setValue(self.tsGlobalYMin)
        self.param_tree.param('Time Series View', 'Time Series Bounds','Max Y-Axis').setValue(self.tsGlobalYMax)
        if len(self.range_time_plot.items) > 0:
            self.range_time_plot.items[1].setRegion((self.tsGlobalXMin, self.tsGlobalXMax))
        if len(self.range_data_plot.items) > 0:
            self.range_data_plot.items[1].setRegion((self.tsGlobalYMin, self.tsGlobalYMax))
        
    def updateXTimeSeries(self):
        newVal = self.param_tree.param('Time Series View', 'Time Series Location','X').value()
        self.xTimeSeries = np.min((self.im,np.max((0,newVal))))
        self.tsYMin = np.min(self.data[self.xTimeSeries,self.yTimeSeries,:])
        self.tsYMax = np.max(self.data[self.xTimeSeries,self.yTimeSeries,:])
        self.restoreTimeSeriesBounds()
        
    def updateYTimeSeries(self):
        newVal = self.param_tree.param('Time Series View', 'Time Series Location','Y').value()
        self.yTimeSeries = np.min((self.jm,np.max((0,newVal))))
        self.tsYMin = np.min(self.data[self.xTimeSeries,self.yTimeSeries,:])
        self.tsYMax = np.max(self.data[self.xTimeSeries,self.yTimeSeries,:])
        self.restoreTimeSeriesBounds()

    def updateTsYMin(self):
        newVal = self.param_tree.param('Time Series View', 'Time Series Bounds','Min Y-Axis').value()
        self.tsCurrentYMin = np.min((self.tsGlobalYMax,np.max((self.tsGlobalYMin,newVal))))

    def updateTsYMax(self):
        newVal = self.param_tree.param('Time Series View', 'Time Series Bounds','Max Y-Axis').value()
        self.tsCurrentYMax = np.min((self.tsGlobalYMax,np.max((self.tsGlobalYMin,newVal))))

    def updateTsXMin(self):
        newVal = self.param_tree.param('Time Series View', 'Time Series Bounds','Min X-Axis').value()
        self.tsCurrentXMin = np.min((self.tsGlobalXMax,np.max((self.tsGlobalXMin,newVal))))

    def updateTsXMax(self):
        newVal = self.param_tree.param('Time Series View', 'Time Series Bounds','Max X-Axis').value()
        self.tsCurrentXMax = np.min((self.tsGlobalXMax,np.max((self.tsGlobalXMin,newVal))))

    def exportImage(self):
        exportMap = self.param_tree.param('Export Image', 'Map').value()
        exportTimeSeries = self.param_tree.param('Export Image', 'Time Series').value()
        if exportMap:
            exporter = pg.exporters.ImageExporter(self.top_layout)
            exporter.export('map.civet.png')
        if exportTimeSeries:
            exporter = pg.exporters.ImageExporter(self.bottom_layout)
            exporter.export('time_series.civet.png')
            
    # end of parameter tree functionality

    def addGlobalRange(self):
        # time will be above main plot
        self.range_time_plot.showAxis('left',False)
        self.range_time_plot.showAxis('bottom',False)
        self.range_time_plot.showAxis('top',True)
        xbounds=[np.min(self.time),np.max(self.time)]
        self.range_time_plot.setLimits(xMin=xbounds[0],xMax=xbounds[1])

        lr = pg.LinearRegionItem(xbounds,pen=pg.mkPen(width=5))
        lr.setOpacity(0.33)
        lr.setZValue(20)
        self.range_time_plot.plot(self.time,self.time,pen=None)
        self.range_time_plot.addItem(lr)

        self.range_time_plot.disableAutoRange()
        self.range_time_plot.hideButtons()

        # data will be to right of main plot
        self.range_data_plot.showAxis('left',False)
        self.range_data_plot.showAxis('bottom',False)
        self.range_data_plot.showAxis('right',True)
        data=self.data[self.xTimeSeries,self.yTimeSeries,:]
        ybounds=[self.mapGlobalMin,self.mapGlobalMax]
        self.range_data_plot.setLimits(yMin=ybounds[0],yMax=ybounds[1])

        lr = pg.LinearRegionItem(ybounds,orientation='horizontal',
                                 pen=pg.mkPen(width=5))
        lr.setOpacity(0.33)
        lr.setZValue(20)
        x=(np.min(data),np.max(data))
        self.range_data_plot.plot((0,1),ybounds,pen=None)
        self.range_data_plot.addItem(lr)

        self.range_data_plot.disableAutoRange()
        self.range_data_plot.hideButtons()

        self.range_time_plot.items[1].sigRegionChanged.connect(self.updateTimeRange)
        self.range_data_plot.items[1].sigRegionChanged.connect(self.updateDataRange)

    def updateTimeRange(self):
        x=self.range_time_plot.items[1].getRegion()
        #update parameter tree values
        self.param_tree.param('Time Series View', 'Time Series Bounds','Min X-Axis').setValue(x[0])
        self.param_tree.param('Time Series View', 'Time Series Bounds','Max X-Axis').setValue(x[1])
        self.updateTsXMin
        self.updateTsXMax
        
    def updateDataRange(self):
        x=self.range_data_plot.items[1].getRegion()
        #update parameter tree values
        self.param_tree.param('Time Series View', 'Time Series Bounds','Min Y-Axis').setValue(x[0])
        self.param_tree.param('Time Series View', 'Time Series Bounds','Max Y-Axis').setValue(x[1])
        self.updateTsYMin
        self.updateTsYMax
        
    def updateTSDecomposition(self):
        if self.param_tree.param('Time Series Analysis', 'Calculate Decomposition').value():
            self.doTSDecomposition = True
        else:
            self.doTSDecomposition = False
            self.param_tree.param('Time Series Analysis', 'Mean').setValue(-999)
            self.param_tree.param('Time Series Analysis', 'Trend').setValue(-999)
            self.param_tree.param('Time Series Analysis', 'Amplitude Annual Cycle').setValue(-999)
        
    def calcTSDecomposition(self):
        x=self.tsTime
        y=self.tsData
        ncoefs = 4
        coefs=libcivet.fit_seasonal(x,y,ncoefs)
        self.tsMean = np.mean(self.tsData)
        self.tsTrend = coefs[3]
        self.tsAmp = np.sqrt(np.power(coefs[0],2)+np.power(coefs[1],2))
        self.param_tree.param('Time Series Analysis', 'Mean').setValue(self.tsMean)
        self.param_tree.param('Time Series Analysis', 'Trend').setValue(self.tsTrend)
        self.param_tree.param('Time Series Analysis', 'Amplitude Annual Cycle').setValue(self.tsAmp)
        return coefs

    def createLUT(self):
        # create colormap (directly calling the cm object returns ndarray of RGBA)
        # the point of this is to extract the matplotlib colormap rgb
        # values, which are then used to create a pyqtgraph colormap
        #x=plt.cm.viridis
        x=plt.cm.nipy_spectral
        #x=plt.cm.rainbow
        #x=plt.cm.gist_rainbow
        #x=plt.cm.bwr
        ncols=x.N
        self.lut=np.asarray(255 * x(range(ncols)),dtype=np.uint8) #(256,4), values [0-1]

        # subset number of contours
        if 1==2:
            cols=x(range(ncols))
            cols2=cols[1:ncols:16,]
            ncols2=cols2.shape[0]
            lut=np.asarray(255 * cols2,dtype=np.uint8) #(256,4), values [0-1]

        # set first level to white
        self.lut[0,:] = (255,255,255,255)

    def restore_map_bounds(self):
        self.westLon  = self.lonMin
        self.eastLon  = self.lonMax
        self.southLat = self.latMin
        self.northLat = self.latMax
        
    def readData(self,input_filename,noFillValue=True,decimalTime=False):
        htag = 'clm2.h0'
        self.varname=input_filename.split(htag)[-1].split('.')[1]

        f1   = netcdf4.Dataset(input_filename, 'r')
        self.data = np.copy(f1.variables[self.varname])
        self.yUnits = f1.variables[self.varname].getncattr('units')
        self.lon  = np.copy(f1.variables['lon'])
        self.lat  = np.copy(f1.variables['lat'])
        self.time = np.copy(f1.variables['time'])
        self.xUnits = f1.variables['time'].getncattr('units')
        self.dataRefYear = np.float(self.xUnits.split()[2].split('-')[0])
        if noFillValue:
            fillValue = np.float(f1.variables[self.varname].getncattr('_FillValue'))
        self.im=self.lon.shape[0]
        self.jm=self.lat.shape[0]
        f1.close()
        # swap time and longitude to give (lon,lat,time)
        self.data = np.swapaxes(self.data,0,2)

        valid_data_ind=np.where(self.data < fillValue)
        # remove fill values
        if noFillValue:
            ind=np.where(self.data >= fillValue)
            self.data[ind] = 0.0

        # convert to decimal year time
        # (ok for no_leap calendar; wrong for gregorian calendar)
        if decimalTime:
            self.time = self.time/365.0 + self.dataRefYear
            self.xUnits = 'year'

        self.dlon=np.abs(self.lon[0]-self.lon[1])
        self.dlat=np.abs(self.lat[0]-self.lat[1])

        # set maximum geographical coordinates
        self.lonMin = np.min((360.0,np.max((0.0,np.min(self.lon)))))
        self.lonMax = np.min((360.0,np.max((0.0,np.max(self.lon)))))
        self.latMin = np.min((90.0,np.max((-90.0,np.min(self.lat)))))
        self.latMax = np.min((90.0,np.max((-90.0,np.max(self.lat)))))

        # set data limits
        self.mapGlobalMin  = np.min(self.data[valid_data_ind])
        self.mapGlobalMax  = np.max(self.data[valid_data_ind])
        self.mapCurrentMin  = self.mapGlobalMin
        self.mapCurrentMax  = self.mapGlobalMax
        self.contour_levels = (self.mapCurrentMin,self.mapCurrentMax)

        self.tsGlobalXMin = np.min(self.time)
        self.tsGlobalXMax = np.max(self.time)
        self.tsGlobalYMin = np.min(self.data[valid_data_ind])
        self.tsGlobalYMax = np.max(self.data[valid_data_ind])
        self.tsLocalYMin = np.min(self.data[self.xTimeSeries,self.yTimeSeries,:])
        self.tsLocalYMax = np.max(self.data[self.xTimeSeries,self.yTimeSeries,:])
        self.tsCurrentXMin = self.tsGlobalXMin
        self.tsCurrentXMax = self.tsGlobalXMax
        self.tsCurrentYMin = self.tsGlobalYMin
        self.tsCurrentYMax = self.tsGlobalYMax
        
        # Update parameter tree values
        self.param_tree.param('Map View', 'Map Scale','Min').setValue(self.mapCurrentMin)
        self.param_tree.param('Map View', 'Map Scale','Min').setDefault(self.mapGlobalMin)
        self.param_tree.param('Map View', 'Map Scale','Max').setValue(self.mapCurrentMax)
        self.param_tree.param('Map View', 'Map Scale','Max').setDefault(self.mapGlobalMax)
        self.restoreMapBounds()
        self.restoreTimeSeriesBounds()
        
    def create_colorbar(self,cbx=None,cby=None):
        if cbx == None:
            cbx = 10
        if cby == None:
            cby = 256
        self.cb=np.zeros((cbx,cby))
        for i in range(cbx):
            self.cb[i,:] = np.linspace(0,cby-1,cby)

        self.colorbar.setImage(self.cb)
        self.colorbar.setLookupTable(self.lut)
        #self.scale_colorbar()

    def scale_colorbar(self,rescale=False):
        # remove current scaling prior to new scaling
        if rescale:
            # assumes plot has been previously scaled
            self.colorbar.translate(0,-self.cbTrans)
            self.colorbar.scale(1,(1./self.cbScale))
            
        self.cbScale = (self.mapCurrentMax-self.mapCurrentMin)/self.cb.shape[1]
        self.cbTrans = self.mapCurrentMin / self.cbScale
        self.colorbar.scale(1,self.cbScale)
        self.colorbar.translate(0,self.cbTrans)

    def applyBilinearInterpolation(self):
        ''' bilinResMultiplier is a multiplier that is
        applied to the native resolution '''
        if self.bilinResMultiplier > 1:
            
            from scipy import interpolate
        
            x = interpolate.interp2d(self.mapLat,self.mapLon,self.mapData, kind='linear')
            im = np.int(self.mapLon.size * self.bilinResMultiplier)
            blon = (self.lonMax - self.lonMin) \
                   * np.arange(im)/np.float(im-1) + self.lonMin

            jm = np.int(self.mapLat.size * self.bilinResMultiplier)
            blat = (self.latMax - self.latMin) \
                   * np.arange(jm)/np.float(jm-1) + self.latMin

            fld =x(blat,blon)
            self.mapData = fld
            self.mapLon = blon
            self.mapLat = blat
        
    def applyGaussian(self):
        if self.gaussWidth > 0:
            hwx = np.max((1,np.min((self.gaussWidth,self.maxGaussWidth))))
            hwy = hwx
            fld = pg.gaussianFilter(self.mapData, (hwx,hwy))
            self.mapData = fld
        
    def applyLogTransform(self):
        # apply to absolute value of data
        if self.logTransform:
            fld = np.abs(np.copy(self.mapData))
            ind=np.where(fld > 0)
            fld[ind] = np.log(fld[ind])
            self.mapData = fld
        
    def plotMap(self,rescale=False):
        # Update contour levels
        self.contour_levels = (self.mapCurrentMin,self.mapCurrentMax)

        # Rescale colorbar
        self.scale_colorbar(rescale=rescale)
        # Subset map image
        self.setMapData()

        # Apply Bilinear Interpolation
        self.applyBilinearInterpolation()
        # Apply Gaussian smoothing
        self.applyGaussian()
        # Apply Logarithmic transform (to absolute value of data)
        self.applyLogTransform()

        self.map_image.setLookupTable(self.lut)
        self.map_image.setImage(self.mapData,levels=self.contour_levels)        

        self.scale_map_boundaries(rescale=rescale)

    def setMapData(self):
        self.mapData=self.data[:,:,self.mapTime]
        self.mapLon=self.lon
        self.mapLat=self.lat
        
    def scale_map_boundaries(self,scale=True,rescale=False):
        # remove current scaling prior to new scaling
        if rescale:
            # assumes plot has been previously scaled
            # order of operations matters, i.e. reverse of scale/translate
            # is -translate/(1/scale)
            self.map_image.translate(-self.mapLonTrans,-self.mapLatTrans)
            self.map_image.scale((1.0/self.mapLonScale),(1.0/self.mapLatScale))

        if scale:
            self.mapLonScale = (self.lonMax-self.lonMin)/(self.mapLon.size-1)
            self.mapLonTrans = self.lonMin / self.mapLonScale
            self.mapLatScale = (self.latMax-self.latMin)/(self.mapLat.size-1)
            self.mapLatTrans = self.latMin / self.mapLatScale
            self.map_image.scale(self.mapLonScale,self.mapLatScale)
            self.map_image.translate(self.mapLonTrans,self.mapLatTrans)
            self.vb_image.setRange(yRange=[self.southLat,self.northLat],
                                   xRange=[self.westLon, self.eastLon])

    def axisLength(self,axisItem):
        # return size in pixels of an axis
        # copied from AxisItem.generateDrawSpecs
        bounds = axisItem.mapRectFromParent(axisItem.geometry())
        if axisItem.orientation == 'left':
            span = (bounds.topRight(), bounds.bottomRight())
        elif axisItem.orientation == 'right':
            span = (bounds.topLeft(), bounds.bottomLeft())
        elif axisItem.orientation == 'top':
            span = (bounds.bottomLeft(), bounds.bottomRight())
        elif axisItem.orientation == 'bottom':
            span = (bounds.topLeft(), bounds.topRight())
        
        ## determine size of this item in pixels
        points = list(map(axisItem.mapToDevice, span))
        if None in points:
            return
        lengthInPixels = pg.Point(points[1] - points[0]).length()
        return lengthInPixels
    
    def replaceAxisTicks(self,axis_dst,axis_src,dstRange):
        # set destination axis strings to span dstRange
        # and plot at source axis values
        orient_dst = axis_dst.orientation
        orient_src = axis_src.orientation
        # could check orientation consistency

        # get tick spacing/values from source axis
        tv_src = axis_dst.tickValues(axis_src.range[0],axis_src.range[1],self.axisLength(axis_dst))

        # extract tick values and sort 
        val_src = np.asarray([])
        for n in tv_src:
            val_src = np.append(val_src, np.asarray(n[1]))
        val_src.sort()
        nticks = val_src.size

        # map from source to destination axis
        vscale = (val_src - axis_src.range[0]) / (axis_src.range[1] - axis_src.range[0])
        val_dst = dstRange[0] + (dstRange[1] - dstRange[0]) * vscale
        
        majTick = [(val_src[n],'{:.2f}'.format(val_dst[n])) for n in range(nticks)]
        #input to setTicks is tuple of (value,namestring)
        setTick = [majTick,[]]        
        axis_dst.setTicks(setTick)

    def autoBtnClicked(self):
        # this is occurring before autoscale, so to add a delay
        # set a flag that will be used by a subsequent function call
        self.doAutoScaleUpdateXY = True
                
    def autoScaleUpdateXY(self):
        x = self.time_series_plot.vb.viewRange()[0]
        y = self.time_series_plot.vb.viewRange()[1]
        self.param_tree.param('Time Series View', 'Time Series Bounds','Min X-Axis').setValue(x[0])
        self.param_tree.param('Time Series View', 'Time Series Bounds','Max X-Axis').setValue(x[1])
        self.param_tree.param('Time Series View', 'Time Series Bounds','Min Y-Axis').setValue(y[0])
        self.param_tree.param('Time Series View', 'Time Series Bounds','Max Y-Axis').setValue(y[1])

        self.doAutoScaleUpdateXY = False
       
    def plotTimeSeries(self,rescale=False,clear=True):
        self.setTimeSeriesData()
        if clear:
            self.time_series_plot.clear()
        #self.time_series_plot.plot(self.tsTimeAll,self.tsDataAll,pen=(0,0,0))
        self.time_series_plot.plot(self.tsTimeAll,self.tsDataAll,
                                   pen=pg.mkPen((0,0,0),width=self.tsPenWidth))
        self.ts_bot_left_axis.setLabel(self.yUnits)
        self.ts_bottom_axis.setLabel(self.xUnits)
        self.time_series_plot.showAxis('left',False)
        self.time_series_plot.showAxis('bottom',False)
        self.time_series_plot.showAxis('top',False)
        self.time_series_plot.showAxis('right',True)
        self.time_series_plot.getAxis('right').setStyle(showValues=False)
        
        self.time_series_plot.disableAutoRange()
        self.time_series_plot.getViewBox().disableAutoRange()
        self.time_series_plot.getViewBox().setLimits(xMin=self.tsGlobalXMin,
                                                     xMax=self.tsGlobalXMax,
                                                     yMin=self.tsGlobalYMin,
                                                     yMax=self.tsGlobalYMax)

        # need an intermediate function to get correct timing of button event
        # calls 
        self.time_series_plot.sigRangeChanged.connect(self.autoScaleUpdateXY)

        #self.scaleTimeSeriesAxes(rescale=True)
        self.updateTimeSeriesBounds()

        # Overplot time series decomposition
        if self.doTSDecomposition:
            coefs = self.calcTSDecomposition()
            tsMean = np.zeros((self.tsTimeAll.size)) ; tsMean[:]=self.tsMean
            if self.param_tree.param('Time Series Analysis', 'Plot Mean').value():
                self.time_series_plot.plot(self.tsTimeAll,tsMean,pen=(255,0,0))

            tsTrend = np.zeros((self.tsTimeAll.size)) ; tsTrend[:]=self.tsTimeAll*coefs[3]+coefs[2]
            if self.param_tree.param('Time Series Analysis', 'Plot Trend').value():
                self.time_series_plot.plot(self.tsTimeAll,tsTrend,pen=(0,255,0))

            tsAmp = np.zeros((self.tsTimeAll.size)) ; tsAmp[:]=libcivet.synth_seasonal(self.tsTimeAll,coefs)
            if self.param_tree.param('Time Series Analysis', 'Plot Annual Cycle').value():
                self.time_series_plot.plot(self.tsTimeAll,tsAmp,pen=(0,0,255))

    def setTimeSeriesData(self):
        self.tsDataAll=self.data[self.xTimeSeries,self.yTimeSeries,:]
        self.tsTimeAll=self.time
        # filter x-axis
        ind =np.where(np.logical_and((self.time >= self.tsCurrentXMin),(self.time <= self.tsCurrentXMax)))[0]
        self.tsData=self.data[self.xTimeSeries,self.yTimeSeries,ind]
        self.tsTime=self.time[ind]
        # filter y-axis
        ind2=np.where(np.logical_and((self.tsData >= self.tsCurrentYMin),(self.tsData <= self.tsCurrentYMax)))[0]
        self.tsData=self.tsData[ind2]
        self.tsTime=self.tsTime[ind2]
         
    def updateTimeSeriesBounds(self):
        self.time_series_plot.setXRange(self.tsCurrentXMin,self.tsCurrentXMax)
        self.time_series_plot.setYRange(self.tsCurrentYMin,self.tsCurrentYMax)


#===============  module function definitions  ===========================

# set command line inputs 
parser = argparse.ArgumentParser(description='CIVET')
parser.add_argument('input_filename', type=str, nargs='?',
                    help='input history file name')
args = parser.parse_args()
infile = args.input_filename

if infile != None:
    pg.mkQApp()
    
    win = CivetGUI()
    
    win.Civet(infile)

    # default screen resolution
    xs0 = 1000
    sf = 0.8
    ys0 = sf * xs0
    screen_resolution = QtGui.QDesktopWidget().screenGeometry()
    xScreen, yScreen = screen_resolution.width(), screen_resolution.height()
    rx = np.min((1.0,xScreen/xs0))
    ry = np.min((1.0,yScreen/ys0))
    rScreen = np.min((rx,ry))
    xs,ys = (rScreen*xs0, rScreen*ys0)
    win.resize(xs,ys)
    win.resizeGUI()
    
    proxy_top = pg.SignalProxy(win.vb_image.scene().sigMouseMoved, rateLimit=60, slot=win.mouseMovedMap)
    
    proxy = pg.SignalProxy(win.time_series_plot.scene().sigMouseMoved, rateLimit=60, slot=win.mouseMovedTimeSeries)
    
    proxy_top_click = pg.SignalProxy(win.vb_image.scene().sigMouseClicked, slot=win.mouseClickedMap)
    
    proxy_bottom_click = pg.SignalProxy(win.time_series_plot.scene().sigMouseClicked, slot=win.mouseClickedTimeSeries)
    
    win.show()
    
    ## Start Qt event loop unless running in interactive mode or using pyside.
    if __name__ == '__main__':
        import sys
        if (sys.flags.interactive != 1) or not hasattr(pg.QtCore, 'PYQT_VERSION'):
            pg.QtGui.QApplication.instance().exec_()

