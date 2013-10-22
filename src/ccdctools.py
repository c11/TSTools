# -*- coding: utf-8 -*-
# vim: set expandtab:ts=4
"""
/***************************************************************************
 CCDCTools
                                 A QGIS plugin
 Plotting & visualization tools for CCDC Landsat time series analysis
                              -------------------
        begin                : 2013-03-15
        copyright            : (C) 2013 by Chris Holden
        email                : ceholden@bu.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
# Import the PyQt and QGIS libraries
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import QgsMapToolEmitPoint

# Initialize Qt resources from file resources.py
import resources_rc

import os

from ccdc_config import CCDCConfig 
from ccdc_controller import Controller
from ccdc_controls import CCDCControls
from ccdc_plot import CCDCPlot
import ccdc_settings as settings

class CCDCTools:

    def __init__(self, iface):
        ### Optional stuff to move elsewhere... #TODO
        # Save reference to the QGIS interface
        self.iface = iface
        self.canvas = self.iface.mapCanvas()

        ### Location info - define these elsewhere
        self.location = os.getcwd()
        self.image_pattern = 'LND*'
        self.stack_pattern = '*stack'

        # self.location = '/home/ceholden/Dropbox/Work/Research/pyCCDC/Dataset/p012r031/images'
        # self.location = '/net/caseq/lcscratch/ceholden/p012r030/images'
        # self.location = '/net/caseq/lcscratch/ceholden/QGIS/p013r029/images/'
		# self.image_pattern = 'LND*'
        # self.stack_pattern = '*stack'

        ### Toolbar - map tool & config
        self.init_toolbar()

    ### Toolbar section
    def init_toolbar(self):
        ### MapTool button
        self.action = QAction(QIcon(':/plugins/ccdctools/icon.png'),
                              'CCDC Tool', self.iface.mainWindow())
        self.action.setCheckable(True)
        self.action.triggered.connect(self.set_tool)
        self.iface.addToolBarIcon(self.action)

        self.tool_ts = QgsMapToolEmitPoint(self.canvas)
        self.tool_ts.setAction(self.action)
        self.tool_ts.canvasClicked.connect(self.plot_request)
        
        ### Configuration button
        self.action_cfg = QAction(QIcon(':/plugins/ccdctools/icon.png'),
            'Configure', self.iface.mainWindow())
        self.action_cfg.triggered.connect(self.handle_show_config)
        self.iface.addToolBarIcon(self.action_cfg)

    def set_tool(self):
        self.canvas.setMapTool(self.tool_ts)
    
    def handle_show_config(self):
        print 'Show/hide config'
        self.config = CCDCConfig(self, 
                                 self.location, 
                                 self.image_pattern,
                                 self.stack_pattern)
        # Connect signals for okay/cancel buttons
        self.config.accepted.connect(self.config_accepted)
        self.config.canceled.connect(self.config_canceled)
        # Execute & show dialog
        self.config.exec_()

    def config_accepted(self):
        print 'Accepted config!'
        # Accept the new values
        self.location = str(self.config.location)
        self.image_pattern = str(self.config.image_pattern)
        self.stack_pattern = str(self.config.stack_pattern)
        # Now, update the time series
        success = self.controller.get_time_series(self.location, 
                                        self.image_pattern,
                                        self.stack_pattern)
        # Close and disconnect the window if successful
        if success:
            self.config.close()
            self.config.accepted.disconnect()
            self.config.canceled.disconnect()
        else:
            message = QMessageBox()
            message.critical(None, 'Configuration Error',
                             'Could not find time series',
                             QMessageBox.Ok)

    def config_canceled(self):
        print 'Canceled config!'
        self.config.accepted.disconnect()
        self.config.canceled.disconnect()
        self.config.close()

    ### Configuration and plot widgets section
    def init_controls(self):
        """
        Initialize and add signals to the left side control widget
        """
        print 'init_controls'
        # Create widget
        self.ctrl_widget = CCDCControls(self.iface)
        # Create dock and add control widget
        self.ctrl_dock = QDockWidget("CCDC Tools", self.iface.mainWindow())
        self.ctrl_dock.setObjectName("CCDC Tools")
        self.ctrl_dock.setWidget(self.ctrl_widget)
        # Add to iface
        self.iface.addDockWidget(Qt.LeftDockWidgetArea, self.ctrl_dock)

    def init_plotter(self):
        """
        Initialize and add signals to the bottom area plotting widget
        """
        # Create widget
        self.plot_widget = CCDCPlot(self.iface)
        # Create dock and add plot widget
        self.plot_dock = QDockWidget('CCDC Plot', self.iface.mainWindow())
        self.plot_dock.setObjectName('CCDC Plot')
        self.plot_dock.setWidget(self.plot_widget)
        # Add to iface
        self.iface.addDockWidget(Qt.BottomDockWidgetArea, self.plot_dock)

    def initGui(self):
        """
        Required method for Qt to load components. Also inits signal controller
        """
        self.init_controls()

        self.init_plotter()
        self.controller = Controller(self.ctrl_widget, self.plot_widget,
                                     self.iface)
                    
    def plot_request(self, pos, button=None):
        print 'Trying to fetch...'
        if self.canvas.layerCount() == 0 or pos is None:
            print 'Could not fetch...'
            return
        layer = self.canvas.currentLayer()
        if (layer == None or layer.isValid() == False or 
            layer.type() != QgsMapLayer.RasterLayer):
            print 'Invalid layer...'
            return

        # Check if position needs to be reprojected to layer CRS
        layerCrs = layer.crs()
        mapCrs = self.canvas.mapRenderer().destinationCrs()

        if not mapCrs == layerCrs and self.canvas.hasCrsTransformEnabled():
            crsTransform = QgsCoordinateTransform(mapCrs, layerCrs)
            try:
                pos = crsTransform.transform(pos)
            except QgsCsException, err:
                print 'Transformation error'
                pass #TODO handle better?

        # If layer has position, get data
        if layer and layer.extent().contains(pos):
            self.controller.fetch_data(pos)
            self.controller.update_display()
            if settings.canvas['show_click']:
                self.controller.show_click(pos)

    def unload(self):
        """
        Handle startup/shutdown/hide/etc behavior
        """
        # Close toolbars
        self.iface.removeToolBarIcon(self.action)
        self.iface.removePluginMenu("", self.action)
        self.iface.removeToolBarIcon(self.action_cfg)
        # Disconnect controller signals
        self.controller.disconnect()
        # Close and remove dock, disconnect widget
        self.iface.removeDockWidget(self.ctrl_dock)
        self.iface.removeDockWidget(self.plot_dock)
        self.ctrl_widget.disconnect()
        self.plot_widget.disconnect()
        self.ctrl_dock.close()
        self.plot_dock.close()
