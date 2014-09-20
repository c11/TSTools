Time Series Tools (TSTools)
-------------------

### About
TSTools is a plugin for QGIS (version 2.0+) that helps visualize remote sensing time series by linking time series dataset models (objects that describe and characterize the time series) with user interface tools designed to harmonize the spatial and temporal dimensions of these large datasets.

While this QGIS plugin was originally designed for use with the "Continuous Change Detection and Classification" algorithm (Zhu and Woodcock 2014), I am working (slowly) to make the backend code and user interface extensible for use with any number of time series algorithms.

The goal is for users to describe their own data set structures, algorithm parameters, and algorithm outputs and then plug these customizations into the TSTools plugin by inheriting from the abstract base class "TimeSeries". This base class acts as an interface descriptor that characterizes what methods and properties are required for use within the user interface.

### Installation
#### Locally
This plugin has not been uploaded to the main QGIS plugin repository so installation will need to be done manually.

In most cases, the QGIS Python plugins folder will be located in your home directory within the ".qgis2/python/plugins" folder. Any plugins you have installed previously will be located here. For more information, see this excellent answer on [Stack Exchange](http://gis.stackexchange.com/questions/26979/how-to-install-a-qgis-plugin-when-offline).

1. Install QGIS and the required Python libraries (see requirements section below)
2. Download the file "tstools.zip" from this repository on Github.
3. Unzip the ZIP file to find the "tstools" folder.
4. Copy this "tstools" folder into your QGIS Python plugins directory (see above for where this is located)
5. Launch QGIS and open the Plugin Manage dialog (Plugins menu -> Manage and Install Plugins)
6. Check the box next to "TSTools" to enable the plugin

Two new icons will be added to the plugins toolbar. These icons have the letters "TS" in capital red colored letters. To initialize a timeseries dataset within the plugin, click the icon without the crosshair symbol. Point this dialog to your timeseries and configure any additional options before clicking "Okay". To retrieve the timeseries for any given pixel, add an image from your timeseries to QGIS using the "Images" tab and click the "TS" icon with the crosshairs to replace your current map tool with the "TSTools" map tool.

An example dataset for this plugin is located here:
https://github.com/ceholden/landsat_stack

#### Virtual machine demo
To help out people who find the installation of this software is not so straightforward (e.g., it is more difficult on Windows than Linux), I have created a virtual machine of the 14.04 LTS [Xubuntu distribution](http://xubuntu.org/) with everything installed. This virtual machine contains a full stack of softwares - GDAL, Python, QGIS, NumPy, SciPy, etc. - that are required to use the plugin. The virtual machine is formatted as a [VirtualBox image](https://www.virtualbox.org/) and I would recommend you to use [VirtualBox](https://www.virtualbox.org/) to run the virtual machine. VirtualBox is a free and open source softare that can create and host virtual machines and is comparable to commercial solutions such as VMWare or Parallels.

The virtual machine has been exported to a [VirtualBox appliance](http://www.virtualbox.org/manual/ch01.html#ovf) and uploaded to my university department's anonymous FTP server:

ftp://ftp-earth.bu.edu/ceholden/TSTools/

Please see the included README for further instructions. A md5sum of the virtual disk appliance is provided for confirming the file transfer integrity.

### Requirements
#### Main dependencies:

    Python (2.7.x tested)
    Numpy (1.7.x tested)
    GDAL (1.10.0 tested)
    
#### Additional dependencies:

- For reading CCDC results:
    Scipy (0.12.0 tested)
- For live plotting with YATSM (CCDC / BFAST inspired clone or hybrid):
    see https://github.com/ceholden/yatsm

#### Developer dependencies:
To help develop this plugin, you will need QGIS, Python, and the Qt developer tools for Python (for building). The Qt dependencies are available on Ubuntu in the "pyqt4-dev-tools" package.

### Features
##### Plot time series and time series model fits by clicking on image
<img src="https://raw.githubusercontent.com/ceholden/TSTools/master/docs/media/beetle_ts_2013.png" align="center" width=500/>

*Time series fit from Zhe Zhu's CCDC*

##### Plot features
+ Click a plot point and open corresponding image in QGIS
+ Adjust X and Y plot limits
+ Turn on or off model results
+ Export image as PNG, EPS, etc.

##### Quickly add/remove time series images from table
<img src="https://raw.githubusercontent.com/ceholden/TSTools/master/docs/media/tstools_imagetable.png" align="center" width=250/>

*Metadata columns coming soon*

##### Control image symbology for all time series images
<img src="https://raw.githubusercontent.com/ceholden/TSTools/master/docs/media/tstools_symbology.png" align="center" width=250/>

##### Add your own time series model with custom initialization requirements
<img src="https://raw.githubusercontent.com/ceholden/TSTools/master/docs/media/tstools_customconfig.png" align="center" width=250/>
