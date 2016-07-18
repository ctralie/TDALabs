# TUMTopoTimeSeries2016

For the summer workshop "Mathematical Methods for High-Dimensional Data Analysis"
http://www-m15.ma.tum.de/Allgemeines/SummerSchool2016

#Installation Instructions

Below are instructions for installing and running these tutorials

##Installing Jupyter Notebook

To run these modules, you will need to have jupyter notebook installed with a *Python 3* backend with numpy, scipy, and matplotlib.  The easiest way to install this is with Anaconda:

https://www.continuum.io/downloads

Once you have downloaded and installed all of these packages, navigate to the root of this repository and type

~~~~~ bash
jupyter notebook
~~~~~

This will launch a browser on your computer that will allow you to run the modules via the local Jupyter backend server

##Installing librosa
After you have the proper Python environment, you will need to install the [librosa library](https://github.com/librosa/librosa) for the third module on audio processing.

##Compiling the TDA Tools Library
To run code which computes TDA (modules 2-4), you will need to compile some C++ code, written by [Rann Bar-On](https://math.duke.edu/people/rann-bar).  This requires you to have CMake and OpenMP installed.  Change into the directory "h2phatclique" and issue the following commands

~~~~~ bash
mkdir build
cd build
cmake ..
make
cp ../bin/phatclique ../..
~~~~~


To test this, type
~~~~~ bash
python TDA.py
~~~~~

at the root of this repository.  This should compute a persistence diagram of a circle, which has only one persistence point.
