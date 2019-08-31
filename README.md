# TDALabs

This series of Jupyter Notebooks serves as a walkthrough of topological data analysis, topological time series analysis (including rhythm analysis in music and periodicity / quasiperiodicity quantification in video), and 3D shape analysis.  The main engine behind all of the code is a Python port of the [ripser](http://github.com/ctralie/ripser) library.

This started off as a tutorial for the summer workshop "[Mathematical Methods for High-Dimensional Data Analysis](http://www-m15.ma.tum.de/Allgemeines/SummerSchool2016)."  Now it is used more generally to support pedagogical activities to support the NSF big data grant DKA-1447491, as well as assisting with the [ICERM Summer Undergraduate Program](https://icerm.brown.edu/summerug/2017/) and the workshop "<a href = "http://www.science.unitn.it/cirm/TDAPH2018.html">Topological Data Analysis and Persistent Homology</a>" in Levico Terme (Trento), Italy.

## Suggested Order
For 4 sessions with 90 minutes each (as in Levico), the following order is suggested (with some wiggle room between sessions)

### Session 1
* Basic Shapes
* Wasserstein And Bottleneck

### Session 2
* SlidingWindow1-Basics
* SlidingWindow2-PersistentHomology
* SlidingWindow3-AudioApplications
* SlidingWindow4-Video

### Session 3
* Approximate Sparse Filtrations
* Image Patches
* Diffusion Maps And TDA

### Session 4
* 3D Shapes
* Whatever is leftover


# Installation Instructions

Below are instructions for installing and running these tutorials

## Checking out code

~~~~~ bash
git clone --recursive https://github.com/ctralie/TDALabs.git
~~~~~

## Installing Jupyter Notebook And Other Dependencies

To run these modules, you will need to have jupyter notebook installed with a *Python 3* backend with numpy, scipy, and matplotlib.  The easiest way to install this is with Anaconda:

<a href = "https://www.anaconda.com/download/">https://www.anaconda.com/download/</a>

Once you have downloaded and installed all of these packages, navigate to the root of this repository and type the following commands, which will install dependencies

~~~~~ bash
pip install cython
pip install ripser
pip install imageio
~~~~~

## Updating code / submodules

At times, some updates may have happened to submodule dependencies (pykhs for the heat kernel signature and dreimac for projective coordinates).  To update these, type

~~~~~ bash
git submodule update --init
git submodule update --remote
~~~~~

## Installing avconv
For loading video in SlidingWindow4-Video, you will need to install the [ffmpeg](https://ffmpeg.org/download.html) binaries

## Running the code

At the root of this directory, type

~~~~~ bash
jupyter notebook
~~~~~

This will launch a browser window, where you can run the modules.  Click on one of the files (e.g. Basic Shapes.ipynb) to begin
