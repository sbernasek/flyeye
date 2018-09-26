FlyEye Analysis
===============


**FlyEye Analysis** is part of the **NU FlyEye** platform for studying gene expression in the developing *Drosophila* eye. The analysis package provides methods for analyzing expression data measured *in vivo* using **FlyEye Silhouette**.

The initial release is limited to basic tools required to replicate [our study](https://github.com/sebastianbernasek/pnt_yan_ratio) of Pnt and Yan expression during retinal patterning in *Drosophila*. Given one or more ``.silhouette`` files, analysis package facilitates:

   - **Timeseries Construction.** Convert measurements to timepoints

   - **Data Management.** Query cells by developmental age and cell type

   - **Dynamic Analysis.** Analyze and visualize expression dynamics

   - **Spatial Analysis.** Detect and quantify spatial expression patterns


Please refer to the [documentation](https://sebastianbernasek.github.io/flyeye/index.html#) page for tips on getting started with analyzing your data.


Dependencies
============

 - Python 3.6+

Required:

 - [Scipy](https://www.scipy.org/)
 - [Pandas](https://pandas.pydata.org/)
 - [PIL](https://pillow.readthedocs.io/en/5.2.x/)
 - [Matplotlib](https://matplotlib.org/)
 - [AstroML](https://pypi.org/project/astroML/)
 - [AstroML Addons](https://github.com/astroML/astroML_addons)


Installation
============

After downloading the [latest distribution](https://github.com/sebastianbernasek/flyeye/archive/v0.1.0-beta.tar.gz), the simplest method is to install via ``pip``:

    pip install flyeye-0.1.0-beta.tar.gz


Examples
========

For an example of a complete project utilizing the entire **NU FlyEye** platform, please refer to the [code](https://github.com/sebastianbernasek/pnt_yan_ratio) used to generate the figures in our manuscript.


Authors
=======

[Sebastian Bernasek](https://github.com/sebastianbernasek)
