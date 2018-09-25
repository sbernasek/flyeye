===================
NU FlyEye: Analysis
===================


**NU FlyEye: Analysis** is part of the **NU FlyEye** platform for studying gene expression in the developing *Drosophila* eye. The analysis package provides methods for analyzing expression data measured *in vivo* using **NU FlyEye: Silhouette**.

Given one or more ``.silhouette`` files, **NU FlyEye: Analysis** facilitates:

   - **Timeseries Construction:** convert measurements to timepoints
   - **Data Management:** query cells by developmental age and cell type
   - **Dynamic Analysis:** analyze and visualize expression dynamics
   - **Spatial Analysis:** detect and quantify spatial expression patterns

The initial release is limited to basic tools required to replicate [our study](https://github.com/sebastianbernasek/pnt_yan_ratio) of Pnt and Yan expression during retinal patterning.



Installation
============

After downloading the [latest distribution](https://github.com/sebastianbernasek/flyeye/archive/v0.1.0-beta.tar.gz), the simplest method is to install via ``pip``:

    pip install flyeye-0.1.0-beta.tar.gz



Using NU FlyEye: Analysis
=========================

See our [documentation](https://sebastianbernasek.github.io/flyeye/index.html) page.



Examples
========

For an example of a complete project utilizing the entire **NU FlyEye** platform, please refer to the [code](https://github.com/sebastianbernasek/pnt_yan_ratio) used to generate the figures in our manuscript.
