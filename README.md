
FlyEye Overview
===========

FlyEye provides methods for analyzing and visualizing expression data measured *in vivo* using the **FlyEye Silhouette** package for OSX. The initial release is limited to basic tools required to replicate [our study](https://github.com/sebastianbernasek/pnt_yan_ratio) of Pnt and Yan expression during retinal patterning in *Drosophila*.



Installation
=========

After downloading the [latest distribution](https://github.com/sebastianbernasek/flyeye/archive/v0.1.0-beta.tar.gz), the simplest method is to install via ``pip``:

    pip install flyeye-0.1.0-beta.tar.gz



FlyEye Input
=========

Package scope begins with one or more ``.silhouette`` files. Each of these files corresponds to a single *Drosophila* eye disc that has been marked with fluorescent reporters, dissected, and imaged. Following segmentation and annotation via **FlyEye Silhouette**, each ``.silhouette`` file contains quantitative expression level measurements for all cells annotated within the respective eye disc.

FlyEye requires that cell type annotation be complete, as unlabeled measurements are ignored upon import. FlyEye supports analysis of individual ``.silhouette`` files, or aggregation of multiple eye discs collected under similar experimental conditions.



FlyEye Modules
=========

FlyEye is organized into four modules:

* Data. Components for managing **FlyEye Silhouette** data. FlyEye provides three levels of organization:

  1. ``Cells`` objects contain one or more expression level measurements

  2. ``Disc`` objects contain all expression level measurements from a single ``.silhouette`` file

  3. ``Experiment`` objects contain multiple ``Disc`` instances collected under similar conditions

* Processing. Methods for converting cell measurements into expression time series.

* Dynamics. Methods for time series visualization.

* Analysis. Methods for quantitative analysis of expression data.



Example Usage
=========

Import an experiment from a directory containing ``.silhouette`` files:

    #!/usr/bin/env python

    from flyeye.data import experiments

    path = './silhouette_data'
    experiment = experiments.Experiment(path)


Select a specific disc:

    disc_id = 2
    disc = experiment.discs[disc_id]


Select a specific cell type:

    cell_type = 'pre'
    cells = disc.select_cell_type(cell_type)


Plot expression dynamics:

    fluorescence_channel = 'green'
    cells.plot_dynamics(fluorescence_channel)


Further Examples
-------------

For detailed usage examples, please refer to the [code](https://github.com/sebastianbernasek/pnt_yan_ratio) used to generate the figures in our manuscript.
