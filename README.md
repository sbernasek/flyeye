
===========
FlyEye
===========

FlyEye provides methods for analyzing and visualizing expression data measured *in vivo* using the **FlyEye Silhouette** package for OSX. The initial release is limited to basic tools required to replicate `our study of Pnt and Yan expression <https://github.com/sebastianbernasek/pnt_yan_ratio>`_ during retinal patterning in Drosophila.


FlyEye Installation
=========

After downloading the [latest distribution](/dist/FlyEye-0.1.0.tar.gz), the simplest method is to install via ``pip``::

    pip install FlyEye-0.1.0.tar.gz


FlyEye Input
=========

Package scope begins with one or more ``.silhouette`` files. Each of these files corresponds to a single *Drosophila* eye disc that has been marked with fluorescent reporters, dissected, and imaged. Following segmentation and annotation via **FlyEye Silhouette**, eacg ``.silhouette`` file contains quantitative expression level measurements for all cells annotated within the respective eye disc.

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

Importing an experiment::

    #!/usr/bin/env python

    from flyeye.data import experiments

    # define path to ``.silhouette`` file directory
    path = './silhouette_data'

    # load experiment
    experiment = experiments.Experiment(path)


Selecting a specific disc::

    # select specific disc
    disc_id = 2
    disc = experiment.discs[disc_id]


Selecting a specific cell type::

    # select specific cells
    cell_type = 'pre'
    cells = disc.select_cell_type(cell_type)


Visualizing expression dynamics::

    fluorescence_channel = 'green'
    cells.plot_dynamics(fluorescence_channel)


Further Examples
-------------

For detailed usage examples, please refer to the `code <https://github.com/sebastianbernasek/pnt_yan_ratio>`_ used to generate the figures in our manuscript.
