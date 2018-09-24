.. _data:

DATA
====

**NU FlyEye: Analysis** is built upon quantitative expression data measured using **NU FlyEye: Silhouette**. These measurements are obtained from ``.silhouette`` files using the ``flyeye.data.silhouette`` submodule. The data may then be queried at three levels:

  1. ``Cells`` objects contain one or more expression level measurements

  2. ``Disc`` objects contain all expression level measurements from a single ``.silhouette`` file

  3. ``Experiment`` objects contain multiple ``Disc`` instances collected under similar conditions

The images contained in ``.silhouette`` files are managed separarely by the ``flyeye.data.image`` submodule.


Silhouette
----------

Interface for managing the ``.silhouette`` filetype.

.. automodule:: flyeye.data.silhouette
   :members:


Cells
-----

Cells are a collection of one or more labeled expression measurements. Cells may be a subsample of a single Disc, or a combination of several discs.

.. automodule:: flyeye.data.cells
   :members:


Discs
-----

Discs are the set of cells that comprise an individual eye disc. This submodule contains all of the methods required to convert an individual disc's cell measurements into developmental timepoints. Measurements are imported from ``.silhouette`` files at this level.

.. automodule:: flyeye.data.discs
   :members:


Experiments
-----------

An Experiment is a set of one or more Disc objects collected under the same experimental conditions. Measurements are often aggregated between discs and compared at this level.

.. automodule:: flyeye.data.experiments
   :members:


Image
------

A submodule for managing, processing, and visualizing images contained in a ``.silhouette`` file.

.. automodule:: flyeye.data.image
   :members:

