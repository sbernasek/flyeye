.. image:: graphics/Northwestern_purple_RGB.png
   :width: 30%
   :align: right
   :alt: nulogo

.. _start:


Getting Started
===============

The fastest way to gain familiarity with **FlyEye Analysis** is to start with a working example. Please feel free to use the data from `our study <https://github.com/sebastianbernasek/pnt_yan_ratio>`_ of Pnt and Yan expression during eye development.

We recommend reading the sections below before working with your own data.



Data Format
-----------

**FlyEye Analysis** supports analysis of expression data contained in ``.silhouette`` files. Each of these files corresponds to a single eye disc that has been marked with fluorescent reporters, dissected, imaged, segmented, and annotated. See the **FlyEye Silhouette** documentation for tips on generating ``.silhouette`` files from your microscopy data.


The ``.silhouette`` filetype includes image, measurement (reporter levels), and annotation data for each layer in a 3D image stack. The contents are structured as:

.. code-block:: bash

   example.silhouette
   ├── feud.json        # stack dimensionality
   ├── feed.json        # measurement annotations
   ├── 0.json           # first layer measurements
   ├── 0.png            # first layer image
   |
   ├── ... N.json       # Nth layer reporter levels (measurements)
   └── ... N.png        # Nth layer image

Each ``<layer_number>.json`` file contains all reporter levels measured during segmentation of the corresponding layer. Measured reporter levels reflect the mean pixel intensity within each nuclear contour, evaluated across all reporter wavelengths. These values are raw measurements; all subsequent normalization and processing by **FlyEye Analysis** are performed in local memory.

The ``feud.json`` file contains all user-assigned contour labels. **FlyEye Analysis*** automatically pairs measurements with their corresponding labels upon import of a ``.sihouette`` file. Unlabeled contours are ignored.

The ``<layer_number>.png`` images are compressed versions of the original microscopy. They provide a clear visual representation of the original images, but they are not suitable for expression quantification.


.. Note::
   The initial release of the **NU FlyEye** platform only supports RGB image stacks. The available reporter colors are thus limited to red, green, and blue. One of these reporter colors must be reserved for a nuclear marker in order to facilitate segmentation via **FlyEye Silhouette**. This leaves at most two reporter colors available for measuring target gene expression in any one experiment.


Data Requirements
-----------------

To analyze a ``.silhouette`` file:

 - R8 cells must be fully annotated within a locally contiguous region. [*]_

 - Only one measurement should be labeled per cell that appears in the 3-D image stack.

 - Progenitors must be labeled 'p' or 'pre'.

 - R8 cells must be labeled 'r8' or 'R8'.


.. [*] Timeseries construction relies upon regularly spaced R8 measurements. This requirement may be relaxed if estimated developmental times are ignored.


.. Note::
   Custom labels for cell types other than progenitors and R8 cells are possible without any modification of the :ref:`flyeye.data <data>` source code.



Data Management
---------------

**FlyEye Analysis** offers three levels of organization for managing cell measurement data. At the highest level, measurements are combined between eye discs collected under similar experimental conditions. We recommend organizing your ``.silhouette`` files in an equivalent manner by creating a separate directory for each experiment:

.. code-block:: bash

   data
   ├── experiment_A
   |   ├── eye0.silhouette
   |   |
   |   └── ... eyeN.silhouette
   |
   └── ... experiment_Z


Loading Data
------------

Measurement data must be loaded as ``data.discs.Disc`` instances prior to analysis. Several important operations are automatically triggered upon instantiation of a ``Disc``:

#. Each cell is assigned a developmental age based on its proximity to the furrow
#. Expression levels are normalized against the level of the reporter used to mark cell nuclei
#. The expression ratio between the two remaining reporters is evaluated
#. Moving average expression trends are evaluated for each labeled cell type

These operations are governed by a handful of user-specified parameters such as furrow velocity and the reporter color used to mark cell nuclei. These parameters must be specified in accordance with your particular dataset.

.. Note::
   **FlyEye Analysis** assumes that one of the three available reporter colors was reserved for a nuclear marker. The expression ratio assigned to each cell is evaluated using the two remaining reporter colors.


To load an individual ``.silhouette`` file:

.. code-block:: python

   from flyeye.data import discs

   path_to_disc = './data/experiment_A/eye0.silhouette'

   disc = discs.Disc.from_silhouette(path_to_disc)


Alternatively, the ``experiments.Experiment`` constructor will automatically load and combine all discs within a specified directory:

.. code-block:: python

   from flyeye.data import experiments

   path_to_experiment = './data/experiment_A'

   experiment = experiments.Experiment(path_to_experiment)


**Your data are now ready for analysis!**
