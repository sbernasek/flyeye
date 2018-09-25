.. image:: graphics/Northwestern_purple_RGB.png
   :width: 30%
   :align: right
   :alt: nulogo


Getting Started
===============

The fastest way to gain familiarity with *FlyEye Analysis* is to start with a working example. Please feel free to use the data from `our study <https://github.com/sebastianbernasek/pnt_yan_ratio>`_ of Pnt and Yan expression during eye development.

To analyze your own segmented microscopy data, please read on!


Input File Structure
--------------------

*FlyEye Analysis* supports analysis of expression data contained in ``.silhouette`` files. Each of these files corresponds to a single *Drosophila* eye disc that has been marked with fluorescent reporters, dissected, imaged, segmented, and annotated. The ``.silhouette`` filetype follows a standardized structure:

.. code-block:: bash

   example.silhouette
   ├── feud.json        # stack dimensionality
   ├── feed.json        # measurement annotations
   ├── 0.json           # first layer measurements
   ├── 0.png            # first layer image
   |
   | ...
   |
   ├── N.json           # Nth layer measurements
   └── N.png            # Nth layer image

The layer images contained in a ``.silhouette`` file are compressed versions of the original full resolution microscopy. They provide a clear visual impression of the original images while maintaining a manageable filesize. However, they are not suitable as a basis for further expression quantification.


Data Requirements
-----------------

*NU FlyEye Analysis* requires:

 - All measurements of interest must be labeled. Unlabeled measurements are ignored upon import.

 - Progenitors must be labeled 'p' or 'pre'. Other names are possible but would require manual modification of the :ref:`flyeye.data <data>` source code.

 - R8 cells must be fully annotated within a contiguous region. The automated conversion of measurements to developmental timepoints is dependent upon regularly spaced R8 measurements.
