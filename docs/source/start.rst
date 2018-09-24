Getting Started
===============

The fastest way to gain familiarity with **NU FlyEye:** *Analysis* is to start with a working example. Please feel free to use the data from `our study <https://github.com/sebastianbernasek/pnt_yan_ratio>`_ of Pnt and Yan expression during eye development.

To analyze your own segmented microscopy data, please read on!


Input File Structure
--------------------

**NU FlyEye:** *Analysis* is only compatible with ``.silhouette`` format expression data measured via **NU FlyEye:** *Silhouette*. Each of these files corresponds to a single *Drosophila* eye disc that has been marked with fluorescent reporters, dissected, imaged, segmented, and annotated. The ``.silhouette`` filetype follows a standardized structure:

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

**NU FlyEye:** *Analysis* supports analysis of individual ``.silhouette`` files, as well as collections of eye discs collected under similar experimental conditions.

** NOTE ** that the layer images contained in a ``.silhouette`` file are compressed versions of the full resolution microscopy data, and are not suitable as a basis for repeated expression quantification.


Data Requirements
-----------------

**NU FlyEye:** *Analysis* requires that cell type annotation be complete, as unlabeled measurements are ignored upon import.
