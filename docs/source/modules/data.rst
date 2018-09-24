DATA
====

**NU FlyEye:** *Analysis* provides three levels of organization for managing cell measurement data:

  1. ``Cells`` objects contain one or more expression level measurements

  2. ``Disc`` objects contain all expression level measurements from a single ``.silhouette`` file

  3. ``Experiment`` objects contain multiple ``Disc`` instances collected under similar conditions


Cells
------

Cells are a collection of one or more expression measurements.

.. automodule:: flyeye.data.cells
   :members:


Disc
------

A Disc is a collection of all labeled expression measurements within an individual eye disc.

.. automodule:: flyeye.data.discs
   :members:


Experiments
-----------

An Experiment is a set of one or more Disc objects collected under the same experimental conditions.

.. automodule:: flyeye.data.experiments
   :members:
