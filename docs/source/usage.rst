=============
Example Usage
=============

**NU FlyEye:** *Analysis* provides a wide range of functionality for analyzing expression data measured using **NU FlyEye:** *Silhouette*. A very brief introduction to some core functionalities is provided here. For detailed usage instructions please see the :ref:`documentation <documentation>`.


Loading an experiment
---------------------

Import an experiment from a directory containing ``.silhouette`` files:


.. code-block:: python

   from flyeye.data import experiments

   path = './data'
   experiment = experiments.Experiment(path)


Selecting a disc
----------------

Select a specific disc:

.. code-block:: python

   disc_id = 2
   disc = experiment.discs[disc_id]



Querying by cell type
---------------------

Select a specific cell type:

.. code-block:: python

   cell_type = 'pre'
   cells = disc.select_cell_type(cell_type)


Timeseries Visualization
------------------------

Plot expression dynamics:

   fluorescence_channel = 'green'
   cells.plot_dynamics(fluorescence_channel)


Additional Examples
-------------------

For detailed usage examples, please refer to the `code <https://github.com/sebastianbernasek/pnt_yan_ratio>`_ used to generate the figures in our manuscript.
