=============
Example Usage
=============

**NU FlyEye: Analysis** provides a wide range of functionality for analyzing the expression data contained in ``.silhouette`` files. A very brief survey of some basic operations is provided below. For detailed usage instructions please see the :ref:`documentation <documentation>`.


Loading an experiment
---------------------

Import an experiment from a directory containing ``.silhouette`` files:


.. code-block:: python

   from flyeye.data import experiments

   path = './data'
   experiment = experiments.Experiment(path)


Querying measurements
---------------------

Select by disc ID:

.. code-block:: python

   disc_id = 2
   disc = experiment.discs[disc_id]

Select by cell type:

.. code-block:: python

   cell_type = 'pre'
   cells = disc.select_cell_type(cell_type)

Select by developmental time:

.. code-block:: python

   tmin = 5
   tmax = 15
   cells = disc.select_by_position(tmin=tmin, tmax=tmax)


Visualizing Dynamics
--------------------

Plot expression dynamics:

.. code-block:: python

   fluorescence_channel = 'green'
   cells.plot_dynamics(fluorescence_channel)


Additional Examples
-------------------

For detailed usage examples, please refer to the `code <https://github.com/sebastianbernasek/pnt_yan_ratio>`_ used to generate the figures in our manuscript.
