.. image:: graphics/Northwestern_purple_RGB.png
   :width: 30%
   :align: right
   :alt: nulogo


=============
Example Usage
=============

**FlyEye Analysis** provides a wide range of functionality for analyzing the expression data contained in ``.silhouette`` files. A very brief survey of some basic operations is provided below. For detailed usage instructions please see the :ref:`documentation <documentation>`.


Loading Silhouette Data
-----------------------

See :ref:`Getting Started <start>` for input data requirements.


To load an individual eye disc:

.. code-block:: python

   from flyeye.data import discs

   disc = discs.Disc.from_silhouette(path_to_silhouette_file)


To load an entire experiment:

.. code-block:: python

   from flyeye.data import experiments

   experiment = experiments.Experiment(path_to_experiment)




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

   reporter_channel = 'green'

   cells.plot_dynamics(reporter_channel)


Plot expression heterogeneity dynamics:

.. code-block:: python

   reporter_channel = 'green'

   fluctuations_channel = reporter_channel + '_flux'

   cells.plot_dynamics(fluctuations_channel)


Additional Examples
-------------------

For detailed usage examples, please refer to the `code <https://github.com/sebastianbernasek/pnt_yan_ratio>`_ used to generate the figures in our manuscript.
