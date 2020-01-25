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

   >>> from flyeye.data import Disc
   >>> disc = Disc.from_silhouette(path_to_silhouette_file)


To load an entire experiment:

.. code-block:: python

   >>> from flyeye.data import Experiment
   >>> experiment = Experiment(path_to_experiment)



Querying measurements
---------------------

Select by disc ID:

.. code-block:: python

   >>> disc_id = 2
   >>> disc = experiment.discs[disc_id]

Select by cell type:

.. code-block:: python

   >>> cell_type = 'pre'
   >>> cells = disc.select_cell_type(cell_type)

Select by developmental time:

.. code-block:: python

   >>> tmin, tmax = 5, 15
   >>> cells = disc.select_by_position(tmin=tmin, tmax=tmax)


Visualizing Dynamics
--------------------

Plot expression dynamics:

.. code-block:: python

   >>> reporter_channel = 'ch0'
   >>> cells.plot_dynamics(reporter_channel)


Plot expression heterogeneity dynamics:

.. code-block:: python

   >>> reporter_channel = 'ch0'
   >>> fluctuations_channel = reporter_channel + '_flux'
   >>> cells.plot_dynamics(fluctuations_channel)


Additional Examples
-------------------

For additional usage examples, please refer to our `study <https://github.com/sebastianbernasek/pnt_yan_ratio>`_ of Pnt and Yan expression during photoreceptor specification.
