from unittest import TestCase
from os.path import join
from pandas import DataFrame
from flyeye.data import Experiment, Disc, ImageStack
from flyeye.data.cells import Cells


class TestDisc(TestCase):
    """
    Tests for Disc class.
    """

    def setUp(self):
        """ Initialize test instance of Disc. """
        self.dpath = 'flyeye/tests/fixtures/disc.silhouette'
        self.disc = Disc.from_silhouette(self.dpath, 'ch0', recompile=False)

    def test_load_imagestack(self):
        """ Load imagestack.  """
        stack = self.disc.load_imagestack()
        self.assertTrue(isinstance(stack, ImageStack))

    def test_select_celltype(self):
        """ Select single cell type.  """
        cells = self.disc.select_cell_type('pre')
        self.assertTrue(isinstance(cells, Cells))

    def test_select_multiple_celltypes(self):
        """ Select multiple cell types.  """
        cells = self.disc.select_cell_type(['r8', 'r1', 'r2'])
        self.assertTrue(isinstance(cells, Cells))


class TestExperiment(TestCase):
    """
    Tests for Experiment class.
    """

    def setUp(self):
        """ Initialize test instance of Disc. """
        self.epath = 'flyeye/tests/fixtures'
        self.exp = Experiment(self.epath, 0, auto_alignment=False)

    def test_get_disc(self):
        """ Get disc instance..  """
        disc = self.exp[0]
        self.assertTrue(isinstance(disc, Disc))

    def test_select_celltype(self):
        """ Select single cell type. """
        cells = self.exp.get_cells('pre')
        self.assertTrue(isinstance(cells, Cells))

    def test_select_multiple_celltypes(self):
        """ Select multiple cell types. """
        cells = self.exp.get_cells(['pre', 'r8'])
        self.assertTrue(isinstance(cells, Cells))

    def test_get_early_neuron_data(self):
        """ Get early neuron dataframe. """
        data = self.exp.get_early_neuron_data()
        self.assertTrue(isinstance(data, DataFrame))

