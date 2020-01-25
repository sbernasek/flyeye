from unittest import TestCase
from os.path import join
from flyeye.data import Experiment, Disc


class TestIO(TestCase):
    """
    Tests for FlyEye/Silhouette interface.
    """

    def setUp(self):
        """ Initialize paths to test data. """
        self.epath = 'flyeye/tests/fixtures'
        self.dpath = join(self.epath, 'disc.silhouette')

    def test_load_disc(self):
        """ Load individual disc from Silhouette file.  """
        disc = Disc.from_silhouette(self.dpath, 'ch0', recompile=False)
        self.assertTrue(isinstance(disc, Disc))

    def test_load_disc_integer_normalization(self):
        """ Load individual disc from Silhouette file, using integer.  """
        disc = Disc.from_silhouette(self.dpath, 0, recompile=False)
        self.assertTrue(isinstance(disc, Disc))

    def test_load_disc_recompile(self):
        """ Load individual disc from Silhouette file, recompiling data.  """
        disc = Disc.from_silhouette(self.dpath, 'ch0', recompile=True)
        self.assertTrue(isinstance(disc, Disc))

    def test_load_experiment(self):
        """ Load experiment from collection of Silhouette files.  """
        exp = Experiment(self.epath,
                         normalization='ch0',
                         auto_alignment=False)
        self.assertTrue(isinstance(exp, Experiment))
