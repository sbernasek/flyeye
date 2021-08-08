from os.path import join
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec

from .discs import Disc
from .cells import Cells

from ..utilities.iteration import Iterator
from ..utilities.string_handling import format_channel
from ..processing.triangulation import Triangulation
from ..processing.alignment import DiscAlignment, ExperimentAlignment


class Experiment:
    """
    Object representing multiple eye discs obtained under a single set of conditions.

    Attributes:

        discs (dict) - {disc ID: data.discs.Disc} pairs

        num_discs (int) - number of discs within experiment

    """

    def __init__(self, dirpath, normalization,
                 auto_alignment=True,
                 align_by='ch1_normalized',
                 **kwargs):
        """
        Instantiate object representing all discs obtained under a single set of conditions.

        Args:

            dirpath (str) - path to directory containing silhouette files

            normalization (str or int) - normalization channel

            auto_alignment (bool) - if True, align discs

            align_by (str or int) - channel used to align discs

            kwargs: keyword arguments for disc instantiation

        """
        self.discs = self.load(dirpath, normalization=normalization, **kwargs)

        # align discs
        if auto_alignment:
            self.align_discs(align_by)
            self.align_to_first_r8()

    def __getitem__(self, idx):
        """ Returns disc indexed by <idx>. """
        return self.discs[idx]

    def __iter__(self):
        """ Iterate over discs. """
        return Iterator(list(self.discs.values()))

    @property
    def num_discs(self):
        """ Number of discs in experiment. """
        return len(self.discs)

    @property
    def num_progenitors(self):
        """ Number of progenitor measurements in experiment. """
        return len(self.get_cells('pre').data)

    @staticmethod
    def load(dirpath, normalization, **kwargs):
        """
        Load discs from silhouette files.

        Args:

            dirpath (str) - path to directory containing silhouette files

            normalization (str or int) - normalization channel

            kwargs: keyword arguments for disc instantiation

        Returns:

            discs (dict) - {disc_id: data.discs.Disc} pairs

        """

        # identify silhouette files
        silhouette_paths = sorted(glob(join(dirpath, '*.silhouette')))

        # load discs
        discs = {}
        for i, path in enumerate(silhouette_paths):
            discs[i] = Disc.from_silhouette(path,
                                            normalization=normalization,
                                            **kwargs)

        return discs

    def set_ratio(self, num, den):
        """
        Add fluorescence ratio to each disc's dataframe, defined by <num>/<den> channels.
        """
        for disc in self.discs.values():
            disc.set_ratio(num, den)

    def align_discs(self, channel):
        """
        Align all discs within experiment.

        Args:

            channel (str) - expression channel by which discs are aligned

        """
        channel = format_channel(channel)
        al = ExperimentAlignment(self, channel=channel)
        self.discs = al.get_aligned_experiment().discs

    def get_pairwise_alignment(self, window_size=10, **kw):
        """
        Compute pairwise quality of alignment between each disc.

        Args:

            window_size (int) - number of cells for smoothing

            kw: keyword arguments for DiscAlignment

        Returns:

            scores (np.ndarray) - mean quality of alignment for each disc

        """

        # compute pairwise alignment between discs
        N = self.num_discs
        scores = np.zeros((N, N))
        for i, d0 in self.discs.items():
            for j, d1 in self.discs.items():
                al = DiscAlignment(d0, d1, window_size=window_size, **kw)
                scores[i, j] = al.score

        # mask diagonal
        mask = np.ones(scores.shape, dtype=bool)
        np.fill_diagonal(mask, 0)

        return scores[mask].reshape(N, N-1).mean(axis=1)

    def apply_lag(self, lag=0):
        """
        Apply time shift to all discs in experiment.

        Args:

            lag (float) - time shift applied to each disc

        """
        _ = [disc.apply_lag(offset=lag) for disc in self.discs.values()]

    def align_to_first_r8(self, disc_id=0):
        """
        Shift all discs s.t. t=0 is the first R8 in the reference disc.

        Args:

            disc_id (int) - index of disc used as reference

        """

        # get time of first R8
        reference = self.discs[disc_id]
        t = sorted(reference.select_cell_type('r8').data.t.values)[1]

        # apply lag
        self.apply_lag(lag=-t)

    def get_cells(self, cell_type='pre', **selection_kw):
        """
        Return Cells object for all specified cells.

        Args:

            cell_type (str or list) - type of cells to select

            selection_kw: keyword arguments for cell position selection

        Returns:

            cells (data.cells.Cells)

        """

        # assign disc_id
        for disc_id, disc in self.discs.items():
            disc.data['disc_id'] = disc_id

        # get all cells
        cells = np.sum(list(self.discs.values()))

        # filter cell selection
        cells = cells.select_cell_type(cell_type)
        cells = cells.select_by_position(**selection_kw)

        # sort inplace
        cells.sort(by='t')

        return cells

    def select_by_concurrency(self,
                              reference_types,
                              N=10,
                              lower_slip=0,
                              upper_slip=0):
        """
        Select cells concurrent with first N identified cells of reference cell type.

        Args:

            reference_types (array like) - reference cell type(s)

            N (int) - number of reference cells defining time window

            lower_slip (float) - extension before first reference cell, hours

            upper_slip (int) - reference cells skipped (excludes outliers)

        Returns:

            data (DataFrame) - cells concurrent with reference cell type

        """

        # aggregate cells from just before/after their identification
        progenitors = Cells()
        references = Cells()

        for disc_id, disc in self.discs.items():

            # select reference cells
            ref = disc.select_cell_type(reference_types)
            ref.data['disc_id'] = disc_id
            n_current = len(ref.data)
            if n_current == 0:
                continue

            # get time of first reference cell
            tmin = ref.data.iloc[upper_slip]['t'] - lower_slip

            # get time of Nth (or last) reference cell
            if n_current >= N:
                tmax = ref.data.iloc[N-1]['t']
            else:
                tmax = ref.data.iloc[-1]['t']

            # select concurrent progenitors and reference cells
            pre = disc.select_cell_type('pre')
            pre.data['disc_id'] = disc_id
            pre = pre.select_by_position(tmin=tmin, tmax=tmax)
            ref = ref.select_by_position(tmin=tmin, tmax=tmax)

            # append cell selections
            progenitors += pre
            references += ref

        # label precursors as multipotent
        progenitors.data['Population'] = 'Multipotent'
        progenitors.data['original_idx'] = progenitors.data.index

        # label neurons as differentiated
        references.data['Population'] = 'Differentiated'
        references.data['original_idx'] = references.data.index

        # label with corresponding reference cell type and append to data
        data = pd.concat((progenitors.data, references.data))
        data['ReferenceType'] = '/'.join([n.upper() for n in reference_types])

        return data

    def get_early_neuron_data(self,
                          N=10,
                          lower_slip=0,
                          upper_slip=1):
        """
        Compile Dataframe of early R cells and concurrent progenitors.

        Args:

            N (int) - number of reference cells defining time window

            lower_slip (float) - extension before first reference cell, hours

            upper_slip (int) - reference cells skipped (excludes outliers)

        Returns:

            data (DataFrame) - measurement data for early R cells and concurrent progenitors

        """

        cell_types = [['r8'], ['r2', 'r5'], ['r3', 'r4'], ['r1', 'r6'], ['r7']]

        data = pd.DataFrame()
        for types in cell_types:
            x = self.select_by_concurrency(types, N, lower_slip, upper_slip)
            data = pd.concat([data, x])
        return data
