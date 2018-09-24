__author__ = 'Sebastian Bernasek'

from os.path import join
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec

from .discs import Disc
from .cells import Cells
from ..processing.triangulation import Triangulation
from ..processing.alignment import DiscAlignment, ExperimentAlignment
from ..analysis.correlation import SpatialCorrelation


class Experiment:
    """
    Object representing multiple eye discs obtained under a single set of conditions.

    Attributes:

        discs (dict) - {disc ID: data.discs.Disc} pairs

        num_discs (int) - number of discs within experiment

    """

    def __init__(self, dirpath,
                 normalization='red',
                 auto_alignment=True,
                 **kwargs):
        """
        Instantiate object representing all discs obtained under a single set of conditions.

        Args:

            dirpath (str) - path to directory containing silhouette files

            normalization (str) - normalization channel

            auto_alignment (bool) - if True, align discs

            kwargs: keyword arguments for disc instantiation

        """

        self.discs = self.load(dirpath, normalization=normalization, **kwargs)

        # align discs
        if auto_alignment:
            self.align_discs()
            self.align_to_first_r8()

        # count discs
        self.num_discs = len(self.discs)

    def __getitem__(self, key):
        """ Return disc. """
        return self.discs[key]

    @staticmethod
    def load(dirpath, normalization='red', **kwargs):
        """
        Load discs from silhouette files.

        Args:

            dirpath (str) - path to directory containing silhouette files

            normalization (str) - normalization channel

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

    def align_discs(self, channel='green'):
        """
        Align all discs within experiment.

        Args:

            channel (str) - expression channel by which discs are aligned

        """
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
        t = sorted(reference.select_cell_type('r8').df.t.values)[1]

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
            disc.df['disc_id'] = disc_id

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

            df (DataFrame) - cells concurrent with reference cell type

        """

        # aggregate cells from just before/after their identification
        progenitors = Cells()
        references = Cells()

        for disc in self.discs.values():

            # select reference cells
            ref = disc.select_cell_type(reference_types)
            n_current = len(ref.df)
            if n_current == 0:
                continue

            # get time of first reference cell
            tmin = ref.df.iloc[upper_slip]['t'] - lower_slip

            # get time of Nth (or last) reference cell
            if n_current >= N:
                tmax = ref.df.iloc[N-1]['t']
            else:
                tmax = ref.df.iloc[-1]['t']

            # select concurrent progenitors and reference cells
            pre = disc.select_cell_type('pre')
            pre = pre.select_by_position(tmin=tmin, tmax=tmax)
            ref = ref.select_by_position(tmin=tmin, tmax=tmax)

            # append cell selections
            progenitors += pre
            references += ref

        # label precursors as multipotent
        progenitors.df['Population'] = 'Multipotent'
        progenitors.df['original_idx'] = progenitors.df.index

        # label neurons as differentiated
        references.df['Population'] = 'Differentiated'
        references.df['original_idx'] = references.df.index

        # label with corresponding reference cell type and append to data
        df = pd.concat((progenitors.df, references.df))
        df['ReferenceType'] = '/'.join([n.upper() for n in reference_types])

        return df

    def get_early_neuron_df(self,
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

            df (DataFrame) - measurement data for early R cells and concurrent progenitors

        """

        cell_types = [['r8'], ['r2', 'r5'], ['r3', 'r4'], ['r1', 'r6'], ['r7']]

        df = pd.DataFrame()
        for types in cell_types:
            data = self.select_by_concurrency(types, N, lower_slip, upper_slip)
            df = pd.concat([df, data])
        return df

    def get_spatial_correlations(self,
                                 cell_type='pre',
                                 channel='green',
                                 y_only=False,
                                 raw=False,
                                 discs_included='all',
                                 **selection_kw):
        """
        Compile SpatialCorrelation instance for all specified cells.

        Args:

            cell_type (str) - type of cells to select

            channel (str) - expression channel for which correlations are desired

            y_only (bool) - if True, only use y-component of data

            raw (bool) - if True, removes normalization

            discs_included (list or str) - included discs, defaults to all

            selection_kw: keyword arguments for cell position selection

        Returns:

           corr (analysis.correlation.SpatialCorrelation)

        """

        # instantiate empty SpatialCorrelation object
        corr = SpatialCorrelation()

        # specify included discs
        if discs_included == 'all':
            discs = self.discs.values()
        else:
            discs = [self.discs[i] for i in discs_included]

        # iterate across all included discs
        for i, disc in enumerate(discs):

            # select cells of specified type within specified time window
            cells = disc.select_cell_type(cell_types=cell_type)
            cells = cells.select_by_position(**selection_kw)

            # concatenate fluctuations to existing correlation object
            corr += SpatialCorrelation(cells, channel, y_only=y_only, raw=raw)

        return corr


