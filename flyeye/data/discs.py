__author__ = 'Sebastian Bernasek'

import numpy as np
import pandas as pd
from scipy.stats import norm

from .silhouette import Silhouette
from .image import ImageStack
from .cells import Cells
from ..dynamics.averages import detrend_signal
from ..processing.triangulation import Triangulation

# treat warnings as exceptions
#import warnings
#warnings.filterwarnings('error')


class Disc(Cells):
    """
    Object representing all cells in a single eye disc.

    Attributes:
    path (str) - unique silhouette filepath
    triangulation (processing.triangulation.Triangulation)

    Inherited Attributes:
    df (pd.DataFrame) - cell measurement data
    """

    def __init__(self,
                 df=None,
                 path=None,
                 normalization='red',
                 flip_about_yz=False,
                 flip_about_xy=False,
                 furrow_velocity=2.,
                 offset=None):
        """
        Instantiate object representing all cells in a single eye disc.

        Args:
        df (pd.DataFrame) - cell measurement data
        normalization (str) - normalization channel
        flip_about_yz (bool) - if True, invert about YZ plane
        flip_about_xy (bool) - if True, invert about XY plane
        path (str) - silhouette filepath
        furrow_velocity (float) - furrow inverse-velocity (hours per column)
        offset (float) - time by which disc is shifted from first R8
        """

        Cells.__init__(self, df, normalization=normalization)

        # set path to silhouette file
        self.path = path

        # standardize precursor labels
        self.standardize_labels()

        # orient disc (flip horizontally)
        if flip_about_yz is True:
            self.flip_about_yz()
            self.sort()

        # flip stack vertically
        if flip_about_xy is True:
            self.flip_about_xy()

        # normalize expression levels
        if normalization is not None:
            self.normalize_expression(base=255)
            self.normalize_by_reference(reference_channel=normalization)
            self.set_ratio()

        # construct triangulation of R8 positions
        self.triangulation = None
        r8_neurons = self.select_cell_type('r8')
        if len(r8_neurons.df) > 0:
            self.triangulation = Triangulation(self, threshold=1.75, furrow_velocity=furrow_velocity)

        # apply triangulation to compute estimated developmental times
        self.apply_time_scaling()

        # shift all cells such that first R8 occurs at time zero
        self.apply_lag(offset=offset)

        # detrend
        self.detrend()

    @staticmethod
    def from_silhouette(path,
                        normalization='red',
                        furrow_velocity=2.,
                        recompile=False,
                        **kwargs):
        """
        Instantiate disc from silhouette file.

        Args:
        path (str) - silhouette filepath
        normalization (str) - normalization channel
        furrow_velocity (float) - furrow inverse-velocity (hours per column)
        recompile (bool) - if True, recompile measurements from all layers
        kwargs: keyword arguments for disc instantiation

        Returns:
        disc (data.discs.Disc)
        """

        # load silhouette file
        silhouette = Silhouette(path, recompile=recompile)

        # instantiate disc
        disc = Disc(df=silhouette.df,
                    path=silhouette.path,
                    normalization=normalization,
                    flip_about_yz=silhouette.flip_about_yz,
                    flip_about_xy=silhouette.flip_about_xy,
                    furrow_velocity=furrow_velocity,
                    **kwargs)

        return disc

    def load_imagestack(self):
        """
        Load image stack from silhouette file.

        Returns:
        stack (data.image.ImageStack)
        """
        return ImageStack(self.path)

    def standardize_labels(self):
        """ Convert all alternate precursor labels to 'pre' """
        # for alternate in ['p', 'prepre', 'v']:
        #     self.df.ix[self.df.label==alternate, 'label'] = 'pre'
        names = ['p', 'prepre', 'v']
        self.df.loc[self.df.label.isin(names), 'label'] = 'pre'

    def flip_about_xy(self):
        """ Flip disc bottom to top. """
        self.df['centroid_y'] = self.df.centroid_y.max() - self.df.centroid_y + self.df.centroid_y.min()

    def flip_about_yz(self):
        """ Flip disc left to right. """
        self.df['centroid_x'] = self.df.centroid_x.max() - self.df.centroid_x + self.df.centroid_x.min()

    def apply_time_scaling(self):
        """ Apply distance-to-time scaling to generate time vector. """
        if self.triangulation is None:
            hours_per_pixel = 1
        else:
            hours_per_pixel = self.triangulation.hours_per_pixel
        self.hours_per_pixel = hours_per_pixel
        self.df['t'] = self.df['centroid_x'] * self.hours_per_pixel

    def apply_lag(self, offset=None):
        """
        Shift disc in time.

        Args:
        offset (str or float) - position of new origin, if 'first_r8' shift such that first R8 occurs at time zero
        """

        if offset is None:
            offset = 'first_r8'

        if offset == 'first_r8':
            r8_neurons = self.select_cell_type('r8')
            if len(r8_neurons.df) > 0:
                offset = -r8_neurons.df.centroid_x.min()
            else:
                offset = 0

        self.df['t'] += offset
        self.df['centroid_x'] += offset / self.hours_per_pixel

    def normalize_expression(self, base=255):
        """
        Convert each channel's fluorescence to a 0-1 scale.

        Args:
        base (float) - maximum fluorescence
        """
        for channel in ['red', 'green', 'blue']:
            self.df[channel] /= base
            self.df[channel+'_std'] /= base

    def normalize_by_reference(self, reference_channel='red'):
        """
        Normalize expression levels by the reference channel.

        Args:
        reference_channel (str) - channel used for normalization
        """
        normalization = self.df[reference_channel]
        for channel in ['red', 'green', 'blue']:
            if channel != reference_channel:
                self.df[channel] /= normalization

    def set_ratio(self):
        """ Add fluorescence ratio to dataframe. """

        if self.normalization == 'red':
            den = 'blue'
        else:
            den = 'red'

        if self.df[den].max() != 0:
            self.df['ratio'] = np.log2(self.df.green/self.df[den])
            self.df['raw_ratio'] = self.df.green/self.df[den]
        else:
            self.df['ratio'] = np.zeros(len(self.df))

    def detrend(self, order=1):
        """
        Add detrended fluctuations for each fluorescence channel within each cell type to disc dataframe.

        Args:
        order (int) - polyorder for local trendfitting
        """

        fluctuations = []
        for cell_type in self.df.label.unique():

            # select cells of specified type
            cells = self.select_cell_type(cell_type)

            data, columns = [], []
            for ch in ['red', 'green', 'blue', 'ratio']:

                # get samples for current cell type and channel
                x = cells.df[ch].values

                # determine window size for lowpass filter
                if cell_type == 'pre':
                    window_size = np.floor(x.size/10)
                else:
                    window_size = np.floor(x.size/5)

                # evaluate residuals by substracting lowpass filtered data
                if window_size > 3:
                    residuals, trend = detrend_signal(x, window_size, order)
                else:
                    trend = x.mean()
                    residuals = x - trend

                # normalize residuals
                if ch == 'ratio':
                    flux_raw = residuals
                else:
                    flux_raw = residuals/trend
                flux = np.abs(flux_raw)

                # append channel to list
                data.extend([trend, residuals, flux_raw, flux])
                columns.extend([ch+'_trend',
                                ch+'_residuals',
                                ch+'_flux_raw',
                                ch+'_flux'])

            # append cells fluctuations dataframe to list
            data = np.array(data).T
            data = data.reshape((len(cells.df.index), len(columns)))
            df = pd.DataFrame(data=data, index=cells.df.index, columns=columns)
            fluctuations.append(df)

        # join measurements with refisuals on measurement index
        self.df = self.df.join(pd.concat(fluctuations))

    def get_multipotent_layers(self, q=.9):
        """
        Determine which layers span multipotent cells. Bounds correspond to lower and upper quantiles taken from a normal fit to all progenitor layers.

        Args:
        q (float) - width parameter (interquantile range), 0 to 1

        Returns:
        bounds (np.ndarray, length 2) - lower and upper bounds
        """
        progenitors = self.select_cell_type('pre')
        layers = progenitors.df.layer.values
        bounds = norm(*norm.fit(layers)).ppf([(1-q)/2, (1+q)/2])
        return np.round(bounds).astype(int)
