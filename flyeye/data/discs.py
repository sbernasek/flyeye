__author__ = 'Sebastian Bernasek'

import warnings
import numpy as np
import pandas as pd
from scipy.stats import norm

from .silhouette import SilhouetteData
from .image import ImageStack
from .cells import Cells, format_channel
from ..dynamics.averages import detrend_signal
from ..processing.triangulation import Triangulation

# treat warnings as exceptions
#import warnings
#warnings.filterwarnings('error')


class DiscProperties:
    """ Properties for Disc class. """

    @property
    def path(self):
        """ Path to silhouette file. """
        return self.silhouette.path

    @property
    def flip_about_yz(self):
        """ If True, disc is inverted about YZ plane. """
        return self.silhouette.flip_about_yz

    @property
    def flip_about_xy(self):
        """ If True, disc is inverted about XY plane. """
        return self.silhouette.flip_about_xy


class Disc(Cells, DiscProperties):
    """
    Object representing all cells in a single eye disc.

    Attributes:

        silhouette (flyeye.SilhouetteData) - data from silhouette file

        furrow_velocity (float) - furrow inverse-velocity (hours per column)

        offset (float) - time by which disc is shifted from first R8

        bit_depth (int) - fluorescence intensity bit depth, log2 scaled

        triangulation (processing.triangulation.Triangulation)

    Inherited Attributes:

        normalization (str or int) - channel used to normalize intensities

        df (pd.DataFrame) - cell measurement data

    Properties:

        path (str) - unique silhouette filepath

        flip_about_yz (bool) - if True, disc is inverted about YZ plane

        flip_about_xy (bool) - if True, disc is inverted about XY plane

    """

    def __init__(self,
                 silhouette,
                 normalization=None,
                 furrow_velocity=2.,
                 offset=None,
                 bit_depth=None):
        """
        Instantiate object representing all cells in a single eye disc.

        Args:

            silhouette (flyeye.SilhouetteData) - data from silhouette file

            normalization (str or int) - channel used to normalize intensities

            furrow_velocity (float) - furrow inverse-velocity (hours per column)

            offset (float) - time by which disc is shifted from first R8

            bit_depth (int) - fluorescence intensity bit depth, log2 scaled

        """

        Cells.__init__(self, silhouette.df, normalization=normalization)

        # set sihouette attribute
        self.silhouette = silhouette

        # standardize precursor labels
        self.standardize_labels()

        # orient disc (flip horizontally)
        if self.flip_about_yz is True:
            self.flip_horizontal()
            self.sort()

        # flip stack vertically
        if self.flip_about_xy is True:
            self.flip_vertical()

        # normalize fluorescence intensities by bit depth
        if bit_depth is not None:
            self.bit_depth = bit_depth
            self.normalize_expression(base=2**self.bit_depth)

        # normalize fluorescence intensities by reference channel
        if self.normalization is not None:
            self.normalize_by_reference(self.normalization)

        # construct triangulation of R8 positions
        self.furrow_velocity = furrow_velocity
        self.triangulation = None
        r8_neurons = self.select_cell_type('r8')
        if len(r8_neurons.df) > 0:
            self.triangulation = Triangulation(self, threshold=1.75, furrow_velocity=furrow_velocity)

        # apply triangulation to compute estimated developmental times
        self.apply_time_scaling()

        # shift all cells such that first R8 occurs at time zero
        self.apply_lag(offset=offset)

        # detrend
        self.detrend(channels=self.normalized_channels)

    @staticmethod
    def from_silhouette(path,
                        normalization=None,
                        furrow_velocity=2.,
                        recompile=False,
                        **kwargs):
        """
        Instantiate disc from silhouette file.

        Args:

            path (str) - silhouette filepath

            normalization (str or int) - normalization channel

            furrow_velocity (float) - furrow inverse-velocity (hours per column)

            recompile (bool) - if True, recompile measurements from all layers

            kwargs: keyword arguments for disc instantiation

        Returns:

            disc (data.discs.Disc)

        """

        # load silhouette file data
        silhouette = SilhouetteData(path, recompile=recompile)

        # instantiate disc
        disc = Disc(silhouette,
                    normalization=normalization,
                    furrow_velocity=furrow_velocity,
                    **kwargs)

        return disc

    def load_imagestack(self):
        """
        Load image stack from silhouette file.

        Returns:

            stack (data.image.ImageStack)

        """
        return ImageStack.from_silhouette(self.path)

    def standardize_labels(self):
        """ Convert all alternate precursor labels to 'pre' """
        names = ['p', 'prepre', 'v']
        self.df.loc[self.df.label.isin(names), 'label'] = 'pre'

    def flip_vertical(self, save=False):
        """
        Flip disc bottom to top.

        Args:

            save (bool) - if True, save new orientation to silhouette file

        """
        ymax, ymin = self.df.centroid_y.max(), self.df.centroid_y.min()
        self.df['centroid_y'] = ymax - self.df.centroid_y + ymin

        # update feed file
        if save:
            self.silhouette.feed['orientation']['flip_about_xy'] = not self.flip_about_xy
            self.silhouette.write_feed()

    def flip_horizontal(self, save=False):
        """
        Flip disc left to right.

        Args:

            save (bool) - if True, save new orientation to silhouette file

        """
        xmax, xmin = self.df.centroid_x.max(), self.df.centroid_x.min()
        self.df['centroid_x'] = xmax - self.df.centroid_x + xmin
        if 't' in self.df.keys():
            tmax, tmin = self.df.t.max(), self.df.t.min()
            self.df['t'] = tmax - self.df.t + tmin

        # update feed file
        if save:
            self.silhouette.feed['orientation']['flip_about_yz'] = not self.flip_about_yz
            self.silhouette.write_feed()

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
                offset = -r8_neurons.df.t.min()
            else:
                offset = 0

        self.df['t'] += offset
        self.df['centroid_x'] += offset / self.hours_per_pixel

    def normalize_expression(self, max_value):
        """
        Convert each channel's fluorescence to a 0-1 scale.

        Args:

            max_value (float) - maximum fluorescence intensity

        """
        for channel in self.channels:
            self.df[channel] /= max_value
            self.df[channel+'_std'] /= max_value

    def normalize_by_reference(self, reference):
        """
        Normalize expression levels by the reference channel.

        Args:

            reference (str or int) - channel used for normalization

        """
        for channel in self.channels:
            if channel != reference:
                self.df[channel+'_norm'] = self.df[channel]/self.df[reference]

    def set_ratio(self, num, den):
        """
        Add fluorescence ratio to dataframe, defined by <num>/<den> channels.
        """

        num, den = format_channel(num), format_channel(den)

        if self.df[den].max() != 0:
            self.df['logratio'] = np.log2(self.df[num]/self.df[den])
            self.df['ratio'] = self.df[num]/self.df[den]
            self.detrend(channels=('logratio',))

        else:
            raise ValueError('Denominator channel is empty.')

    def detrend(self, channels, order=1):
        """
        Add detrended fluctuations for each fluorescence channel within each cell type to disc dataframe.

        Args:

            channels (iterable) - channels to be detrended

            order (int) - polyorder for local trendfitting

        """

        fluctuations = []
        for cell_type in self.df.label.unique():

            # select cells of specified type
            cells = self.select_cell_type(cell_type)

            data, columns = [], []
            for ch in channels:

                # skip missing channels
                if ch not in self.df.columns:
                    continue

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
                    trend = x.mean() * np.ones(x.size)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        residuals = x - trend

                # normalize residuals
                if ch == 'logratio':
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

        # join measurements with residuals on measurement index
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
