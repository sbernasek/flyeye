from copy import deepcopy
from functools import reduce
from operator import add
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re

from ..utilities.string_handling import format_channel, standardize_channels
from ..dynamics.visualization import TimeseriesPlot, IntervalPlot
from ..dynamics.resampling import DiscResampler


class CellProperties:
    """ Properties of Cells object. """

    @property
    def channels(self):
        """ List of unique fluorescence channels. """
        return [s for s in self.data.columns if len(re.findall('ch[0-9]+$', s))]

    @property
    def normalized_channels(self):
        """ List of normalized channel names. """
        return [ch+'_normalized' for ch in self.channels]

    @property
    def num_channels(self):
        """ Number of unique fluorescence channels. """
        return len(self.channels)

    @property
    def is_sorted(self):
        return (self.data.centroid_x.diff() < 0).sum() == 0

    @property
    def cell_type_counts(self):
        return self.data.groupby('label')['label'].count()

    @property
    def cell_types(self):
        """ Unique cell types. """
        return self.cell_type_counts.index.tolist()


class Cells(CellProperties):
    """

    Object represents a population of cells. Each cell is completely described by a single record in a DataFrame of cell measurements. These measurements include cell positions, expression levels, and cell type annotations. Object may contain cells of one or more cell types.

    Attributes:

        data (pd.DataFrame) - cell measurement data

        normalization (str or int) - channel used to normalize intensities

    """

    def __init__(self, data=None, normalization=None):
        """
        Instantiate population of cells.

        Args:

            data (pd.DataFrame) - cell measurement data

            normalization (str or int) - channel used to normalize intensities

        """

        # store measurements
        if data is None:
            data = pd.DataFrame()
        self.data = standardize_channels(data)

        # store normalization
        self.normalization = format_channel(normalization)

        # standardize levels
        if len(self.data) > 0:
            self.sort()

    def __add__(self, cells):
        """ Concatenate second Cell instance. """
        cells = Cells(pd.concat((self.data, cells.data), sort=True), self.normalization)
        cells.sort(by='t')
        return cells

    def sort(self, by='centroid_x'):
        """
        Sort cell measurements in place.

        Args:

            by (str) - key on which measurements are sorted

        """
        self.data = self.data.sort_values(by=by, ascending=True)

    def apply_lag(self, lag):
        """
        Shift cells in time.

        Args:

            lag (float) - shift (NOTE: x-positions are unaffected)

        """
        self.data['t'] += lag

    def select_cell_type(self, cell_types):
        """
        Select subset of cells corresponding to a specified label.

        Args:

            cell_types (str or list) - type of cells to be selected (e.g. pre, r8)

        Returns:

            cells (data.cells.Cells)

        """

        # convert string to list
        if type(cell_types) == str:
            cell_types = [cell_types]

        # add both precursor labels if needed
        if 'pre' in cell_types:
            cell_types.append('p')
        elif 'p' in cell_types:
            cell_types.append('pre')

        # select cells
        data = self.data[self.data.label.apply(lambda x: x in cell_types)]

        # instantiate cells object
        cells = Cells(data, self.normalization)

        return cells

    def select_by_position(self,
                           xmin=-np.inf,
                           xmax=np.inf,
                           ymin=-np.inf,
                           ymax=np.inf,
                           zmin=-np.inf,
                           zmax=np.inf,
                           tmin=-np.inf,
                           tmax=np.inf):
        """
        Select subset of cells within specified spatial bounds.

        Args:

            xmin, xmax (float) - x-coordinate bounds

            ymin, ymax (float) - y-coordinate bounds

            zmin, zmax (float) - z-coordinate (layer number) bounds

            tmin, tmax (float) - time interval bounds

        Returns:

            cells (data.cells.Cells) - copied subset of cells

        """

        # initialize filter to include all cells
        data = deepcopy(self.data)

        # apply sequential filters
        data = data[data['centroid_x'].between(xmin, xmax)]
        data = data[data['centroid_y'].between(ymin, ymax)]
        data = data[data['layer'].between(zmin, zmax)]
        data = data[data['t'].between(tmin, tmax)]

        # instantiate subpopulation
        cells = Cells(data, self.normalization)

        return cells

    def get_nuclear_diameter(self):
        """
        Returns median nuclear diameter. Diameters are approximated as that of a circle with equivalent area to each nuclear contour.

        Returns:

            nuclear_diameter (float) - median diameter

        """
        return (2*np.sqrt(self.data.pixel_count/np.pi)).median()

    @staticmethod
    def get_binned_mean(x, values,
                        bins=None,
                        bin_width=1):
        """
        Bin cells and compute mean for each bin.

        Args:

            x (pd.Series) - coordinate on which to bin values

            values (pd.Series) - values to be aggregated

            bins (np array) - edges for specified bins

            bin_width (float) - width of bins used if no bins specified

        Returns:

            bin_centers (np.ndarray) - bin centers

            means (np array) - mean value within each bin

        """

        if bins is None:
            bins = np.arange(x.min(), x.max(), bin_width)
        bin_centers = [bins[i] + (bins[i+1] - bins[i])/2 for i in range(0, len(bins)-1)]
        means, _, _ = st.binned_statistic(x, values, statistic='mean', bins=bins)
        return bin_centers, means

    def plot_dynamics(self, channel,
                         ax=None,
                         scatter=False,
                         average=True,
                         interval=False,
                         marker_kw={},
                         line_kw={},
                         interval_kw={},
                         ma_kw={}):
        """
        Plot expression dynamics for specified channel.

        Args:

            channel (str) - expression channel

            ax (mpl.axes.AxesSubplot) - if None, create axes

            scatter (bool) - if True, add markers for each measurement

            average (bool) - if True, add moving average

            interval - if True, add confidence interval for moving average

            marker_kw (dict) - keyword arguments for marker formatting

            line_kw (dict) - keyword arguments for line formatting

            interval_kw (dict) - keyword arguments for interval formatting

            ma_kw (dict) - keyword arguments for interval construction

        Returns:

            ax (mpl.axes.AxesSubplot)

        """

        # sort values inplace
        self.sort('t')

        # instantiate TimeseriesPlot
        x, y = self.data.t.values, self.data[channel].values
        tsplot = TimeseriesPlot(x, y, ax=ax)

        # plot dynamics
        tsplot.plot(scatter=scatter,
                     average=average,
                     interval=interval,
                     marker_kw=marker_kw,
                     line_kw=line_kw,
                     interval_kw=interval_kw,
                     ma_kw=ma_kw)

        return tsplot.ax

    def plot_resampled_dynamics(self, channel,
                         ax=None,
                         average=True,
                         interval=False,
                         marker_kw={},
                         line_kw={},
                         interval_kw={},
                         resampling_kw={}):
        """
        Plot expression dynamics for specified channel, resampling from discrete subpopulations of cells.

        Args:

            channel (str) - expression channel

            ax (mpl.axes.AxesSubplot) - if None, create axes

            average (bool) - if True, add moving average

            interval - if True, add confidence interval for moving average

            line_kw (dict) - keyword arguments for line formatting

            interval_kw (dict) - keyword arguments for interval formatting

            resampling_kw (dict) - keyword arguments for disc resampler

        Returns:

            ax (mpl.axes.AxesSubplot)

        """

        # sort values inplace
        self.sort('t')

        # resample discs and cells within them
        time = DiscResampler(self, 't', **resampling_kw).mean
        resampler = DiscResampler(self, channel, **resampling_kw)
        mean = resampler.mean
        lower, upper = resampler.confidence_interval

        # construct interval plot
        interval_plot = IntervalPlot(time, lower, upper, mean, ax=ax)

        # plot dynamics
        interval_plot.plot(average=average,
                     interval=interval,
                     line_kw=line_kw,
                     interval_kw=interval_kw)

        return interval_plot.ax

    def scatterplot(self, x, y,
                      color='grey',
                      s=5,
                      alpha=0.5,
                      fraction=False,
                      ax=None):
        """
        Create XY scatterplot of two fluorescence channels.

        Args:

            x, y (str or int) - channels used for x and y axes

            color (str) - marker color

            s (float) - marker size

            alpha (float) - transparency of markers

            fraction (bool) - if True, annotate fraction above midline

            ax (mpl.axes.AxesSubplot) - if None, create figure

        Returns:

            ax (mpl.axes.AxesSubplot)

        """

        # get string representation of channel names
        x, y = format_channel(x), format_channel(y)

        # create figure
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))

        # format axes
        ax.set_xlim(0, 2), ax.set_ylim(0, 2)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.tick_params(labelsize=10)
        ax.grid(True)

        # scatter data
        ax.scatter(self.data[x], self.data[y], c=color, s=s, alpha=alpha, lw=0)

        # add fraction above midline
        if fraction:
            ratio = self.data[y]/self.data[x]
            self.annotate_fraction(ax, ratio, p=2.5)

        return ax

    @staticmethod
    def annotate_fraction(ax, ratio, p=2.5):
        """
        Add fraction of cells above midline.

        Args:

            ax (mpl.axes.AxesSubplot)

            ratio (array like) - vector of ratios

            p (float) - text position relative to center line

        """
        fraction = sum(ratio >= 1) / len(ratio)
        ax.text(p, p+0.5, '{:2.1%}'.format(fraction), ha='right', fontsize=8)
        ax.text(p, p-0.5, '{:2.1%}'.format(1-fraction), ha='left', fontsize=8)

    def plot_spectrogram(self, channel,
                     periods=None,
                     ymax=None,
                     ax=None, **kwargs):
        """
        Plot Lomb Scargle periodogram.

        Args:

            channel (str or int) - expression channel

            periods (array like) - spectral frequencies to be tested

            ymax (float) - max spectral power

            ax (mpl.axes.AxesSubplot)

            kwargs: spectrogram visualization keywords

        Returns:

            ax (mpl.axes.AxesSubplot)

        """

        # get string representation of channel name
        channel = format_channel(channel)

        # compile spectrogram
        precursors = self.select_cell_type('pre')
        spectrogram = Spectrogram(precursors.data.centroid_y.values, precursors.data[channel].values, periods=periods)

        # plot power spectrum
        ax = spectrogram.simple_visualization(ax=ax, ymax=ymax, **kwargs)

        return ax
