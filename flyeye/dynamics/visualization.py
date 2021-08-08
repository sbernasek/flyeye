import matplotlib.pyplot as plt
from .averages import savgol, get_rolling_mean, get_rolling_mean_interval


def plot_mean(x, y, ax,
              label=None,
              ma_type='sliding',
              window_size=100,
              resolution=1,
              line_color='k',
              line_width=1,
              line_alpha=1,
              linestyle=None,
              markersize=2,
              smooth=False,
              **kw):
    """
    Plot moving average.

    Args:

        x, y (array like) - timeseries data

        ax (matplotlib.axes.AxesSubplot) - axis which to which line is added

        label (str) - data label

        ma_type (str) - type of average used, either sliding, binned, or savgol

        window_size (int) - size of window

        resolution (int) - sampling interval

        line_color, line_width, line_alpha, linestyle - formatting parameters

        smooth (bool) - if True, apply secondary savgol filter

    Returns:

        line (matplotlib.lines.Line2D)

    """

    # get moving average (skip first point to avoid outliers)
    if ma_type == 'savgol':
        x_av = x[1:]
        y_av = savgol(y, window_size=window_size, polyorder=1)[1:]
    else:
        if ma_type == 'binned':
            resolution = window_size
        x_av = get_rolling_mean(x, window_size=window_size, resolution=resolution)
        y_av = get_rolling_mean(y, window_size=window_size, resolution=resolution)

    # get line and dashstyles
    if linestyle == None:
        linestyle = 'solid'
    dashstyles = {'solid': (None, None), 'dashed': (2.0, 2.0)}
    dashstyle = dashstyles[linestyle]

    # apply secondary smoothing (for visualization)
    if smooth:
        sw = int(window_size/5)
        if sw > 1:
            y_av = savgol(y_av, window_size=sw, polyorder=1)

    # plot line
    line = ax.plot(x_av, y_av,
                   linestyle=linestyle, dashes=dashstyle,
                   lw=line_width, color=line_color, alpha=line_alpha,
                   label=label, markersize=markersize, **kw)

    return line


def plot_mean_interval(x, y, ax,
                       ma_type='sliding',
                       window_size=100,
                       resolution=10,
                       nbootstraps=1000,
                       confidence=95,
                       color='grey',
                       alpha=0.25,
                       error_bars=False,
                       lw=0.):
    """
    Adds confidence interval for line average (sliding window or binned) to existing axes.

    Args:

        x, y (array like) - data

        ax (axes) - axis which to which line is added

        ma_type (str) - type of average used, either 'sliding' or 'binned'

        window_size (int) - size of sliding window or bin (num of cells)

        interval_resolution (int) - sampling resolution for confidence interval

        nbootstraps (int) - number of bootstraps

        confidence (float) - confidence interval, between 0 and 100

        color, alpha - formatting parameters

    """

    if ma_type == 'binned':
        resolution = window_size
    x_av = get_rolling_mean(x, window_size=window_size, resolution=resolution)
    y_av = get_rolling_mean(y, window_size=window_size, resolution=resolution)
    interval = get_rolling_mean_interval(y, window_size=window_size, resolution=resolution, nbootstraps=nbootstraps, confidence=confidence)
    y_lower, y_upper = interval.T
    _ = ax.fill_between(x_av, y_lower, y_upper,
                        color=color, alpha=alpha, lw=lw)

    if error_bars == True:
        ax.errorbar(x_av, y_av, yerr=[y_av-y_lower, y_upper-y_av], fmt='-o', color=color)


class TimeseriesPlot:
    """
    Object describes a 1D timeseries.

    Attributes:

        x (np.ndarray) - independent variable

        y (np.ndarray) - dependent variable

        ax (matplotlib.axes.AxesSubplot)

    """

    def __init__(self, x, y, ax=None):
        """
        Instantiate a 1D timeseries.

        Args:

            x (np.ndarray) - independent variable

            y (np.ndarray) - dependent variable

            ax (matplotlib.axes.AxesSubplot)

        """

        self.x = x
        self.y = y

        # set axis
        if ax is None:
            ax = self.create_figure()
        self.ax = ax

    def create_figure(self):
        """ Instantiate figure. """
        fig, ax = plt.subplots(ncols=1, figsize=(3, 2))
        ax.set_xlim(self.x.min(), self.x.max())
        ax.set_ylim(0, 1.1*self.y.max())
        ax.set_xlabel('Time (h)'),
        ax.set_ylabel('Expression (a.u.)')
        return ax

    def scatter(self,
                 color='k',
                 alpha=1,
                 s=1,
                 rasterized=False,
                 **additional):
        """
        Scatterplot markers for x and y data.

        Args:

            color (str) - marker color

            alpha (float) - marker alpha

            s (float) - marker size

            rasterized (bool) - if True, rasterize markers

        """
        marker_kw = dict(color=color, s=s, alpha=alpha, lw=0, rasterized=rasterized)
        _ = self.ax.scatter(self.x, self.y, **marker_kw, **additional)

    def average(self,
            ma_type='savgol',
            window_size=100,
            resolution=1,
            smooth=True,
            color='k',
            alpha=1,
            lw=1,
            linestyle=None, 
            **additional
            ):
        """
        Plot moving average of x and y data.

        Args:

            ma_type (str) - type of average, 'savgol', 'sliding', or 'binned'

            window_size (int) - size of sliding window or bin (num of cells)

            resolution (int) - sampling resolution for confidence interval

            smooth (bool) - if True, apply secondary savgol filter

            color, alpha, lw, linestyle - formatting parameters

        """

        ma_kw = dict(ma_type=ma_type, window_size=window_size, resolution=resolution, smooth=smooth)
        line_kw = dict(line_color=color, line_alpha=alpha, line_width=lw, linestyle=linestyle)

        if len(self.y) > window_size:
            _ = plot_mean(self.x, self.y, ax=self.ax, **ma_kw, **line_kw, **additional)

    def interval(self,
            ma_type='sliding',
            window_size=100,
            resolution=25,
            nbootstraps=1000,
            confidence=95,
            color='k',
            alpha=0.5,
            **additional):
        """
        Plot confidence interval for moving average of x and y data.

        Args:

            ma_type (str) - type of moving average, 'sliding' or 'binned'

            window_size (int) - size of sliding window or bin (num of cells)

            resolution (int) - sampling resolution for confidence interval

            nbootstraps (int) - number of bootstraps

            confidence (float) - confidence interval, between 0 and 100

            color, alpha - formatting parameters

        """

        # define moving average keyword arguments
        ma_kw = dict(ma_type=ma_type,
                     window_size=window_size,
                     resolution=resolution,
                     nbootstraps=nbootstraps,
                     confidence=confidence)

        # define interval shading keyword arguments
        shade_kw = dict(color=color, alpha=alpha)

        # plot confidence interval
        if len(self.y) > window_size:
            plot_mean_interval(self.x,
                                  self.y,
                                  ax=self.ax,
                                  **ma_kw,
                                  **shade_kw)

    def plot(self,
             scatter=False,
             average=True,
             interval=False,
             marker_kw={},
             line_kw={},
             interval_kw={},
             ma_kw={}):
        """
        Plot timeseries data.

        Args:

            scatter (bool) - if True, add datapoints

            average (bool) - if True, add moving average

            interval (bool) - if True, add moving average interval

            marker_kw (dict) - keyword arguments for marker formatting

            line_kw (dict) - keyword arguments for line formatting

            interval_kw (dict) - keyword arguments for interval formatting

            ma_kw (dict) - keyword arguments for moving average

        """

        # add scattered data
        if scatter:
            self.scatter(**marker_kw)

        # add moving average
        if average:
            self.average(**ma_kw, **line_kw)

        # add confidence interval for moving average
        if interval:
            self.interval(**ma_kw, **interval_kw)


class IntervalPlot(TimeseriesPlot):
    """
    Object describes the 95% confidence interval for a 1D timeseries.

    Attributes:

        x (np.ndarray) - independent variable

        y_lower (np.ndarray) - lower bound for dependent variable

        y_upper (np.ndarray) - upper bound for dependence variable

        y (np.ndarray) - mean or median value for dependent variable

        ax (matplotlib.axes.AxesSubplot)

    """

    def __init__(self, x, y_lower, y_upper, y=None, ax=None):
        """
        Instantiate a 1D timeseries.

        Args:

            x (np.ndarray) - independent variable

            y_lower (np.ndarray) - lower bound for dependent variable

            y_upper (np.ndarray) - upper bound for dependence variable

            y (np.ndarray) - median value for dependent variable

            ax (matplotlib.axes.AxesSubplot)

        """

        self.x = x
        self.y = y
        self.y_lower = y_lower
        self.y_upper = y_upper

        # set axis
        if ax is None:
            ax = self.create_figure()
        self.ax = ax

    def average(self,
            smooth=True,
            color='k',
            alpha=1,
            lw=1,
            linestyle=None, **addtl):
        """
        Plot moving average of x and y data.

        Args:

            smooth (bool) - if True, apply first-order savgol filter

            color, alpha, lw, linestyle - formatting parameters

        """
        line_kw = dict(color=color, lw=lw, alpha=alpha, linestyle=linestyle)
        self.ax.plot(self.x, self.y, **line_kw)

    def interval(self,
            color='k',
            alpha=0.5,
            **additional):
        """
        Plot confidence interval for moving average of x and y data.

        Args:

            color, alpha - formatting parameters

        """
        shade_kw = dict(color=color, alpha=alpha)
        self.ax.fill_between(self.x, self.y_lower, self.y_upper, **shade_kw)

    def plot(self,
             average=True,
             interval=False,
             line_kw={},
             interval_kw={}):
        """
        Plot timeseries data.

        Args:

            average (bool) - if True, add moving average

            interval (bool) - if True, add moving average interval

            line_kw (dict) - keyword arguments for line formatting

            interval_kw (dict) - keyword arguments for interval formatting

        """

        # add moving average
        if average:
            self.average(**line_kw)

        # add confidence interval for moving average
        if interval:
            self.interval(**interval_kw)
