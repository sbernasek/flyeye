import numpy as np
from matplotlib import pyplot as plt
from astroML.time_series import lomb_scargle, lomb_scargle_bootstrap


class Spectrogram:
    """
    Object for spectral decomposition of a 1-D timeseries via Lomb-Scargle periodogram. Internal functions are based on AstroML library.

    Attributes:

        t (np.ndarray) - timepoints

        y (np.ndarray) - values

        dy (np.ndarray) - estimated measurement error

        periods (np.ndarray) - oscillation periods tested (same units as t)

        omegas (np.ndarray) - oscillation frequencies tested

        PS (np.ndarray) - spectral power of each frequency

        power (np.ndarray) - max. spectral power

        dominant_period (float) - oscillation period of max. spectral power

        dominant_frequency (float) - oscillation frequency of max. spectral power

    """

    def __init__(self, t, y,
                 dy=None,
                 periods=None):
        """
        Instantiate object for spectral decomposition of a 1-D timeseries via Lomb-Scargle periodogram.

        Args:

            t (np.ndarray) - timepoints

            y (np.ndarray) - values

            dy (np.ndarray) - estimated measurement error

            periods (np.ndarray) - oscillation periods tested (same units as t)

        """

        # order timeseries by time
        sort_indices = np.argsort(t)
        self.t = t[sort_indices]
        self.y = y[sort_indices]

        # determine or estimate measurement error
        if dy is not None:
            self.dy = dy[sort_indices]
        else:
            self.dy = np.std(self.y)/np.sqrt(self.y.size) * 1

        # define spectral components
        if periods is None:
             periods = 10 ** np.linspace(np.log10(40), np.log10(200), 1000)
        self.periods = periods
        self.omegas = 2*np.pi/periods

        # compute periodogram
        self.evaluate_periodogram()
        self.power = self.PS.max()
        self.dominant_period = self.periods[self.PS.argmax()]
        self.dominant_frequency = self.omegas[self.PS.argmax()]

    @staticmethod
    def _periodogram(t, y, dy, omegas):
        """
        Evaluate periodogram.

        Args:

            t (np.ndarray) - timepoints

            y (np.ndarray) - values

            dy (np.ndarray) - estimated measurement error

            omega (np.ndarray) - spectral frequencies tests

        Returns:

            PS (np.ndarray) - normalized power spectrum

        """
        kw = dict(generalized=True, subtract_mean=True)
        PS = lomb_scargle(t, y, dy, omegas, **kw)
        return PS

    def evaluate_periodogram(self):
        """ Evaluate periodogram. """
        self.PS = self._periodogram(self.t, self.y, self.dy, self.omegas)

    @staticmethod
    def _compute_thresholds(t, y, dy, omegas,
                                     confidence=None,
                                     nbootstraps=1000):
        """
        Determine periodicity significance thresholds. Thresholds are obtained by repeatedly subsampling data, then compiling a distribution of maximum power levels detected in each sample.

        Args:

            t (np.ndarray) - timepoints

            y (np.ndarray) - values

            dy (np.ndarray) - estimated measurement error

            omega (np.ndarray) - spectral frequencies tests

            confidence (array like) - confidence levels to assess, length C

            nbootstraps (int) - number of boostrap samples

        Returns:

            thresholds (np.ndarray) - spectral power thresholds, length C

        """

        if confidence is None:
            confidence = [95, 99, 99.9]

        # get seed for random state
        seed = np.random.randint(0, 1000)

        # assemble spectral powers for null models (random subsamples)
        null_ps = lomb_scargle_bootstrap(t, y, dy, omegas, generalized=True, N_bootstraps=nbootstraps, random_state=seed)

        # compute significance thresholds
        thresholds = {c: np.percentile(null_ps, c) for c in confidence}

        return thresholds

    def compute_thresholds(self, confidence=None, nbootstraps=1000):
        """
        Determine significance thresholds.

        Args:

            confidence (array like) - confidence levels to assess, length C

            nbootstraps (int) - number of boostrap samples

        Returns:

            thresholds (np.ndarray) - spectral power thresholds, length C

        """
        return self._compute_thresholds(self.t, self.y, self.dy, self.omegas, confidence=confidence, nbootstraps=nbootstraps)

    @staticmethod
    def _plot_samples(ax, t, y, dy):
        """
        Plot timeseries samples.

        Args:

            ax (matplotlib.axes.AxesSubplot)

            t (array like) - timeponts

            y (array like) - values

            dy (array like) - value uncertainties

        Returns:

            ax (matplotlib.axes.AxesSubplot)

        """
        ax.errorbar(t, y, dy, fmt='.k', lw=1, ecolor='gray')
        ax.set_xlabel('distance')
        ax.set_ylabel('samples')
        ax.set_xlim(t.min(), t.max()+t.ptp()/10)
        return ax

    def plot_samples(self, ax):
        """
        Plot timeseries samples.

        Args:

            ax (mpl.axes.AxesSublot)

        """
        self._plot_samples(ax, self.t, self.y, self.dy)

    @staticmethod
    def _plot_spectrogram(ax, xvals, PS,
                          thresholds={},
                          xaxis='period',
                          xunits='px',
                          ymax=None,
                          color='k',
                          labelsize=10):
        """
        Plot spectrogram.

        Args:

            ax (mpl.axes.AxesSublot)

            xvals (np.ndarray) - x values (periods or frequencies)

            PS (np.ndarray) - normalized spectral powers

            thresholds (dict) - keys are confidence levels, values are thresholds

            xaxis (str) - quantity for x axis, either 'period' or 'frequency'

            xunits (str) - units for x-axis label

            ymax (float) - upper limit for y axis (max spectral power)

            color (str) - line color

            labelsize (int) - axis labelsize

        Returns:

            ax (mpl.axes.AxesSublot)

        """

        # plot spectral powers
        ax.plot(xvals, PS, '-', c=color, lw=1, zorder=1)

        # plot significant thresholds
        for sig, threshold in thresholds.items():
            ax.plot([xvals.min(), xvals.max()], [threshold, threshold], linestyle='dashed', c='red', linewidth=0.5, dashes=(2, 2))

        # determine ylim
        if ymax is None:
            if len(thresholds.values()) == 0:
                ymax = 1.25*PS.max()
            else:
                ymax = max(1.25*max(thresholds.values()), PS.max())

        # format plot
        ax.set_xlim(xvals.min(), xvals.max())
        ax.set_ylim(0., ymax)
        ax.set_xlabel(r'{:s} ({:s})'.format(xaxis.capitalize(), xunits), fontsize=labelsize)
        ax.set_ylabel('Power', fontsize=labelsize)

        return ax

    def plot_spectrogram(self, ax,
                         confidence=[99],
                         nbootstraps=None,
                         xaxis='period',
                         annotate=True,
                         **kwargs):
        """
        Plot spectrogram.

        Args:

            ax (mpl.axes.AxesSublot)

            confidence (array like) - confidence levels to be added

            nbootstraps (int) - number of bootstrap samples

            xaxis (str) - quantity for x axis, either 'period' or 'frequency'

            annotate (bool) - if True, label peak spectral power

            kwargs: plot formatting keyword arguments

        Returns:

            ax (mpl.axes.AxesSublot)

        """

        # compute significance thresholds
        if nbootstraps is not None:
            thresholds = self.compute_thresholds(confidence, nbootstraps)
        else:
            thresholds = {}

        # determine x axis
        if xaxis == 'period':
            xvals = self.periods
            xunits = 'px'
        elif xaxis == 'frequency':
            xvals = self.omegas
            xunits = 'rad/s'

        # plot spectrogram
        ax =  self._plot_spectrogram(ax,
                                     xvals,
                                     self.PS,
                                     thresholds=thresholds,
                                     xaxis=xaxis,
                                     xunits=xunits,
                                     **kwargs)

        # annotate peak spectral power if it falls above lowest threshold
        if annotate:
            lowest_threshold = thresholds[min(thresholds.keys())]
            if self.power >= lowest_threshold:
                self.add_power(ax, xaxis=xaxis)

        return ax

    def add_power(self, ax, xaxis='period'):
        """
        Annotate peak spectral power.

        Args:

            ax (mpl.axes.AxesSublot)

            xaxis (str) - quantity for x axis, either 'period' or 'frequency'

        """

        # get label position
        if xaxis == 'period':
            xp = self.periods.max() - .03*(self.periods.ptp())
            xunits = 'px'
            xpeak = self.dominant_period
        else:
            xp = self.omegas.max() - .03*(self.omegas.ptp())
            xunits = 'rad/s'
            xpeak = self.dominant_frequency
        yp = self.power + .1

        # add label
        txt = '*{:0.2f} at {:0.0f} {:s}'.format(self.power, xpeak, xunits)
        ax.text(xp, yp, txt, color='k', ha='right', va='bottom')

    def simple_visualization(self,
                             ax=None,
                             confidence=[99],
                             nbootstraps=1000,
                             xaxis='period',
                             annotate=True,
                             **kwargs):

        """
        Plot spectrogram.

        Args:

            ax (mpl.axes.AxesSublot) - if None, create figure

            confidence (array like) - confidence levels to be added

            nbootstraps (int) - number of bootstrap samples

            xaxis (str) - quantity for x axis, either 'period' or 'frequency'

            annotate (bool) - if True, label peak spectral power

            kwargs: plot formatting keyword arguments

        Returns:

            ax (mpl.axes.AxesSublot)

        """

        # define axes
        if ax is None:
            fig = plt.figure(figsize=(3, 3))
            fig.subplots_adjust(left=0.1, right=0.9, hspace=0.25)
            ax = fig.add_subplot(211)

        # add spectrogram
        self.plot_spectrogram(ax, confidence, nbootstraps, xaxis, annotate, **kwargs)

        return ax

    def full_visualization(self,
                           confidence=[99],
                           nbootstraps=1000,
                           xaxis='period'):
        """
        Plot spectrogram and corresponding signal on adjacent axes.

        Args:

            confidence (array like) - confidence levels to be added

            nbootstraps (int) - number of bootstrap samples

            xaxis (str) - quantity for x axis, either 'period' or 'frequency'

        """

        # Create figure
        fig = plt.figure(figsize=(5, 3.75))
        fig.subplots_adjust(left=0.1, right=0.9, hspace=0.25)
        ax0 = fig.add_subplot(211)
        ax1 = fig.add_subplot(212, xscale=xscale)

        # set axis scale
        if xaxis == 'period':
            xscale = 'linear'
        else:
            xscale = 'log'

        # plot data
        self.plot_samples(ax0)

        # plot spectrogram
        self.plot_spectrogram(ax1, confidence, nbootstraps, xaxis)

        plt.tight_layout()
