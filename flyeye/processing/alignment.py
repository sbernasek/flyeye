import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from functools import reduce
from operator import add
from math import ceil, floor
from scipy.interpolate import interp1d

from ..dynamics.averages import get_running_mean, get_rolling_mean


class Alignment:
    """
    Object for alignment of a 1-D timeseries with another.

    Attributes:

        t, x (np.ndarray) - timeseries to be aligned

        t_ref, x_ref (np.ndarray) - reference to which timeseries is aligned

        window_size (int) - window size used for local smoothing

        lag (float) - computed time shift for optimal alignment

        score (float) - computed metric for quality of alignment

    """

    def __init__(self, t, x,
                 t_ref=None, x_ref=None,
                 metric='crosscorrelation',
                 window_size=10):
        """
        Instantiate object for alignment of one timeseries with another.

        Args:

            t (np.ndarray) - timepoints of timeseries to be aligned

            x (np.ndarray) - values of timeseries to be aligned

            t_ref (np.ndarray) - timepoints of reference

            x_ref (np.ndarray) - values of reference

            metric (str) - name of alignment criterion

            window_size (int) - window size used for local smoothing

        """

        self.t, self.x = self._initialize_position(t, x)
        self.window_size = window_size
        if t_ref is None and x_ref is None:
            t_ref, x_ref = t, x
        self.t_ref, self.x_ref = self._initialize_position(t_ref, x_ref)
        lag, score = self.align(metric=metric, window_size=window_size)
        self.lag = lag
        self.score = score

    def _get_alignnment_object(self):
        """ Return new Alignment instance. """
        return Alignment(self.t, self.x, self.t_ref, self.x_ref)

    @staticmethod
    def _sort(t, x):
        """
        Sort timeseries by time.

        Args:

            t (np.ndarray) - timepoints

            x (np.ndarray) - values

        Returns:

            ts, xs (np.ndarray) - time-sorted timepoints and values

        """
        sort_ind = np.argsort(t)
        return t[sort_ind], x[sort_ind]

    @staticmethod
    def _align_zero(x):
        """ Shift values so minimum is zero. """
        x -= x.min()
        return x

    @classmethod
    def _initialize_position(cls, t_, x_):
        """ Sort vector and shift to zero start position. """
        t, x = cls._sort(t_, x_)
        return cls._align_zero(t), x

    @staticmethod
    def _get_regular_timepoints(t, dt=0.1):
        """ Construct array of regularly sampled timepoints. """
        return np.arange(t.min(), t.max(), step=dt)

    @staticmethod
    def _interpolate(t, x, t_regular):
        """ Interpolate values onto regularly sampled timepoints. """
        interp = interp1d(t, x)
        return interp(t_regular)

    @classmethod
    def _smooth(cls, t, x,
                window_size=10,
                dt=0.1):
        """ Return moving average interpolated onto regular timepoints. """

        # evaluate moving average
        ts = get_rolling_mean(t, window_size=window_size)
        xs = get_rolling_mean(x, window_size=window_size)

        # interpolate onto regulat timepoints
        tsr = cls._get_regular_timepoints(ts, dt=dt)
        xsr = cls._interpolate(ts, xs, tsr)

        return tsr, xsr

    @staticmethod
    def _zscore(x):
        """ Z-score values. """
        return (x - x.mean())/x.std()

    @staticmethod
    def _unit_scale(x):
        """ Re-scale by maximum value. """
        return x/x.max()

    @classmethod
    def _get_crosscorrelation(cls, x, y):
        """ Evaluate normalized crosscorrelation between z-scored vectors. """
        xz, yz = cls._zscore(x), cls._zscore(y)
        return np.correlate(xz, yz, mode='full') / min(len(x), len(y))

    @classmethod
    def _get_crossvariogram(cls, x, y):
        """ Evaluate nonnormalized correlation between unit-scaled vectors. """
        xu, yu = cls._unit_scale(x), cls._unit_scale(y)
        return np.correlate(xu, yu, mode='full') / min(len(x), len(y))

    @staticmethod
    def _slide_with_extrapolation(f, x, y):
        """ Apply function over slices in which y rolls across x. Outer-bounds of y are filled with its min/max values. """
        n, m = len(x), len(y)
        x_ref = np.hstack((np.mean(y[0:5])*np.ones(n), y, np.mean(y[-5:])*np.ones(n)))
        return np.array([f(x, x_ref[-(i+1+n):-(i+1)]) for i in range(n+m-1)])

    @staticmethod
    def _slide_full(func, x, y):
        """ Apply function over slices in which y rolls over x. Behavior is quivalent to 'full' mode in np.correlate. """
        n, m = len(x), len(y)
        evals = []
        for i in range(n+m-1):
            if i < min(n, m):
                func_eval = func(y[-(i+1):], x[:(i+1)])
                evals.append(func_eval)
            elif m >= n:
                func_eval = func(y[(-i-1):(n-i)-1], x[(i - m-n+1):])
                evals.append(func_eval)
            else:
                func_eval = func(x_ref[:-(i - n - m+1)], x[-(m-i)+1:-(-i-1)])
                evals.append(func_eval)
        return np.array(evals)

    @classmethod
    def _get_inverse_distance(cls, x, y):
        """ Evaluate inverse distance between unit-scaled vectors. """
        get_inverse_distance = lambda x1, x2: 1/((abs(x1-x2)).sum())
        scores = cls._slide_with_extrapolation(get_inverse_distance, cls._unit_scale(x), cls._unit_scale(y))
        return scores

    @classmethod
    def _get_inverse_squared_distance(cls, x, y):
        """ Evaluate squared inverse distance between unit-scaled vectors. """
        get_inverse_squared_distance = lambda x1, x2: 1/(((x1-x2)**2).sum())
        scores = cls._slide_with_extrapolation(get_inverse_squared_distance, cls._unit_scale(x), cls._unit_scale(y))
        return scores

    @staticmethod
    def _get_lag_vector(t, t_ref, **kwargs):
        """
        Returns ordered vector of lag times. Lags are shifted by minimums in order to account for first timepoint being shifted above zero by smoothing.
        """
        return np.hstack((t_ref[::-1][:-1]-t.min(), t_ref.min()-t.min(), -t[1:]+ t_ref.min()))

    @staticmethod
    def _maximize(lags, scores):
        """ Returns lag that maximizes scoring metric, along with corresponding score. """
        index = scores.argmax()
        return lags[index], scores[index]

    def get_smoothed_scores(self,
                            metric='crosscorrelation',
                            window_size=10,
                            dt=0.1):
        """ Returns smoothed alignment scores. """

        # interpolate rolling averages onto regular grid
        t_ref, x_ref = self._smooth(self.t_ref, self.x_ref, window_size=window_size, dt=dt)
        t, x = self._smooth(self.t, self.x, window_size=window_size, dt=dt)

        # evaluate scores
        if metric == 'crosscorrelation':
            scores = self._get_crosscorrelation(x, x_ref)
        elif metric == 'crossvariogram' or metric == 'variogram':
            scores = self._get_crossvariogram(x, x_ref)
        elif metric == 'inverse_distance':
            scores = self._get_inverse_distance(x, x_ref)
        elif metric == 'inverse_squared_distance':
            scores = self._get_inverse_squared_distance(x, x_ref)
        else:
            raise ValueError('Alignment metric not recognized.')

        return t, t_ref, scores

    def align(self,
              metric='crosscorrelation',
              window_size=50,
              dt=0.1):
        """
        Run alignment to determine optimal lag.

        Args:

            metric (str) - name of alignment criterion

            window_size (int) - window size for local smoothing

            dt (float) - resolution

        Returns:

            lag (float) - time shift that maximizes quality of alignment

            score (float) - quality of alignment

        """

        # get smoothed metric vectors
        t, t_ref, scores = self.get_smoothed_scores(metric=metric, window_size=window_size, dt=dt)

        # find lag that maximizes metric
        lags = self._get_lag_vector(t, t_ref, window_size=window_size)
        lag, score = self._maximize(lags, scores)

        return lag, score

    def plot_scores(self,
                    metric='crosscorrelation',
                    ax=None,
                    window_size=None,
                    dt=0.1):
        """
        Plot quality of alignment versus lag time.

        Args:

            metric (str) - name of alignment criterion

            ax (matplotlib.axes.AxesSubplot)

            window_size (int) - window size for local smoothing

            dt (float) - resolution

        """

        # create axes if none provided
        if ax is None:
            fig, ax = plt.subplots(figsize=(4, 2))

        if window_size is None:
            window_size = self.window_size

        # get correlations
        t, t_ref, scores = self.get_smoothed_scores(metric=metric, window_size=window_size, dt=dt)
        lags = self._get_lag_vector(t, t_ref, window_size=window_size)

        # plot correlation as a function of lag
        ax.plot(lags, scores, '-k', linewidth=3, markersize=0)

        # add maximum
        ax.scatter(*self._maximize(lags, scores), s=25, color='red', zorder=99)

        # format plot
        ax.set_xlabel('Lag', fontsize=12)
        ax.set_ylabel(metric.title(), fontsize=12)
        if metric=='crosscorrelation' or metric=='crossvariogram':
            ax.set_ylim(-0.7, 1.2)
            ax.set_yticks([-0.5, 0, 0.5, 1])
            ax.axhline(1, linestyle='--', color='k')
            ax.axhline(0, linestyle='--', color='k')

    def plot_alignment(self,
                       scatter=True,
                       trend=True,
                       ax=None,
                       window_size=100):
        """
        Plot aligned time series.

        Args:

            scatter (bool) - if True, show individual samples

            trend (bool) - if True, show moving average

            ax (matplotlib.axes.AxesSubplot)

            window_size (int) - window size for local smoothing

        """

        # add datapoints to plot
        if ax is None:
            fig, ax = plt.subplots(figsize=(4, 2))

        if scatter:
            ax.plot(self.t + self.lag, self.x, '.r', alpha=0.1)
            ax.plot(self.t_ref, self.x_ref, '.k', alpha=0.1)

        # add moving averages to plot
        if trend:
            smooth = lambda x: get_running_mean(x, window_size=window_size)
            red_line, = ax.plot(smooth(self.t + self.lag), smooth(self.x), '-r', alpha=1)
            black_line, = ax.plot(smooth(self.t_ref), smooth(self.x_ref), '-k', alpha=1)

            # format plot
            ax.legend([black_line, red_line], ['reference', 'aligned'], loc=0, frameon=False)


class CellsAlignment(Alignment):
    """
    Object for alignment of one group of cells with another.

    Attributes:

        data (pd.DataFrame) - aligned cells data

    Inherited attributes:

        t, x (np.ndarray) - timeseries to be aligned

        t_ref, x_ref (np.ndarray) - reference to which timeseries is aligned

        window_size (int) - window size used for local smoothing

        lag (float) - computed time shift for optimal alignment

        score (float) - computed metric for quality of alignment

    """

    def __init__(self, data,
                 reference_data=None,
                 channel='ch1_normalized',
                 basis='t',
                 metric='crosscorrelation',
                 window_size=None):
        """
        Instantiate object for alignment of one group of cells with another.

        Args:

            data (pd.DataFrame) - aligned cells data

            reference_data (pd.DataFrame) - reference cells data

            channel (str or int) - cells attribute used as alignment values

            basis (str) - cells attribute  used as alignment timepoints

            metric (str) - name of alignment criterion

            window_size (int) - window size used for local smoothing (# of cells)
        """

        # extract timeseries
        self.data = data
        t, x = data[basis].values, data[channel].values
        t_ref, x_ref = reference_data[basis].values, reference_data[channel].values

        # set window size
        if window_size is None:
            window_size = 10

        # perform alignment
        Alignment.__init__(self, t, x, t_ref, x_ref, metric=metric, window_size=window_size)


class DiscAlignment(CellsAlignment):
    """
    Object for alignment of one Disc instance with another.

    Attributes:

        disc (data.discs.Disc) - copy of aligned disc instance

        hours_per_pixel (float) - distance to time scaling factor

        space_lag (float) - computed time shift for optimal alignment (x-units)

        time_lag (float) - computed time shift for optimal alignment (t-units)

    Inherited attributes:

        t, x (np.ndarray) - timeseries to be aligned

        t_ref, x_ref (np.ndarray) - reference to which timeseries is aligned

        window_size (int) - window size used for local smoothing

        lag (float) - computed time shift for optimal alignment

        score (float) - computed metric for quality of alignment

        data (pd.DataFrame) - aligned cells data

    """

    def __init__(self, disc, reference_disc,
                 channel='ch1_normalized',
                 metric='crosscorrelation',
                 window_size=None):
        """
        Instantiate object for alignment of one Disc instance with another.

        Args:

            disc (data.discs.Disc) - aligned disc

            reference_disc (data.discs.Disc) - reference disc

            channel (str) - cells attribute used as alignment values

            metric (str) - name of alignment criterion

            window_size (int) - window size used for local smoothing (# cells)

        """

        self.hours_per_pixel = disc.hours_per_pixel
        self.disc = deepcopy(disc)

        # select precursor dataframes
        precursor_data = disc.select_cell_type('pre').data
        reference_precursor_data = reference_disc.select_cell_type('pre').data

        # if no precursors were found, check for 'wpre'
        if len(precursor_data) == 0:
            precursor_data = disc.select_cell_type('wpre').data
            reference_precursor_data = reference_disc.select_cell_type('wpre').data

        # initialize alignment object (discs must be aligned on time)
        CellsAlignment.__init__(self, precursor_data, reference_data=reference_precursor_data, channel=channel, basis='t', metric=metric,  window_size=window_size)

        # update time lags
        self.space_lag = self.lag / disc.hours_per_pixel
        self.time_lag = self.lag

    def __call__(self):
        """ Return aligned disc """
        return self.get_aligned_disc()

    @staticmethod
    def apply_to_disc(disc, space_lag, time_lag):
        """
        Apply lag to copy of disc.

        Args:

            disc (data.discs.Disc)

            space_lag (float) - lag to be applied (x-units)

            time_lag (float) - lag to be applied (t-units)

        Returns:

            disc (data.discs.Disc) - copy of disc with lag applied

        """
        disc = deepcopy(disc)
        disc.data.centroid_x += (space_lag - disc.data.centroid_x.min())
        disc.data.t += (time_lag - disc.data.t.min())
        return disc

    def get_aligned_disc(self):
        """ Return aligned copy of disc. """
        disc = deepcopy(self.disc)
        disc = self.apply_to_disc(disc, space_lag=self.space_lag, time_lag=self.time_lag)
        return disc


class ExperimentAlignment:
    """
    Object for the time-alignment of all discs within an experiment.

    Attributes:

        experiment (data.experiments.Experiment) - copy of experiment

        scores (dict) - keys are disc IDs, values are quality of alignment

    """

    def __init__(self, experiment, **kwargs):
        """
        Instantiate object for the time-alignment of all discs within an experiment. The first disc serves as the reference.

        Args:

            experiment (data.experiments.Experiment)

            kwargs: keyword arguments for DiscAlignment instantiation

        """

        self.experiment = deepcopy(experiment)
        self.align_discs(**kwargs)
        self.scores = {disc: alignment.score for disc, alignment in self.alignments.items()}

    @classmethod
    def _align_discs(cls, experiment, **kwargs):
        """
        Align all discs with first disc.

        Args:

            experiment (data.experiments.Experiment)

            kwargs: keyword arguments for DiscAlignment instantiation

        Returns:

            experiment (data.experiments.Experiment) - copy of experiment

            alignments (dict) - {disc ID: Alignment} pairs

        """
        reference_disc = experiment.discs[0]
        alignments = {}
        for disc_id in range(len(experiment.discs)):
            alignment = DiscAlignment(experiment.discs[disc_id], reference_disc, **kwargs)
            alignments[disc_id] = alignment
            experiment.discs[disc_id] = alignment.get_aligned_disc()
        return experiment, alignments

    def align_discs(self, **kwargs):
        """
        Align all discs with first disc.

        kwargs: keyword arguments for DiscAlignment instantiation
        """
        experiment, alignments = self._align_discs(self.experiment, **kwargs)
        self.experiment = experiment
        self.alignments = alignments

    def get_aligned_experiment(self):
        """ Return aligned experiment. """
        return self.experiment


class MultiExperimentAlignment:
    """
    Object for the time-alignment of multiple experiments. The first experiment serves as the reference.

    Attributes:

        experiments (list) - copies of Experiment instances

    """

    def __init__(self, *experiments, **kwargs):
        """
        Instantiate object for the time-alignment of multiple experiments. The first experiment serves as the reference.

        Args:

            experiments (iterable) - data.experiments.Experiment instances

            kwargs: keyword arguments for the alignment

        """
        self.experiments = [deepcopy(experiment) for experiment in experiments]
        self.align_experiments(**kwargs)

    @classmethod
    def _align_experiments(cls, experiments,
                           channel='ch1_normalized', **kwargs):
        """
        Align all experiments with the first experiment.

        Args:

            experiments (iterable) - data.experiments.Experiment instances

            channel (str) - fluorescence channel by which to align experiments

            kwargs: keyword arguments for the alignment

        Returns:

            experiments (list) - list of copied Experiment instances

        """
        reference = cls.aggregate_discs(experiments[0]).select_cell_type('pre')

        ref_initial_offset = reference.data.t.min()

        for experiment in experiments[1:]:
            aggregate_discs = cls.aggregate_discs(experiment).select_cell_type('pre')

            # shift to zero
            initial_offset = aggregate_discs.data.t.min()

            alignment = CellsAlignment(aggregate_discs.data, reference.data, channel=channel, basis='t', **kwargs)
            experiment.apply_lag(alignment.lag+ref_initial_offset-initial_offset)
        return experiments

    def align_experiments(self, **kwargs):
        """
        Align all experiments with first experiment.

        kwargs: keyword arguments for alignment
        """
        experiments = self._align_experiments(self.experiments, **kwargs)
        self.experiments = experiments

    @staticmethod
    def aggregate_discs(experiment):
        """
        Aggregate all discs within an Experiment into single Cells object.

        Args:

            experiment (data.experiments.Experiment)

        Returns:

            cells (data.cells.Cells)

        """
        return reduce(add, experiment.discs.values())

    def get_aligned_experiments(self):
        """
        Return aligned experiments.

        Returns:

            experiments (list) - list of copied Experiment instances

        """
        return self.experiments

