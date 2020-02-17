from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps
from matplotlib import colors
import matplotlib.gridspec as gs

from ..dynamics.averages import get_rolling_mean
from ..dynamics.visualization import plot_mean


class Triangulation:
    """
    Object for estimating the median distance between adjacent columns of R8 cells within an individual eye disc. Distance estimate is obtained by constructing a Delaunay graph connecting all annotated R8 neurons, filtering the edges by length and angle relative to the horizontal axis, then evaluating the median x-component of remaining edges.

    The median inter-column distance is multiplied by the estimated MF velocity (0.5 columns/hr) to generate a distance-to-time scaling factor.

    Attributes:

        params (dict) - triangulation parameters, {name: value}

        xycoords (np.ndarray) - R8 cell positions

        delaunay (scipy.spatial.tri) - Delaunay triangulation

        distances (np.ndarray) - distances between adjacent R8 cells

        edges (np.ndarray) - edge vertices

        hours_per_pixel (float) - distance to time scaling factor

        disc (data.discs.Disc)

    """

    def __init__(self, disc,
                 furrow_velocity=2,
                 threshold=None,
                 min_angle=30,
                 max_angle=60,
                 include_x=True,
                 include_y=False):
        """
        Instantiate object for estimating the median distance between adjacent columns of R8 cells.

        Args:

            disc (data.discs.Disc)

            furrow_velocity (float) - furrow inverse-velocity (hours per column)

            threshold (float) - max. quantile of included distances, 0 to 100

            min_angle, max_angle (float) - min/max angle of included edges

            include_x (bool) - if True, include x-distance

            include_y (bool) - if True, include y-distance

        """

        self.triangulate(disc, furrow_velocity, threshold, min_angle=min_angle, max_angle=max_angle, include_x=include_x, include_y=include_y)

        self.params = {'furrow_velocity': furrow_velocity,
                        'threshold': threshold,
                        'min_angle': min_angle,
                        'max_angle': max_angle,
                        'include_x': include_x,
                        'include_y': include_y}

    def __call__(self, disc):
        """
        Apply distance to time scaling to a disc.

        Args:

            disc (data.discs.Disc)

        Returns:

            disc (data.discs.Disc) - disc with estimated developmental times

        """
        disc = self._apply_time_scaling(disc, self.hours_per_pixel)
        return disc

    def get_disc(self):
        """ Return disc. """
        return self.disc

    @staticmethod
    def _apply_time_scaling(disc, hours_per_pixel):
        """
        Update developmental times.

        Args:

            disc (data.discs.Disc)

            hours_per_pixel (float) - distance to time scaling factor

        Returns:

            disc (data.discs.Disc) - disc with estimated developmental times

        """
        disc.data['t'] = disc.data.centroid_x * hours_per_pixel
        return disc

    @staticmethod
    def _get_delaunay(xycoords):
        """ Return Delaunay triangulation for xy points. """
        return Delaunay(xycoords)

    @staticmethod
    def _get_edges(delaunay,
                   min_angle=30,
                   max_angle=60,
                   include_x=True,
                   include_y=False):
        """ Get distances between adjacent R8 neurons. """

        # get indices of vertex neighbors
        indices, indptr = delaunay.vertex_neighbor_vertices

        # iterate through all R8 cells
        distances, edges = [], []
        for k, (x1, y1) in enumerate(delaunay.points):

            # get all neighbors of current R8
            neighbors = delaunay.points[indptr[indices[k]:indices[k+1]]]

            # if neighbor is within a 30-60 degree angle from horizontal, include x-coordinate of its distance
            for x2, y2 in neighbors:
                theta = np.arctan(abs((y2-y1)/(x2-x1)))

                # if neighbor is within the specified angle from horizontal, include distance
                if theta >= min_angle*np.pi/180 and theta <= max_angle*np.pi/180:
                    distances.append(np.sqrt(((x2-x1)**2)*include_x + ((y2-y1)**2)*include_y))
                    edges.append([[x1, x2], [y1, y2]])

        return np.array(distances), np.array(edges)

    @staticmethod
    def _filter_edges_by_length(distances, edges, threshold=1.75):
        """
        Filter edges by length. Length threshold is computed as a multiple of the median edge length.

        Args:

            distances (np.ndarray) - edge lengths

            edges (np.ndarray) - edges

            threshold (float) - maximum multiple of median length

        Returns:

            distances (np.ndarray) - filtered edge lengths

            edges (np.ndarray) - filtered edges

        """
        indices = np.where(distances < threshold*np.median(distances))[0]
        return distances[indices], edges[indices]

    def triangulate(self, disc,
                    furrow_velocity=2,
                    threshold=None,
                    min_angle=30,
                    max_angle=60,
                    include_x=True,
                    include_y=False):
        """
        Run triangulation.

        Args:

            disc (data.discs.Disc)

            furrow_velocity (float) - furrow inverse-velocity (hr/column)

            threshold (float) - max quantile of included distances, 0 to 100

            min_angle, max_angle (float) - min/max angle of included edges

            include_x (bool) - if True, include x-distance

            include_y (bool) - if True, include y-distance

        """

        # get coordinates
        xycoords = disc.data[disc.data.label=='r8'][['centroid_x', 'centroid_y']].values
        self.xycoords = xycoords

        # get triangulation
        self.delaunay = self._get_delaunay(xycoords)

        # get edges
        self.distances, self.edges = self._get_edges(self.delaunay, min_angle=min_angle, max_angle=max_angle, include_x=include_x, include_y=include_y)

        # filter edges
        if threshold is not None:
            self.distances, self.edges = self._filter_edges_by_length(distances=self.distances, edges=self.edges, threshold=threshold)

        # compute mean distance to time scaling
        self.hours_per_pixel = furrow_velocity/np.mean(self.distances)

        # update time vector
        self.disc = self._apply_time_scaling(disc, self.hours_per_pixel)

    @staticmethod
    def _get_log2_fold_change(values):
        return np.log2(values/values.mean())

    @classmethod
    def _add_edges_to_plot(cls,
                           distances,
                           edges,
                           ax,
                           hours_per_pixel,
                           cmap=cmaps.coolwarm):
        """ Add delaunay edges to existing axes. """

        # get scores for colormap
        scores = cls._get_log2_fold_change(distances)

        # plot lines
        for edge, score in zip(edges, scores):
            times = [x * hours_per_pixel for x in edge[0]]
            ax.plot(times, edge[1], '-', linewidth=2, alpha=1, color=cmap.to_rgba(score), zorder=1)

        ax.set_yticks([])
        ax.set_xlabel('time (hr)')

    def add_edges_to_plot(self, ax, cmap=cmaps.coolwarm):
        """ Add delaunay edges to existing axes. """
        cmap = self.get_colormap(cmap)
        self._add_edges_to_plot(distances=self.distances, edges=self.edges, ax=ax,
                                hours_per_pixel=self.hours_per_pixel, cmap=cmap)
        ax.set_ylim(self.xycoords[:, 1].min(), self.xycoords[:, 1].max())

    @staticmethod
    def get_colormap(cmap=cmaps.coolwarm):
        norm = colors.Normalize(vmin=-1, vmax=1)
        colormap = cmaps.ScalarMappable(norm=norm, cmap=cmap)
        return colormap

    @classmethod
    def _plot_histogram(cls, values, ax=None, dist_type='x'):
        """ Plot colored histogram for a set of values. """

        # get colormap
        cmap = cls.get_colormap(cmap=cmaps.coolwarm)

        # histogram values
        counts, bin_edges = np.histogram(values)
        bin_centers = [(edge+bin_edges[i+1])/2 for i, edge in enumerate(bin_edges[:-1])]
        scores = np.log2(bin_centers / np.mean(values))
        patches = ax.bar(bin_edges[:-1], counts, width=(bin_edges[1]-bin_edges[0]), color=[cmap.to_rgba(score) for score in scores])

        # format plot
        ax.set_yticks([])
        ax.set_xlim(0, 100)
        ax.text(95, .9*max(counts), s=' Mean: {:0.1f} px'.format(np.mean(values)), ha='right', va='top')
        ax.text(95, .9*max(counts), s='\n Med.: {:0.1f} px'.format(np.median(values)), ha='right', va='top')
        ax.text(95, .9*max(counts), s='\n\n N = {:d}'.format(int(len(values)/2)), ha='right', va='top')
        _ = ax.set_xlabel(dist_type+'-distance')
        _ = ax.set_ylabel('edges')

        return counts, bin_edges, patches

    def plot_histogram(self, ax):
        """ Histogram inter-R8 distances. """
        dist_type = self.params['include_x']*'x'+self.params['include_y']*'y'
        self._plot_histogram(self.distances, ax=ax, dist_type=dist_type)

    @staticmethod
    def _plot_expression(ax, disc, channel, hours_per_pixel,
                         window_size=100,
                         color='black',
                         alpha=1):
        """ Plot expression trajectories. """

        cells = disc.select_cell_type('pre')
        cells.plot_dynamics(channel, ax=ax,
                            line_kw={'color': color, 'alpha': alpha, 'lw': 2})

        # format x axis
        ax.set_ylim(0, 2)
        ax.set_yticks([])
        ax.set_xlim(-15, 55)
        ax.set_xticks(np.arange(-10, 60, step=10))
        ax.set_xticklabels([str(int(round(label, 0))) for label in ax.get_xticks()])

    def plot_expression(self, ax, channel, **kwargs):
        """ Plot expression trajectory. """
        self._plot_expression(ax, self.disc, channel, self.hours_per_pixel, **kwargs)

    def overlay_epression(self, ax, channel, **kwargs):
        """ Plot expression trajectory on twin y-axis. """
        ax_alt = ax.twinx()
        self.plot_expression(ax_alt, channel, **kwargs)

    def show(self,
             gs_parent=None,
             include_expression=True,
             channel=None,
             is_subplot=False, **kwargs):
        """
        Plot inter-R8 distance distribution, Delaunay triangulation, and expression.
        """

        # retriangulate
        self.triangulate(self.disc, **self.params)

        # create axes
        if gs_parent is None:
            fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(6, 2))
        else:
            gs_child = gs.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1.5], subplot_spec=gs_parent, hspace=0)
            ax0 = plt.subplot(gs_child[0])
            ax1 = plt.subplot(gs_child[1])
            is_subplot = True

        # add colorbar
        self.add_colorbar(ax1, is_subplot=is_subplot)

        # plot edges, expression, and histogram
        self.add_edges_to_plot(ax1, cmap=cmaps.coolwarm)
        if include_expression and channel is not None:
            self.overlay_epression(ax1, channel, **kwargs)
        self.plot_histogram(ax=ax0)

        plt.tight_layout()

        return ax0, ax1

    @staticmethod
    def add_colorbar(ax, is_subplot=False, fraction=0.2):
        mappable = ax.scatter([1e6, 1e6], [0, 0], alpha=1, c=[-1, 1], cmap=cmaps.coolwarm)
        cbar = plt.colorbar(mappable=mappable, ax=ax, fraction=fraction)

        # simplify colorbar for subplots
        if is_subplot is False:
            cbar.set_ticks([-1, 0.35])
            cbar.ax.tick_params(length=0)
            cbar.ax.set_yticklabels(['compressed', 'stretched'], rotation='vertical', fontsize=6, ha='left', va='bottom')
            cbar.set_label('Log2(F.C. wrt mean)', fontsize=8)

        else:
            cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
            cbar.ax.tick_params(length=0, labelsize=8)
            cbar.set_label('Log2(F.C. wrt mean)', fontsize=8)


class ExperimentTriangulation:
    """
    Object for estimating the median distance between adjacent columns of R8 cells for each disc within an experiment. Distance estimate is obtained by constructing a Delaunay graph connecting all annotated R8 neurons, filtering the edges by length and angle relative to the horizontal axis, then evaluating the median x-component of remaining edges.

    The median inter-column distance is multiplied by the estimated MF velocity (0.5 columns/hr) to generate a distance-to-time scaling factor.

    Attributes:

        experiment (data.experiments.Experiment)

        tri (dict) - {disc ID: Triangulation} pairs

    """

    def __init__(self, experiment, **kwargs):
        """
        Instantiate triangulation objects for all discs in an experiment.discs

        Args:

            experiment (data.experiments.Experiment)

            kwargs: triangulation keyword arguments

        """
        discs = experiment.discs
        self.tri = self._get_triangulations(discs, **kwargs)
        self.experiment = self._apply_triangulations(experiment, self.tri)

    def __call__(self):
        """ Return experiment with updated developmental times. """
        return self.experiment

    @staticmethod
    def _apply_triangulations(experiment, tri):
        """
        Return experiment with updated developmental times.

        Args:

            experiment (data.experiments.Experiment) - experiment to be updated

            tri (dict) - {disc ID: Triangulation} pairs

        Returns:

            experiment (data.experiments.Experiment) - updated experiment

        """
        experiment.discs = {i: t.get_disc() for i, t in tri.items()}
        return experiment

    @staticmethod
    def _get_triangulations(discs, **kwargs):
        """
        Return dictionary of Triangulation objects.

        Args:

            discs (dict) - {disc ID: Disc} pairs

            kwargs: keyword arguments for Triangulation

        Returns:

            tri (dict) - {disc ID: Triangulation} pairs

        """
        return {i: Triangulation(disc, **kwargs) for i, disc in discs.items()}

    def plot_expression(self, ax, channel,
                        color='black',
                        **kwargs):
        """ Plot expression for all triangulations. """
        for triangulation in self.triangulations.values():
            triangulation.plot_expression(ax, channel=channel, color=color, **kwargs)

    def show_triangulations(self):
        """ Visualize all triangulations. """
        gs_parent = gs.GridSpec(nrows=len(self.triangulations), ncols=1)
        fig = plt.figure(figsize=(6, 1.5*len(self.triangulations)))
        for i, (triangulation, gs0) in enumerate(zip(self.triangulations.values(), gs_parent)):
            ax0, ax1 = triangulation.show(gs_parent=gs0)
            if i != len(self.triangulations)-1:
                ax0.set_xticks([]), ax0.set_xlabel('')
                ax1.set_xticks([]), ax1.set_xlabel('')
        plt.tight_layout()
        return fig

    def show_alignment(self, channel,
                       xoffsets=None,
                       ax=None,
                       scatter=False,
                       legend=True,
                       window_size=100,
                       ma_type='sliding',
                       color_wheel='cmyk',
                       figsize=(4, 3)):
        """ Plot alignment of all discs. """

        # assume zero offset for each channel
        if xoffsets is None:
            xoffsets = np.zeros(len(self.triangulations))

        # create axes if none provided
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        # add reference lines
        ax.plot([0, 0], [0, 3], '--k')
        ax.plot([-15, 55], [0.25, 0.25], '--k')

        handles, labels = [], []
        for i, triangulation in self.triangulations.items():

            # get line color for disc
            color = color_wheel[i % len(color_wheel)]

            # get cells
            disc_cells_data = triangulation.data[np.logical_or(triangulation.data.label=='pre', triangulation.data.label=='p')]
            if scatter is True:
                ax.plot(disc_cells_data.t + xoffsets[i], disc_cells_data[channel], '.', alpha=0.1, color=color)

            # add line average
            line = plot_mean(disc_cells_data.t + xoffsets[i], (disc_cells_data[channel]),
                                   ax=ax, label='Disc {:d}'.format(i), ma_type=ma_type, window_size=window_size,
                                   line_color=color, line_width=3, line_alpha=0.5)
            handles.append(line[0]), labels.append('Disc {:d}'.format(i))

        if legend is True:
            ax.legend(handles=handles, labels=labels, loc=0, frameon=False)

        ax.tick_params(labelsize=16)
        ax.set_ylim(0, 2)
        ax.set_ylabel('Fluorescence (a.u.)', fontsize=16)
        ax.set_xlim(-5, 55)
        ax.set_xlabel('Time (hr)', fontsize=16)
