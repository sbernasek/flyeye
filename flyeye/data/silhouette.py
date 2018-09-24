__author__ = 'Sebastian Bernasek'

from os.path import join, abspath, exists
import json
import pandas as pd


class Silhouette:
    """
    Interface to a FlyEye Silhouette file.

    Attributes:

        path (str) - path to Silhouette file

        feed (dict) - feed file containing layer IDs

        feud (dict) - feud file containing cell type labels

    """
    def __init__(self, path):
        """
        Instantiate interface to silhouette file.

        Args:

            path (str) - path to silhouette file

        """
        self.path = abspath(path)

        # load feed and feud files
        self.feed = self.read_json('feed.json')
        self.feud = self.read_json('feud.json')

        # read suggested disc orientation from feed file
        self.load_orientation()

    def read_json(self, filename):
        """
        Read contents of specified JSON file.

        Args:

            filename (str) - filename

        Returns:

            out (dict) - file contents

        """
        filepath = join(self.path, filename)
        with open(filepath, 'r') as f:
            out = json.load(f)
        return out

    def load_orientation(self):
        """ Load suggested disc orientation from feed file. """

        # read orientation
        self.flip_about_yz = bool(self.feed['orientation']['flip_about_yz'])
        self.flip_about_xy = bool(self.feed['orientation']['flip_about_xy'])


class SilhouetteData(Silhouette):
    """
    Interface to data within a FlyEye Silhouette file.

    Upon instantiation, individual cell measurements are aggregated into a data.cells.Cells compatible DataFrame.

    Measurement data must be read on a layer-by-layer basis the first time a Silhouette object is instantiated. Following this initial reading, the aggregated measurement data are serialized and stored within the silhouette file. These serialized measurements may then be accessed directly during future use. The recompile flag indicates whether the serialized measurements should be ignored upon instantiation.

    Attributes:

        df (pd.DataFrame) - cell measurement data

        flip_about_yz (bool) - if True, invert about YZ plane

        flip_about_xy (bool) - if True, invert about XY plane

    Inherited attributes:

        path (str) - path to Silhouette file

        feed (dict) - feed file containing layer IDs

        feud (dict) - feud file containing cell type labels

    """

    def __init__(self, path, recompile=False):
        """
        Instantiate interface to silhouette file data.

        Args:

            path (str) - path to silhouette file

            recompile (bool) - if True, recompile measurements from all layers

        """
        super().__init__(path)
        self.load(recompile=recompile)

    def compile_measurements(self):
        """ Compile measurements from all layers (slow access). """
        labels = self.read_labels()
        self.df = self.read_contours(labels)

    def save_measurements(self):
        """ Save serialized measurements for fast access. """
        self.df.to_json(join(self.path, 'measurements.json'))

    def load_measurements(self):
        """ Load serialized measurements (fast access). """
        self.df = pd.read_json(join(self.path, 'measurements.json'))

    def load(self, recompile=False):
        """
        Read all contour and orientation data from silhouette file.

        Args:

            recompile (bool) - if True, recompile measurements from all layers

        """

        # check whether measurements are available
        measurements_available = exists(join(self.path, 'measurements.json'))

        # load available measurements if recompile flag is false
        if measurements_available and not recompile:
            self.load_measurements()

        # otherwise, recompile and save measurements
        else:
            self.compile_measurements()
            self.save_measurements()

    def read_labels(self):
        """
        Load segment labels from silhouette file.

        Returns:

            labels (dict) - {layer_id: {contour_id: label}} entries for each layer

        """

        # compile labels for all layers
        labels = {}
        for layer in self.feud['layers']:

            # compile {contour_id: contour_label} dictionary for current layer
            annotations = {}
            for contour in layer['contours']:
                label = contour.get('label', None)
                if label is None or label.strip() == '':
                    continue
                annotations[contour['id']] = contour['label']

            # store labels for current layer
            labels[layer['id']] = annotations

        return labels

    @staticmethod
    def parse_contour(contour):
        """
        Convert contour to list format.

        Args:

            contour (dict) - contour from silhouette file

        Returns:

            ctr_list (list) - values in data.cells.Cells compatible list format

        """

        # extract values
        centroid = contour['centroid']
        pixel_count = contour['pixel_count']
        color_avg = contour['color_avg']
        color_std = contour['color_std']

        # reorganize values
        ctr_list = [centroid[0], centroid[1],
                    pixel_count,
                    color_avg['g'], color_std['g'],
                    color_avg['r'], color_std['r'],
                    color_avg['b'], color_std['b']]

        return ctr_list

    def read_contours(self, all_labels={}):
        """
        Read contours from silhouette file.

        Args:

            all_labels (dict) - {layer_id: {contour_id: label}} for each layer

        Returns:

            df (pd.DataFrame) - data.cells.Cells compatible dataframe of contours

        """

        # read contours from all layers
        contours = []
        for layer_id in self.feed['layer_ids']:

            # load labels for current layer
            labels = all_labels.get(layer_id, None)

            # skip layers without any labels
            if labels is None:
                continue

            # read layer from silhouette file
            layer = self.read_json('{:d}.json'.format(layer_id))

            # read all contours within layer
            for contour in layer['contours']:

                # get label for current contour
                label = labels.get(contour['id'], None)

                # skip unlabeled contours
                if label is None:
                    continue

                # convert to list format
                ctr_list = self.parse_contour(contour)
                ctr_list.extend([layer_id, label])

                # store contours from current layer
                contours.append(ctr_list)

        # compile dataframe
        columns = ['centroid_x',
                   'centroid_y',
                   'pixel_count',
                   'green', 'green_std',
                   'red', 'red_std',
                   'blue', 'blue_std',
                   'layer', 'label']
        df = pd.DataFrame(contours, columns=columns)

        return df

