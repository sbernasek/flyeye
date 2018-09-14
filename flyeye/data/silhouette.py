__author__ = 'Sebastian Bernasek'

from os.path import join, abspath, exists
import json
import pandas as pd


class Silhouette:
    """
    Interface to a FlyEye Silhouette file. Upon instantiation, individual cell measurements are aggregated into a data.cells.Cells compatible DataFrame.

    Measurement data must be read on a layer-by-layer basis the first time a Silhouette object is instantiated. Following this initial reading, the aggregated measurement data are serialized and stored within the silhouette file. These serialized measurements may then be accessed directly during future use. The recompile flag indicates whether the serialized measurements should be ignored upon instantiation.

    Attributes:
    path (str) - path to Silhouette file
    df (pd.DataFrame) - cell measurement data
    flip_about_yz (bool) - if True, invert about YZ plane
    flip_about_xy (bool) - if True, invert about XY plane
    """

    def __init__(self, path, recompile=False):
        """
        Instantiate interface to silhouette file.

        Args:
        path (str) - path to silhouette file
        recompile (bool) - if True, recompile measurements from all layers
        """
        self.path = abspath(path)
        self.load(recompile=recompile)

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

    def write_json(self, filename, content):
        """
        Write serialized contents to silhouette file.

        Args:
        filename (str) - filename
        content (dict) - serialized contents
        """
        filepath = join(self.path, filename)
        with open(filepath, 'w') as file:
            json.dumps(content, file)

    def compile_measurements(self):
        """ Compile measurements from all layers (slow access). """
        labels = self.read_labels()
        self.df = self.read_contours(labels)

    def save_measurements(self):
        """ Save serialized measurements for fast access. """
        #serialized_measurements = self.df.to_json(orient='records')
        #self.write_json('measurements.json', serialized_measurements)
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

        # read suggested disc orientation from feed file
        self.load_orientation()

    def load_orientation(self):
        """ Load suggested disc orientation from feed file. """

        # load feed file
        feed = self.read_json('feed.json')

        # read orientation
        self.flip_about_yz = bool(feed['orientation']['flip_about_yz'])
        self.flip_about_xy = bool(feed['orientation']['flip_about_xy'])

    def read_labels(self):
        """
        Load segment labels from silhouette file.

        Returns:
        labels (dict) - {layer_id: {contour_id: label}} entries for each layer
        """

        # load feud file
        feud = self.read_json('feud.json')

        # compile labels for all layers
        labels = {}
        for layer in feud['layers']:

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

        # load feed file
        feed = self.read_json('feed.json')

        # read contours from all layers
        contours = []
        for layer_id in feed['layer_ids']:

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

