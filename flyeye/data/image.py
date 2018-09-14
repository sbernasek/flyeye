__author__ = 'Sebastian Bernasek'

from os.path import join
import matplotlib.pyplot as plt
import scipy
import numpy as np
from glob import glob
from matplotlib import patches
from scipy.ndimage import filters
from copy import deepcopy


class Image:
    """
    Object representing an 8-bit image.

    Attributes:
    im (np.ndarray[uint8]) - 8-bit RGB pixel values. Array is shaped (X,Y,3).
    crops (dict) - stored crop boundaries
    """

    def __init__(self, im,
                 fliplr=True,
                 rotation=0,
                 default_crop=None):
        """
        Instantiate image.

        Args:
        im (np.ndarray[uint8]) - 8-bit RGB pixel values
        fliplr (bool) - if True, flip about Y axis
        rotation (float) - rotate image
        default_crop (array like, length 4) - crop boundary supplied as fractions of each axis, eg (xmin, xmax, ymin, ymax) where values are on a 0 to 1 scale
        """

        # store image
        self.im = im

        # flip image about Y axis
        if fliplr == True:
            self.fliplr()

        # rotate image
        self.rotate(rotation)

        # set default crop
        if default_crop is not None:
            xlow, xhigh, ylow, yhigh = default_crop
            xmin, xmax, ymin, ymax = self.get_crop_indices(xlow=xlow, xhigh=xhigh, ylow=ylow, yhigh=yhigh)
            self.im = self.im[ymin:ymax, xmin:xmax, :]

        # instantiate crop library
        self.crops = {}

    def __getitem__(self, color):
        """ Return scalar image of single color channel. """
        return ScalarImage.from_8bit(self.im[:, :, 'rgb'.index(color)])

    def save(fig, filename,
             crop=None,
             dpi=300,
             fmt='png',
             transparent=True):
        """
        Save rendered image.
        """
        fig.savefig(filename, dpi=dpi, format=fmt, transparent=transparent)

    @staticmethod
    def from_tiff(path, **kwargs):
        """
        Load from ndimage readbale image file.

        Args:
        path (str) - image path
        kwargs: keyword arguments for Image instantiation

        Returns:
        im (data.image.Image)
        """
        im = scipy.ndimage.imread(path, flatten=False, mode='RGB')
        return Image(im, **kwargs)

    @staticmethod
    def _fliplr(im):
        """ Flip image array about Y axis. """
        return np.fliplr(im)

    def fliplr(self):
        """ Flip image about Y axis. """
        self.im = self._fliplr(self.im)

    @staticmethod
    def _rotate(im, angle=0):
        """ Roatate image array. """
        return scipy.ndimage.rotate(im, angle)

    def rotate(self, angle):
        """ Roatate image. """
        self.im = self._rotate(self.im, angle=angle)

    @staticmethod
    def _get_crop_indices(im, xlow=0, xhigh=1, ylow=0, yhigh=1):
        """ Compute image index boundaries from fractional bounds. """
        im_shape = im.shape[0:2]
        ymin, ymax = int(im_shape[0]*ylow), int(im_shape[0]*yhigh)
        xmin, xmax = int(im_shape[1]*xlow), int(im_shape[1]*xhigh)
        return xmin, xmax, ymin, ymax

    def get_crop_indices(self, **kwargs):
        """ Compute image index boundaries from fractional bounds. """
        return self._get_crop_indices(self.im, **kwargs)

    def add_crop(self, name='crop', **kwargs):
        """ Add named crop. """
        xmin, xmax, ymin, ymax = self.get_crop_indices(**kwargs)

        self.crops[name] = {'bounds': (xmin, xmax, ymin, ymax),
                           'image': self.im[ymin:ymax, xmin:xmax, :]}

    def build_colorfilter(self, scheme='rgb', channels='rgb', reference='r'):
        """ Build filter for converting between RGB and MG colorschemes. """

        colorfilter = np.identity(3)

        # MG schemes
        if scheme == 'mg':

            if reference == 'r':
                # channel blue to red
                colorfilter[0, 0] = 0
                colorfilter[2, 0] = 1

            else:
                # channel red to blue
                colorfilter[2, 2] = 0
                colorfilter[0, 2] = 1

            if 'm' not in channels:
                colorfilter[2, 0] = 0
                colorfilter[0, 2] = 0
                colorfilter[2, 2] = 0
                colorfilter[0, 0] = 0

            if 'g' not in channels:
                colorfilter[1, 1] = 0

        # RGB schemes
        else:
            for i, c in enumerate(scheme):
                if c not in channels:
                    colorfilter[:, i] = 0

        return colorfilter

    @staticmethod
    def _apply_colorfilter(im, colorfilter):
        """ Apply colorfilter to image array. """
        return np.dot(im, colorfilter).astype(np.uint8)

    def add_rectangle(self, crop, ax,
                      fill=False,
                      linestyle='dashed',
                      color='white',
                      linewidth=3,
                      **kwargs):
        """ Add rectangular patch to axis. """
        bounds = self.crops[crop]['bounds']
        xmin, xmax, ymin, ymax = bounds
        dx, dy = xmax-xmin, ymax-ymin
        box = patches.Rectangle((xmin, ymin), dx, dy,
                                fill=fill, linestyle=linestyle, linewidth=linewidth, color=color, **kwargs)
        ax.add_patch(box)

    @staticmethod
    def _rescale(im, method='range', lbound=35, ubound=99., blur=0.0):
        """ Rescale image colors for maximum contrast. """
        im = deepcopy(im)
        if method == 'range':
            im_scaled = (im - im.min(axis=(0, 1)))/(im.max(axis=(0, 1))-im.min(axis=(0, 1)))

        elif method == 'asymmetric':
            im = filters.gaussian_filter(im, sigma=blur)
            im[im==0] = 1e-10
            log_im = np.log10(im)
            lower, upper = np.percentile(log_im, q=(lbound, ubound))
            log_im[log_im<lower] = lower
            log_im[log_im>upper] = upper
            with np.errstate(divide='ignore'):
                im_scaled = (log_im-lower)/(upper-lower)

        return im_scaled

    def rescale(self, method='range', **kwargs):
        """ Rescale image and all crops for maximum contrast. """
        self.im = self._rescale(self.im, method=method, **kwargs) * 256
        for crop, crop_dict in self.crops.items():
            crop_dict['image'] = self._rescale(crop_dict['image'], method=method, **kwargs) * 256

    @staticmethod
    def _filter(im, f):
        """ Apply function to an image array. """
        if f is None:
            return im
        else:
            return f(im)

    def apply_filter(self, filter_func=None):
        """ Apply filter to image and all crops. """
        self.im = self._filter(self.im, filter_func) * 256
        for crop, crop_dict in self.crops.items():
            crop_dict['image'] = self._filter(crop_dict['image'], filter_func) * 256

    def render(self,
               crop=None,
               ax=None,
               rectangles=None,
               rect_kwargs={},
               scheme='rgb',
               channels='rgb',
               reference='r',
               flipud=False,
               size=(2, 2)):
        """
        Render an image or one of its crops.
        """

        # get colorfilter
        colorfilter = self.build_colorfilter(scheme=scheme, channels=channels, reference=reference)

        # apply filter to selected crop
        if crop is None:
            im = self._apply_colorfilter(self.im, colorfilter)
        else:
            im = self._apply_colorfilter(self.crops[crop]['image'], colorfilter)

        # create axes and render image
        if ax is None:
            fig, ax = plt.subplots(figsize=size)

        # flip upside down
        if flipud:
            im = np.flipud(im)

        ax.imshow(im)
        ax.set_xticks([]), ax.set_yticks([])

        # add rectangles
        if rectangles is not None:
            for rect_crop in rectangles:
                self.add_rectangle(rect_crop, ax=ax, **rect_kwargs)

        plt.tight_layout()
        return ax


class ScalarImage:
    """
    Object representing a floating-point scalar image.

    Attributes:
    im (np.ndarray[np.float64]) - pixel values. Array is shaped (X,Y).
    fig (matplotlib.figure.Figure) - image rendering
    """

    def __init__(self, im):
        """
        Instantiate image.

        Args:
        im (np.ndarray[np.float64]) - 8-bit RGB pixel values
        """
        self.im = im
        self.fig = None

    @staticmethod
    def from_8bit(im):
        """ Instantiate from 8-bit image. """
        return ScalarImage(im.astype(np.float64) / (2**8))

    def __add__(self, image):
        """ Return pixelwise sum. """
        return ScalarImage(self.im+image.im)

    def __sub__(self, image):
        """ Return pixelwise difference. """
        return ScalarImage(self.im-image.im)

    def __mul__(self, image):
        """ Return pixelwise product. """
        return ScalarImage(self.im*image.im)

    def __truediv__(self, image):
        """ Return pixelwise ratio. """
        with np.errstate(divide='ignore'):
            im = np.log2(self.im/image.im)
        return ScalarImage(im)

    def save(self, name,
             directory='./',
             crop=None,
             dpi=300,
             fmt='png',
             transparent=True):
        """
        Save rendered image.

        Args:
        name (str) - file name
        diretory (str) - target directory
        crop (str) - key for crop to be saved, if None save whole figure
        dpi (int) - resolution
        fmt (str) - image format
        transparent (bool) - if True, remove background
        """
        path = join(directory, name+'.'+fmt)
        self.fig.savefig(path, dpi=dpi, format=fmt, transparent=transparent)

    def render(self,
               size=(2, 2),
               cmap=None,
               vmin=0,
               vmax=1,
               ax=None):
        """
        Render image.

        Args:
        size (tuple) - figure size
        cmap (matplotlib.cm.ColorMap) - colormap, defaults to Greys
        vmin, vmax (float) - colorscale bounds
        ax (matplotlib.axes.AxesSubplot) - if None, create figure
        """

        # create figure
        if ax is None:
            fig, ax = plt.subplots(figsize=size)
            self.fig = fig

        # set colormap
        if cmap is None:
            cmap = plt.cm.Greys

        # visualize image
        ax.imshow(self.im, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks([]), ax.set_yticks([])


class ImageStack:
    """
    Object representing an imagestack. Imagestack is compiled from PNG files within a Silhouette file.

    Attributes:
    path (str) - path to silhouette file
    imagestack (np.ndarray[uin8]) - Stack 8-bit RGB pixel values. Array is shaped (N,X,Y,3) where N is the number of layers and X,Y are the dimensions of a single cross section.
    """

    def __init__(self, path):
        """
        Instantiate image stack from silhouette file.

        Args:
        path (str) - path to silhouette file
        """
        self.path = path
        self.imagestack = self.from_silhouette(path)

    @staticmethod
    def from_silhouette(path, fmt='png'):
        """
        Load imagestack from sequentially numbered images of each layer.

        Args:
        path (str) - path to silhouette file
        fmt (str) - image format for all layers

        Returns:
        stack (data.image.ImageStack)
        """

        files = glob(join(path, '*.'+fmt))
        layers_dict = {}
        for fname in files:
            im_filename = fname.split('/')[-1].split('.')[0]
            if '_' not in im_filename:
                layer = int(im_filename)
                image = scipy.ndimage.imread(fname, flatten=False, mode='RGB')
                layers_dict[layer] = image

        # replace missing layers
        for i in np.arange(60):
            if i not in layers_dict.keys():
                layers_dict[i] = layers_dict[i-1]

        imagestack = np.stack(layers_dict[i] for i in np.arange(60))
        return imagestack

    def get_layer(self, layer_id, **kwargs):
        """
        Return image of individual layer.

        Args:
        layer_id (int) - layer number
        kwargs: keyword arguments for Image instantiation

        Returns:
        im (data.image.Image)
        """
        image = Image(self.imagestack[layer_id], **kwargs)
        return image

    def get_max_projection(self, min_layer=0, max_layer=60, **kwargs):
        """
        Return maximum intensity projection across specified layers.

        Args:
        min_layer (int) - lower layer bound
        max_layer (int) - upper layer bound
        kwargs: keyword arguments for Image instantiation

        Returns:
        im (data.image.Image)
        """
        projection = self.imagestack[min_layer: max_layer].max(axis=0)
        image = Image(projection, **kwargs)
        return image

    def get_mean_projection(self, min_layer=0, max_layer=60, **kwargs):
        """
        Return mean intensity projection across specified layers.

        Args:
        min_layer (int) - lower layer bound
        max_layer (int) - upper layer bound
        kwargs: keyword arguments for Image instantiation

        Returns:
        im (data.image.Image)
        """
        projection = self.imagestack[min_layer: max_layer].mean(axis=0)
        image = Image(projection, **kwargs)
        return image
