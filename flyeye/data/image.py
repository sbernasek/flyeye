from os.path import join
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from matplotlib import patches
from matplotlib.colors import Normalize
import matplotlib.image as mpimg
from scipy.ndimage import filters, rotate
from scipy.ndimage import grey_closing, median_filter
from copy import deepcopy

from .silhouette import Silhouette


class Smoothing:
    """
    Smoothing operation applied to images.

    Default behavior entails three iterations of grey closing followed by a 3-pixel wide median filter.

    Attributes:

        niters (int) - number of grey closing iterations

        size (int) - pixel size of median filter in XY dimensions

    """

    def __init__(self, niters=3, size=3):
        """
        Instantiate smoothing filter.

        Args:

            niters (int) - number of grey closing iterations

            size (int) - pixel size of median filter in XY dimensions

        """
        self.niters = niters
        self.size = size

    def __call__(self, im):
        """
        Apply smoothing to image.

        Args:

            im (np.ndarray(np.float32)) - pixel values. Shape may be (X,Y), (X,Y,3), or (N,X,Y,3) depending on image type.

        Returns:

            im (np.ndarray(np.float32))

        """
        im = self.grey_closing(im, niters=self.niters)
        im = self.median_filter(im, size=self.size)
        return im

    @staticmethod
    def grey_closing(im, niters=3):
        """
        Apply grey closing operation to image.

        Args:

            im (np.ndarray(np.float32)) - pixel values. Shape may be (X,Y), (X,Y,3), or (N,X,Y,3) depending on image type.

            niters (int) - number of grey closing iterations

        Returns

            im (np.ndarray(np.float32))

        """

        # get filter size
        shape_to_size = {2: 3, 3: (3,3,1), 4: (1,3,3,1)}
        filter_size = shape_to_size[len(im.shape)]

        # apply grey closing
        for _ in range(niters):
            im = grey_closing(im, size=filter_size)

        return im

    @staticmethod
    def median_filter(im, size=3):
        """
        Apply grey closing operation to image.

        Args:

            im (np.ndarray(np.float32)) - pixel values. Shape may be (X,Y), (X,Y,3), or (N,X,Y,3) depending on image type.

            size (int) - pixel size of median filter in XY dimensions

        Returns

            im (np.ndarray(np.float32))

        """

        # get filter size
        shape_to_size = {2: size, 3: (size,size,1), 4: (1,size,size,1)}
        filter_size = shape_to_size[len(im.shape)]

        return median_filter(im, size=filter_size)


class ScalarField:
    """
    Object representing a floating-point 2D scalar field that can take on positive or negative values.

    Attributes:

        im (np.ndarray[float32]) - scalar values, (X,Y)

        fig (matplotlib.figure.Figure) - image rendering

    """

    def __init__(self, im):
        """
        Instantiate scalar field.

        Args:

            im (np.ndarray[float32]) - scalar values, (X,Y)

        """
        self.im = im
        self.fig = None

    def __add__(self, image):
        """ Return pixelwise sum. """
        return ScalarField(self.im+image.im)

    def __sub__(self, image):
        """ Return pixelwise difference. """
        return ScalarField(self.im-image.im)

    def __mul__(self, image):
        """ Return pixelwise product. """
        return ScalarField(self.im*image.im)

    def __truediv__(self, image):
        """ Return pixelwise ratio. """
        with np.errstate(divide='ignore'):
            im = np.log2(self.im/image.im)
        return ScalarField(im)

    @property
    def max(self):
        """ Maximum pixel value. """
        return self.im.max()

    @property
    def min(self):
        """ Minimum pixel value. """
        return self.im.min()

    @property
    def abs_max(self):
        """ Maximum absolute pixel value. """
        return np.abs([self.min, self.max]).max()
    
    @staticmethod
    def from_8bit(im):
        """ Instantiate from 8-bit image. """
        return ScalarField(im.astype(np.float64) / 255)

    @staticmethod
    def uint8_to_float32(im_uint8):
        """
        Convert 8-bit array to floating point values on a 0-1 interval.

        Args:

            im_uint8 (np.ndarray[uint8])

        Returns:

            im_float32 (np.ndarray[float32])

        """
        return im_uint8.astype(np.float64) / 255

    @staticmethod
    def _fliplr(im):
        """ Flip image array about Y axis. """
        return np.fliplr(im)

    def fliplr(self):
        """ Flip image about Y axis. """
        self.im = self._fliplr(self.im)

    @staticmethod
    def _rotate(im, angle=0):
        """ Rotate image array. """
        return rotate(im, angle)

    def rotate(self, angle):
        """ Apply rotation to image. """
        self.im = self._rotate(self.im, angle=angle)

    def apply_colormap(self, cmap, vmin=0, vmax=2):
        """
        Render image to RGB format using provided colormap.

        Args:

            cmap (matplotlib.colors.ColorMap)

            vmin, vmax (int) - color map scale

        Returns:

            im (np.ndarray[float]) - RGB image

        """
        norm = Normalize(vmin, vmax)
        return cmap(norm(self.im))

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

    def apply_filter(self, func):
        """
        Apply filter to image.

        Args:

            func (function) - function operating on image array

        """
        self.im = func(self.im)

    def smooth(self, **kwargs):
        """
        Smooth image.

        kwargs: keyword arguments for smoothing operation
        """

        # instantiate smoothing operation
        smoothing = Smoothing(**kwargs)

        # apply smoothing operation
        self.apply_filter(smoothing)

    def save(self, name,
             directory='./',
             dpi=300,
             fmt='png',
             transparent=True):
        """
        Save rendered image.

        Args:

            name (str) - file name

            directory (str) - target directory

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


class Image(ScalarField):
    """
    Object contains a floating point three channel image.

    Attributes:
    im (np.ndarray[float32]) - RGB pixel values, (X,Y,3)
    """

    def __init__(self, im):
        """
        Instantiate image.

        Args:

            im (np.ndarray[float32]) - RGB pixel values, (X,Y,3)

        """
        super().__init__(im)

    def __getitem__(self, channel):
        """
        Return scalar field of a single color channel.

        Args:

            channel (str) - channel to extract, one of 'r', 'g', or 'b'

        Returns:

            image (ScalarField)

        """
        return ScalarField(self.im[:, :, 'rgb'.index(channel)])

    @staticmethod
    def load_image(path):
        """ Load image from <path>. """

        img_format = path.rsplit('.')[-1]
        if img_format.lower() != 'png':
            raise UserWarning('PIL library is required for non-png formats.')

        # read image (non-PNG formats require pillow library)
        im = mpimg.imread(path)

        return im

    @staticmethod
    def from_tiff(path, **kwargs):
        """
        Read image file using matplotlib.image.imread.

        All formats except PNG require pillow or PIL. PNG images are converted from 8-bit to floating point format by default, while other image formats must be manually converted.

        Args:

            path (str) - image path

            kwargs: keyword arguments for Image instantiation

        Returns:

            im (data.image.Image)

        """

        im = Image.load_image(path)

        # convert to floating point format (float32)
        if im.dtype == np.uint8:
            im = uint8_to_float32(im)

        return Image(im, **kwargs)

    def build_colorfilter(self,
                          scheme='rgb',
                          channels='rgb',
                          reference='r'):
        """
        Build filter for converting between RGB and MG colorschemes.

        Args:

            scheme (str) - colorscheme, 'rgb' or 'mg'

            channels (str) - included channels, one or more of 'rgbmg'

            reference (str) - reference channel (used for rgb to mg conversion)

        Returns:

            colorfilter (np.ndarray[np.float32]) - multiplicative colorfilter

        """

        colorfilter = np.identity(3, dtype=np.float32)

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
        """
        Apply colorfilter to image array.

        Args:

            im (np.ndarray[float32]) - RGB pixel values, (X,Y,3)

            colorfilter (np.ndarray[np.float32]) - multiplicative colorfilter

        Returns:

            im_filtered (np.ndarray[float32]) - new color values, (X,Y,3)

        """
        return np.dot(im, colorfilter)

    def render(self,
               ax=None,
               size=(2, 2),
               scheme='rgb',
               channels='rgb',
               reference='r'):
        """
        Render image.

        Args:

            ax (matplotlib.axes.AxesSubplot) - if None, create figure

            size (tuple) - image size, used if no axis provided

            scheme (str) - colorscheme, 'rgb' or 'mg'

            channels (str) - included channels, one or more of 'rgbmg'

            reference (str) - reference channel (used for rgb to mg conversion)

        Returns:

            ax (matplotlib.axes.AxesSubplot)

        """

        # apply colorfilter
        colorfilter = self.build_colorfilter(scheme, channels, reference)
        im = self._apply_colorfilter(self.im, colorfilter)

        # create axes and render image
        if ax is None:
            fig, ax = plt.subplots(figsize=size)

        # render image
        ax.imshow(im)
        ax.set_xticks([]), ax.set_yticks([])

        return ax


class ImageStack(Image):
    """
    Object representing a 3D stack of RGB image layers. Stack is compiled from sequentially numbered images of each layer within a silhouette file.

    Attributes:

        im (np.ndarray[float32]) - RGB pixel values. Array is shaped (N,X,Y,3) where N is the number of layers and (X,Y) are the dimensions of a single cross section.
    """

    def __init__(self, im):
        """
        Instantiate image stack.

        Args:

            im (np.ndarray[float32]) - RGB pixel values, (N,X,Y,3)

        """
        super().__init__(im)

    def __getitem__(self, layer_id):
        """ Returns RGB image of individual layer. """
        return self.get_layer(layer_id)

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

        # load feed file containing layer indices
        silhouette = Silhouette(path)

        # iterate across layer
        ims = []
        for layer in sorted(silhouette.feed['layer_ids']):

            im_path = join(path, '{:d}.{:s}'.format(layer, fmt))

            # read image
            im = Image.load_image(im_path)

            # convert to floating point format (float32)
            if im.dtype == np.uint8:
                im = uint8_to_float32(im)

            # append to imagestack
            ims.append(im)

        # compile image stack
        im = np.stack(ims)

        # flip about yz
        if silhouette.flip_about_yz:
            im = im[:, ::-1, :, :]

        # flip about xy
        if silhouette.flip_about_xy:
            im = im[::-1, :, :, :]

        return ImageStack(im)

    def get_layer(self, layer_id):
        """
        Return RGB image of individual layer.

        Args:

            layer_id (int) - layer number

        Returns:

            image (data.image.Image)

        """
        return Image(self.im[layer_id])

    def project_max(self, min_layer=0, max_layer=-1):
        """
        Return maximum intensity projection across specified layers.

        Args:

            min_layer (int) - lower layer bound

            max_layer (int) - upper layer bound

        Returns:

            image (data.image.Image)

        """
        projection = self.im[min_layer: max_layer].max(axis=0)
        return Image(projection)

    def project_mean(self, min_layer=0, max_layer=-1):
        """
        Return mean intensity projection across specified layers.

        Args:

            min_layer (int) - lower layer bound

            max_layer (int) - upper layer bound

        Returns:

            image (data.image.Image)

        """
        projection = self.im[min_layer: max_layer].mean(axis=0)
        return Image(projection)
