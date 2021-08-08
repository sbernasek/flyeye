from distutils.core import setup
from setuptools import find_packages
from os import path

# read the contents of your README file
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='flyeye',
    version='0.6',
    author='Sebastian Bernasek',
    author_email='sebastian@u.northwestern.com',
    packages=find_packages(exclude=('tests',)),
    scripts=[],
    url='https://sbernasek.github.io/flyeye/',
    license='MIT',
    description='Analysis package for FlyEye Silhouette data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires='>=3',
    install_requires=[
        "scipy >= 1.1.0",
        "pandas >= 1.0",
        "matplotlib >= 2.0.0"],
    extras_require = {
        'frequency_detection': ["astroML == 0.3", "astroML-addons == 0.2.2"],
        },
    tests_require=['nose'],
    test_suite='nose.collector'
)
