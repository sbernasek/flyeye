from distutils.core import setup
from setuptools import find_packages

setup(
    name='FlyEye',
    version='0.1.0',
    author='Sebastian Bernasek',
    author_email='sebastian@u.northwestern.com',
    packages=find_packages(exclude=('tests',)),
    scripts=[],
    url='https://github.com/sebastianbernasek/flyeye',
    license='MIT',
    description='Analysis package for FlyEye Silhouette data.',
    long_description=open('README.md').read(),
    install_requires=[
        "scipy >= 1.1.0",
        "pandas == 0.23.4",
        "matplotlib >= 2.0.0",
        "astroML == 0.3"
    ],
)
