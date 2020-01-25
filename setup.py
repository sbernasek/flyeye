from distutils.core import setup
from setuptools import find_packages

setup(
    name='flyeye',
    version='0.1b',
    author='Sebastian Bernasek',
    author_email='sebastian@u.northwestern.com',
    packages=find_packages(exclude=('tests',)),
    scripts=[],
    url='https://sebastianbernasek.github.io/flyeye/',
    license='MIT',
    description='Analysis package for FlyEye Silhouette data.',
    long_description=open('README.md').read(),
    python_requires='>=3',
    install_requires=[
        "scipy >= 1.1.0",
        "pandas == 0.23.4",
        "matplotlib >= 2.0.0",
        "astroML == 0.3",
        "astroML-addons == 0.2.2",
    ],
    tests_require=['nose'],
    test_suite='nose.collector'
)
