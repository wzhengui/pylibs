from setuptools import setup

setup(
    name='pylib_essentials',
    version='0.0.6',
    description='The essentials of the experimental package of pylibs, which require minimum dependencies',
    license='MIT',
    packages=['pylib_essentials'],
    package_data={},
    install_requires=[
      'setuptools',
      'numpy',
      'scipy',
      'pandas',
      'netCDF4',
      'matplotlib>=3.0.0',
      'pyproj>=3.0.0',
    ],
    extras_require={
      'mpi': ['mpi4py>=3.0.0'],
      'shapefile': ['pyshp>=2.1.0'],
      'eof': ['eofs>=1.4.0'],
      'cloudpickle': ['cloudpickle'],
    }
)

