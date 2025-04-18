from setuptools import setup

setup(
    name='pylib_experimental',
    version='0.0.1',
    description='Experimental functions',
    license='MIT',
    packages=['pylib_experimental'],
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

