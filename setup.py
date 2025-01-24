from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

    setup(
        name='pylibs_ocean',
        packages=[ 'pylibs/src', 'pylibs/scripts'],
        package_dir={'pylibs/src':'src','pylibs/scripts':'scripts'},
        py_modules=['pylib'],
        version='1.2.1',  # Ideally should be same as your GitHub release tag varsion
        package_data={'pylibs/scripts': ['prj.npz','sflux_template.npz','Harmonic_Analysis/*','schismcheck','schismview']},
        description='python tools for ocean reserach',
        long_description='python libraries and utilities for data processing including the pre/post-processing about SCHISM models',
        author='Zhengui Wang',
        author_email='wzhengui@gmail.com',
        url='https://github.com/wzhengui/pylibs',
        classifiers=[],
        install_requires=[
    	    'setuptools',
    	    'pandas',
    	    'numpy',
    	    'netCDF4>=1.5.8',
    	    'matplotlib>=3.0.0',
    	    'scipy'
        ],
        extras_require={
          'mpi': ['mpi4py>=3.0.0'],
          'shapefile': ['pyshp>=2.1.0'],
          'projection': ['pyproj>=3.0.0'],
          'eof': ['eofs>=1.4.0'],
          'cloudpickle': ['cloudpickle==2.2.1'],
          'tiff': ['tifffile==2022.5.4','imagecodecs==2024.6.1'],
          'urllib': ['urllib3==2.2.1'],
          'basemap': ['basemap==1.4.1','basemap-data-hires==1.3.2'],
          'all': ['mpi4py>=3.0.0','pyshp>=2.1.0','pyproj>=3.0.0','eofs>=1.4.0','cloudpickle==2.2.1',
                  'tifffile==2022.5.4','imagecodecs==2024.6.1','urllib3==2.2.1','basemap==1.4.1','basemap-data-hires==1.3.2']
        }
    )
