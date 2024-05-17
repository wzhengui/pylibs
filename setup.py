from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

    setup(
        name='pylibs_ocean',
        packages=[
          'pyUtility',
    	  'pyScripts',
        ],
        py_modules=['pylib'],
        version='1.0.1',  # Ideally should be same as your GitHub release tag varsion
        package_data={'pyScripts': ['prj.npz','sflux_template.npz','Harmonic_Analysis/*','schismcheck','schismview']},
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
        }
    )
