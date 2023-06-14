from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

    setup(
        name='pylibs',
        packages=[
          'pyUtility',
    	    'pyScripts',
        ],
        py_modules=['pylib'],
        version='0.1.12',  # Ideally should be same as your GitHub release tag varsion
        package_data={'pyScripts': ['prj.npz','sflux_template.npz','Harmonic_Analysis']},
        description='python libraries and utilities for pre/post-processing SCHISM models',
        long_description='python libraries and utilities for pre/post-processing SCHISM models',
        author='Zhengui Wang',
        author_email='wangzg@vims.edu',
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
        }
    )
