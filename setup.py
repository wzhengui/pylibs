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
        version='0.1.11',  # Ideally should be same as your GitHub release tag varsion
        package_data={'pyUtility': ['prj.npz']},
        description='python libraries and utilities for pre/post-processing SCHISM models',
        long_description='python libraries and utilities for pre/post-processing SCHISM models',
        author='Zhengui Wang',
        author_email='wangzg@vims.edu',
        classifiers=[],
        install_requires=[
    	    'setuptools',
    	    'mpi4py>=3.0.0',
    	    'pandas',
    	    'numpy',
    	    'pyproj>=3.3.0',
    	    'netCDF4>=1.5.8',
    	    'pyshp>=2.1.3',
    	    'matplotlib>=3.0.0',
    	    'scipy'
        ]
    )
