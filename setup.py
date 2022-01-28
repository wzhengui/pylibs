import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name='pylibs',
    packages=[
        'pylib_Utility',
        'pylib_Scripts',
    ],
    py_modules=['pylib'],
    version='0.1.11',  # Ideally should be same as your GitHub release tag varsion
    package_data={'pylib_Utility': ['prj.npz']},
    description='python libraries and utilities for pre/post-processing SCHISM models',
    long_description='python libraries and utilities for pre/post-processing SCHISM models',
    author='Zhengui Wang',
    author_email='wangzg@vims.edu',
    classifiers=[],
    install_requires=[
        'pathlib',
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
