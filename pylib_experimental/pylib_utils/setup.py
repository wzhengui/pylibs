from setuptools import setup

with open('requirements.txt', 'r') as f:
    requirements = f.read().splitlines()

setup(
    name='pylib_utils',
    version='0.0.1',
    description="A subset of functions from pylib's mylib.py",

    license='MIT',
    packages=['pylib_utils'],
    package_data={},
    install_requires=requirements,
)

