from setuptools import setup

with open('requirements.txt', 'r') as f:
    requirements = f.read().splitlines()

setup(
    name='pylib_essentials',
    version='0.0.1',
    description='The essentials of the experimental package of pylibs, which require minimum dependencies',
    license='MIT',
    packages=['pylib_essentials'],
    package_data={},
    install_requires=requirements,
)

