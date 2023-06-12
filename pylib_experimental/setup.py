from setuptools import setup

with open('requirements.txt', 'r') as f:
  requirements = f.read().splitlines()

  setup(
    name='pylib_experimental',
    version='0.0.1',
    description='An experimental package of pylibs',
    license='MIT',
    packages=[
      'pylib_essentials',
    ],
    package_data={},
    install_requires=requirements,
  )

