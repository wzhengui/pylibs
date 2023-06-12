from setuptools import setup

with open('requirements.txt', 'r') as f:
  requirements = f.read().splitlines()

  setup(
    name='pylib_essentials',
    version='0.0.1',
    description="An essential package only including pylib's core functions",
    license='MIT',
    packages=[
    ],
    package_data={},
    install_requires=requirements,
  )

