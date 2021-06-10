from setuptools import find_packages, setup

setup(
    name='jolymer',
    packages=find_packages(),
    # package_dir = {'':'jolib'},
    version='0.0.1',
    description='My first Python library',
    author='Johannes Allwang',
    include_package_data=True,
    license='MIT',
    package_dir={ "": "src" }
)
