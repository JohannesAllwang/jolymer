from setuptools import find_packages, setup

setup(
    name='jolymer',
    packages=find_packages(),
    # package_dir = {'':'jolib'},
    version='0.0.2',
    description='Library to assist the handling of experimental data',
    author='Johannes Allwang',
    include_package_data=True,
    license='MIT',
    package_dir={ "": "src" },
    package_data={
        # Include the script in the package
        "": ["scripts/smooth-saxs-curve.sh"],
    },
)
