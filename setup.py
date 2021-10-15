from distutils.core import setup
from setuptools import find_packages

setup(
    name='sdo_hmi_rvs',
    version='0.0.4',  # Ideally should be same as your GitHub release tag varsion
    description="Python package to independently derive 'sun-as-a-star' radial velocity variations.",
    author='Tamar Ervin',
    author_email='tamarervin@gmail.com',
    url='https://github.com/tamarervin/sdo_hmi_rvs',
    download_url='https://github.com/tamarervin/sdo_hmi_rvs/archive/refs/tags/v0.0.4.tar.gz',
    keywords=['RVs', 'Solar'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
)


