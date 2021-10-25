import logging
import sys
from io import open
from os import path

try:
    from setuptools import setup, find_packages
except ImportError:
    logging.exception('Please install or upgrade setuptools or pip.')
    sys.exit(1)

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='sdo_hmi_rvs',
    version='0.0.7',
    description="Python package to independently derive 'sun-as-a-star' radial velocity variations.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/tamarervin/sdo_hmi_rvs',
    download_url='https://github.com/tamarervin/sdo_hmi_rvs/archive/refs/tags/v0.0.7.tar.gz',
    keywords=['RVs', 'Solar'],
    author='Tamar Ervin',
    author_email='tamarervin@gmail.com',
    classifiers=[
                "Programming Language :: Python :: 3",
                "License :: OSI Approved :: Apache Software License",
                "Operating System :: OS Independent",
            ],
    packages=['sdo_hmi_rvs',
              'sdo_hmi_rvs.tools',
              'sdo_hmi_rvs.source',
              'sdo_hmi_rvs.examples',
              # 'sdo_hmi_rvs.tools.calculation_funcs',
              # 'sdo_hmi_rvs.tools.coord_funcs',
              # 'sdo_hmi_rvs.tools.lbc_funcs',
              # 'sdo_hmi_rvs.tools.plotting_funcs',
              # 'sdo_hmi_rvs.tools.rvs',
              # 'sdo_hmi_rvs.tools.settings',
              # 'sdo_hmi_rvs.tools.utilities',
              # 'sdo_hmi_rvs.tools.settings_template',
              ],

    package_data={
        # If any package contains *.txt, *.rst or *.fits files, include them:
        '': ['*.txt', '*.yaml', '*.yml', '*.fits', '*.pdf', '*.dat', '*.csv']
    },
    include_package_data=True,

    python_requires='>=3.7',
    install_requires=[
        'astropy',
        'matplotlib',
        'numpy',
        'sunpy',
        'pandas',
        'scipy',
        'scikit-image',
        'scikit-learn'
    ]
)




