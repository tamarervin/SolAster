"""
Tamar Ervin
Date: September 2, 2021

Looking at level 0 data.

"""

import sys

sys.path.append('/Volumes/Portable/NEID_Solar_analysis')

import os
import glob
import datetime
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy import constants
from sunpy.net import attrs as a

from astropy.io import fits

from scipy.integrate import simps

import NEIDcode

import tamar.tools.ccf_funcs as ccffuncs
from tamar.tools.settings import CsvDir, Config

# get list of fits files
all_files = glob.glob(os.path.join(CsvDir.LEVEL0, '*.fits'))
dates = [d.split('/')[-1] for d in all_files]
dates = [d.split('_')[-1] for d in dates]
dates = [d.split('T')[0] for d in dates]

# open file
f = all_files[0]
fd = fits.open(f)
flux = fd['SCIFLUX'].data
wave = fd['SCIWAVE'].data
head = fd[0].header
ccfhead = fd['PYR'].header
ccfwt = np.zeros(122)
main_header = fits.getheader(f)
date_jd = fits.getheader(f)['OBSJD']

sun_head = fd[20]
data = sun_head.data
dates = [d[0] for d in data]
voltage = [d[1] for d in data]
flux = [d[2] for d in data]