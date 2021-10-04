"""
Tamar Ervin

"""

import sys

sys.path.append('/Users/tervin/NEID_Solar_analysis')

import os
import glob
import datetime
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy import constants
from sunpy.net import attrs as a

from astropy.io import fits

import NEIDcode

import tamar.tools.ccf_funcs as ccffuncs
from tamar.tools.settings import CsvDir, Config

order_rvs = os.path.join(CsvDir.NEID_CALC, 'order_rvs.csv')
df = pd.read_csv(order_rvs)

# read in order RVs
rv_orders = []
order_numbers = np.arange(57, 170)
for i in order_numbers:
    i = str(i)
    i = i.zfill(3)
    rv_orders.append(df[i])

# get median value for each order
med, use_med = [], []
for j in rv_orders:
    med.append(np.nanmedian(j))
    if np.nanmedian(j) != 0:
        use_med.append(np.nanmedian(j))

med_val = np.nanmedian(use_med)

# calculate offsets
delta_rv = []
for j in med:
    if j != 0:
        delta_rv.append(j - med_val)
    else:
        delta_rv.append(0)

delta_csv = os.path.join(CsvDir.NEID_DATA, 'delta_rv_orders.csv')
data = {
    'order': np.arange(57, 170),
    'offset': delta_rv
}
df = pd.DataFrame(data)
df.to_csv(delta_csv)