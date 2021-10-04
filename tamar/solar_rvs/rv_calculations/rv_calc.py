"""
Tamar Ervin
Date: June 24, 2021

Calculation of Solar RV variations based on components calculated in
tamar/solar_rvs/rv_components.py

Goal: Comparison with RV variations calculated in Haywood et al. 2016

This takes a matter of seconds!!

"""

import os
import pandas as pd
from tamar.tools.settings import CsvDir

# csv file with rv components
csv_file = os.path.join(CsvDir.MILBOURNE, 'long_bad_dates.csv')

# scaling coefficient values -- from Haywood 2016
A = 2.45
B = 1.85
RV0 = 99.80

# create pandas dataframe
component_df = pd.read_csv(csv_file)

# get column names
components = component_df.columns.values

# get velocities lists
v_phot = component_df.v_phot.values
v_conv = component_df.v_conv.values

# RV calculation
RV = A*v_phot + B*v_conv + RV0

# add RV column
component_df["rv_model"] = RV
component_df.to_csv(csv_file, index=False)







