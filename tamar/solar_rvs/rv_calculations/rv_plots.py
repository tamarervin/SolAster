"""
Tamar Ervin
Date: June 24, 2021

plotting of RVs and other relevant observables and components

"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from tamar.tools.settings import CsvDir

# csv file with rv information
csv_file = os.path.join(CsvDir.NEID_CALC, 'rv_calcs.csv')

# create pandas dataframe
component_df = pd.read_csv(csv_file)

# get column names
components = component_df.columns.values
print("Column Names:", components)

# to get a component name
# name = component_df.name.values

# get components
# date_str, v_quiet, v_phot, v_disc, v_conv, filling_factor, unsigned_obs_flux, unsigned_rad_flux, rsun =\
#     [None] * len(components)

# for i, component in enumerate(components):
#     exec("%s = %n" % (component, component_df.iloc[i].values))

date_obs = component_df.date_obs.values
v_quiet = component_df.v_quiet.values
v_phot = component_df.v_phot.values
v_disc = component_df.v_disc.values
v_conv = component_df.v_conv.values
filling_factor = component_df.filling_factor.values
unsigned_obs_flux = component_df.unsigned_obs_flux.values
unsigned_rad_flux = component_df.unsigned_rad_flux.values
rv_sun = component_df.rv_sun.values

# make some plots!!

# make a dates range

# RVs
date_range = range(0, len(date_obs))
title = "Solar RV Variation in First Part of Haywood"
color = 'purple'
x = date_range
y = rv_sun
xlabel = "Dates since: " + date_obs[0]
ylabel = "Solar RV Variation"
plt.figure()
plt.title(title)
plt.scatter(x, y, color=color)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.show()

# v photometric
date_range = range(0, len(date_obs))
title = "Filling Factor Haywood"
color = 'purple'
x = date_range
y = filling_factor * 100
xlabel = "Dates since: " + date_obs[0]
ylabel = "Solar RV Variation"
plt.figure()
plt.title(title)
plt.scatter(x, y, color=color)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.show()

# v convective
date_range = range(0, len(date_obs))
title = "Convective Velocity in First Part of Haywood"
color = 'purple'
x = date_range
y = v_conv
xlabel = "Dates since: " + date_obs[0]
ylabel = "Solar RV Variation"
plt.figure()
plt.title(title)
plt.scatter(x, y, color=color)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.show()


