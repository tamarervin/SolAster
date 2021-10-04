#!/usr/bin/env python
"""
Tamar Ervin
Date: June 30, 2021

Optimization of scaling factors for RV
calculation.
"""

import os
import sys

sys.path.append('/Users/tervin/NEID_Solar_analysis')

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from tamar.tools.settings import CsvDir


# csv file with rv components
# csv_file = os.path.join(CsvDir.NEID_CALC, 'rvs_from_fits.csv')
csv_file = os.path.join(CsvDir.NEID_CALC, 'paper_rvs.csv')

# create pandas dataframe
component_df = pd.read_csv(csv_file)

# get dates list
date_jd = component_df.date_jd.values
inds = np.argsort(date_jd)

rv_sun = component_df.rv_sun.values[inds]
rv_error = component_df.rv_error.values[inds]

non_nan = np.logical_not(np.isnan(rv_sun))

rv_med = np.nanmedian(np.abs(rv_sun))

good_sun = np.logical_and(np.abs(rv_sun) > rv_med - 2, np.abs(rv_sun) < rv_med + 2)
good_error = np.logical_and(np.abs(rv_error) < .4, np.abs(rv_error) < .225)
good = np.logical_and(good_sun, good_error)

good_rvs = np.logical_and(good, non_nan)

# get velocities lists
v_phot = component_df.v_phot.values[inds][good_rvs]
v_conv = component_df.v_conv.values[inds][good_rvs]
rv_model = component_df.rv_model_weather_coeff.values[inds][good_rvs]
rv_sun = rv_sun[good_rvs]
rv_error = rv_error[good_rvs]
v_quiet = component_df.v_quiet.values[inds][good_rvs]
v_disc = component_df.v_disc.values[inds][good_rvs]
vconv_large = component_df.vconv_large.values[inds][good_rvs]
vconv_small = component_df.vconv_small.values[inds][good_rvs]
vphot_bright = component_df.vphot_bright.values[inds][good_rvs]
vphot_spot = component_df.vphot_spot.values[inds][good_rvs]

# get magnetic observables
f = component_df.f.values[inds][good_rvs]
Bobs = component_df.Bobs.values[inds][good_rvs]
f_bright = component_df.f_bright.values[inds][good_rvs]
f_spot = component_df.f_spot.values[inds][good_rvs]
f_plage = component_df.f_plage.values[inds][good_rvs]
f_network = component_df.f_plage.values[inds][good_rvs]

# make arrays
X = np.zeros(shape=(len(rv_sun), 2))
X[:, 0] = v_phot
X[:, 1] = v_conv
y = rv_sun

# apply regression
reg = LinearRegression().fit(X, y)
print('R^2 for prediction: ' + str(reg.score(X, y)))

# get scaling factors
A = reg.coef_[0]
B = reg.coef_[1]
RV0 = reg.intercept_

# Calculate matrix of predicted class probabilities.
predProbs = reg.predict_proba(X)
X_design = np.hstack([np.ones((X.shape[0], 1)), X])
V = np.diagflat(np.product(predProbs, axis=1))
# get covariance
covLogit = np.linalg.inv(np.dot(np.dot(X_design.T, V), X_design))
print("Covariance matrix: ", covLogit)

# Standard errors
print("Standard errors: ", np.sqrt(np.diag(covLogit)))

print(A, B, RV0)

# # calculate the Model RV for the milbourne data using these coefficients
RV = A*v_phot + B*v_conv + RV0
# component_df["rv_model"] = RV
# component_df.to_csv(csv_file, index=False)

# RV = v_phot + v_conv
# component_df["rv_sum"] = RV
# component_df.to_csv(csv_file, index=False)

# mixing regression methods
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import VotingRegressor

reg1 = GradientBoostingRegressor(random_state=1)
reg2 = RandomForestRegressor(random_state=1)
reg3 = LinearRegression()

Xfit = X[:20]
yfit = y[:20]

reg1.fit(Xfit, yfit)
reg2.fit(Xfit, yfit)
reg3.fit(Xfit, yfit)

ereg = VotingRegressor([('gb', reg1), ('rf', reg2), ('lr', reg3)])
ereg.fit(Xfit, yfit)

xt = X
pred1 = reg1.predict(xt)
pred2 = reg2.predict(xt)
pred3 = reg3.predict(xt)
pred4 = ereg.predict(xt)

import matplotlib.pyplot as plt
plt.figure()
plt.plot(pred1, 'gd', label='GradientBoostingRegressor')
plt.plot(pred2, 'b^', label='RandomForestRegressor')
plt.plot(pred3, 'ys', label='LinearRegression')
plt.plot(pred4, 'r*', ms=10, label='VotingRegressor')
plt.plot(rv_sun, 'ko', alpha=0.5, label='NEID RV')

plt.tick_params(axis='x', which='both', bottom=False, top=False,
                labelbottom=False)
plt.ylabel('Predicted RV')
plt.xlabel('Samples')
plt.legend(loc="best")
plt.title('Regressor predictions and their average compared to NEID RVs')

plt.show()
#
RV_reg = pred1
component_df["grad_reg"] = RV_reg
component_df.to_csv(csv_file, index=False)

### using support vector machine regression
# from sklearn import svm
# regr = svm.LinearSVR()
# fit = regr.fit(X, y)
#
# # get scaling factors
# A = fit.coef_[0]
# B = fit.coef_[1]
# RV0 = fit.intercept_
#
# RV_svm = A*v_phot + B*v_conv + RV0
#
# import tamar.tools.plotting_funcs as plot
# import matplotlib.pyplot as plt
# plot.vert_comp_timeseries(range(0, 57), [RV, RV_svm], title='Compare Methods', xlabel='Observations', ylabel_list=['Linear', 'LinearSVM'])
# plt.show()