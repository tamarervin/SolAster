# Calculation of Model RVs

After we calculate the velocity components, we can move to the calculation
of SDO derived RV variations. To create the model RVs, we use linear regression 
to determine scaling coefficients for the two velocity components 
($v_{phot}$, $v_{conv}$) along with an offset value.  

Because we do not know the relative weights of these two velocity components, we
rely on regression to determine the best parameters to fit our ground based data. 
