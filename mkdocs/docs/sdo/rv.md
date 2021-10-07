# Calculation of Model RVs

The model 'sun-as-a-star' RV variation is a linear combination of the two velocity 
components $\rm (v_{phot}, v_{conv})$ along with an offset factor of $\rm RV_0$: 

$RV = A \ v_{phot} + B \ v_{conv} + RV_0$

Our ability to take high resolution measurements of the Solar disk is what 
allows us to decompose the overall RV variation into these two components. 

Because we do not know the relative weights of these two velocity components, we
rely on regression to determine the best parameters to fit our ground based data.

We calculate scaling factors A and B, and then offset $\rm RV_0$ using linear regression. The need for these scaling
coefficients arises due to the fact that SDO/HMI images the Sun in one wavelength (6173 A)
while ground based measurements average over thousands of lines. The offset accounts for the zero 
point of the instrument.  

An example of the linear regression model used can be found [here](../examples/docs_rv_calcs/).
