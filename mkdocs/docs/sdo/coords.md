# Coordinate Transformations

SDO/HMI Data Products are returned in Helioprojective Cartesian pixel coordinates. In 
order to calculate the spacecraft and solar rotation velocity corrections, we need to
convert this to the Heliographic Carrington Coordinate System (sun centered coordinate system). 
In addition to being sun-centered, Heliographic Carrington coordinates also rotate with
respect to the Carrington rotation period. Transformations are based off Thompson (2006). 

Coordinate functions are [here](https://github.com/shalverson/NEID_Solar_analysis/blob/master/tamar/tools/coord_funcs.py).

## Pixel Array Computation

The first step in the coordinate transformation is computation of the pixel-wise
coordinate arrays with respect to solar center. We do this by taking into account the
reference pixels from the header file. 

## Calculation of Mu and Radius Array

We then calculate the array of $\mu$ and radius values for that image, where $\mu = cos(\theta)$. 
First, we determine pixel size in physical coordinates (solar radii) and then calculate
the center-to-limb angles. We get the radius by determining the distance in solar radii 
to each pixel. 

## Coordinate Transform

Finally, we are able to apply the coordinate transformation based on a calculated rotation 
matrix. We get an array of coordinates in the westward, northward, and radial direction. 

