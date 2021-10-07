# Data Used

This project requires the use of data from the Helioseismic and Magnetic 
Imager (HMI) aboard NASA's Solar Dynamics Observatory (SDO). 

HMI has been continuously imaging the Sun since 2012 in the 6173 A line. HMI
takes images every 45 seconds and provides information regarding the velocity, continuum, and magnetic 
field of the Solar photosphere. 

Additional information regarding the instrument can be found on the [HMI webpage](http://hmi.stanford.edu/).  

The three data products used for these calculations are:  

* **Dopplergrams**: measurements of solar surface velocity  
  
* **Filtergrams**: continuum intensity measurements of the Solar photosphere  
  
* **Magnetograms**: measurements of the line-of-sight Solar magnetic field  


## Acquiring the Data

To get the HMI data products, we utilize Sunpy's querying techniques. 
This allows us to search by date, physical observable, and instrument.  

An example of querying Sunpy for data can be found [here](../examples/docs_sunpy_example/).  
