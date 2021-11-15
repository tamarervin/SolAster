# SDO HMI RVs

Pipeline to independently derive 'sun-as-a-star' radial velocity variations.

# Installation

* package installation using pip  
* install pip  
* install package   
``pip install sdo_hmi_rvs``

# Documentation

**Documentation Site:**  https://tamarervin.github.io/sdo_hmi_rvs/

# Build conda environment

* update dependencies in conda_env.yml [file](conda_env.yml)   
* run the following from the folder containing the .yml file
    * ``conda env create -f conda_env.yml``  
* to add new dependencies, update conda_env.yml [file](conda_env.yml)  
* run the following from the folder containing the .yml file  
    * ``conda env update -f conda_env.yml``  
  
# Examples
Examples are hosted [here](https://github.com/tamarervin/sdo_hmi_rvs/tree/main/sdo_hmi_rvs/examples):  

1. [Sunpy Example](https://github.com/tamarervin/sdo_hmi_rvs/blob/main/sdo_hmi_rvs/examples/sunpy_example.ipynb): 
outlines how to use basic Sunpy functions and usages for this package  
   
2. [Component Calculations](https://github.com/tamarervin/sdo_hmi_rvs/blob/main/sdo_hmi_rvs/examples/component_calculations.ipynb): 
outlines the corrections and component calculation pipeline  
   * creates CSV with calculations of magnetic observables and velocity components  
  
3. [RV Calculations](https://github.com/tamarervin/sdo_hmi_rvs/blob/main/sdo_hmi_rvs/examples/rv_calculation.ipynb):
outlines calculation of full model RV from velocity components  
   * requires input CSV with velocity components from [example 2](https://github.com/tamarervin/sdo_hmi_rvs/blob/main/sdo_hmi_rvs/examples/component_calculations.ipynb)  
   * an example CSV file with calculations is stored [here](https://github.com/tamarervin/sdo_hmi_rvs/blob/main/sdo_hmi_rvs/products/csv_files/calcs/example_calcs.csv)
    
4. [Full Pipeline](https://github.com/tamarervin/sdo_hmi_rvs/blob/main/sdo_hmi_rvs/examples/full_pipeline.ipynb):
full end-to-end pipeline to calculate 'sun-as-a-star' RVs and magnetic observables 
   
