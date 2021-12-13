# SDO HMI RVs

Pipeline to independently derive 'Sun-as-a-star' radial velocity variations.

# Documentation

**Documentation Site:**  https://tamarervin.github.io/SolAster/

# Build conda environment

* update dependencies in conda_env.yml [file](conda_env.yml)   
* run the following from the folder containing the .yml file
    * ``conda env create -f conda_env.yml``  
* to add new dependencies, update conda_env.yml [file](conda_env.yml)  
* run the following from the folder containing the .yml file  
    * ``conda env update -f conda_env.yml``
    
# Installation

* package installation using pip  
* install pip  
* install package   
``pip install SolAster``  
  
# References  

* Ervin et al. (2021) - In Preparation  
* [Milbourne et al. (2019)](https://doi.org/10.3847/1538-4357/ab064a)  
* [Haywood et al. (2016)](https://doi.org/10.1093/mnras/stw187)  
* Based on a technique developed by [Meunier, Lagrange & Desort (2010)](https://doi.org/10.1051/0004-6361/200913551) 
  for SoHO/MDI images.  


# Examples
Examples are hosted [here](https://github.com/tamarervin/SolAster/tree/main/sdo_hmi_rvs/examples):  

1. [Sunpy Example](https://github.com/tamarervin/SolAster/blob/main/sdo_hmi_rvs/examples/sunpy_example.ipynb): 
outlines how to use basic Sunpy functions and usages for this package  
   
2. [Component Calculations](https://github.com/tamarervin/SolAster/blob/main/sdo_hmi_rvs/examples/component_calculations.ipynb): 
outlines the corrections and component calculation pipeline  
   * creates CSV with calculations of magnetic observables and velocity components  
  
3. [RV Calculations](https://github.com/tamarervin/SolAster/blob/main/sdo_hmi_rvs/examples/rv_calculation.ipynb):
outlines calculation of full model RV from velocity components  
   * requires input CSV with velocity components from [example 2](https://github.com/tamarervin/SolAster/blob/main/SolAster/examples/component_calculations.ipynb)  
   * an example CSV file with calculations is stored [here](https://github.com/tamarervin/SolAster/blob/main/SolAster/products/csv_files/calcs/example_calcs.csv)
    
4. [Full Pipeline](https://github.com/tamarervin/SolAster/blob/main/SolAster/examples/full_pipeline.ipynb):
full end-to-end pipeline to calculate 'sun-as-a-star' RVs and magnetic observables 
   
# Issues or Suggestions

* for any issues or bug fixes, please fill out an issue report on the [GitHub page](https://github.com/tamarervin/SolAster/issues)  

# Contact

**Tamar Ervin**: <tamarervin@gmail.com>

**Sam Halverson**: <samuel.halverson@jpl.nasa.gov>
