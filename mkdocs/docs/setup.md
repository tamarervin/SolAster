# Setup

In order to run examples and scripts from SolAster, the Inputs class must be updated. The inputs
class includes information about the instrument for comparison, date frame, and calculation cadence. This class
is included at the start of each example script and includes all required user inputs.  

* csv_name: name of the CSV file to store calculations (str)  
* inst: name of instrument to use for calcuations of RV model (str: either NEID' or 'HARPS-N')  
* cadence: querying cadence in seconds (int)  
* start_date: start date for calculations (datetime)  
* end_date: end date for calculations (datetime)  
* diagnostic_plots: whether you would like to plot diagnostic plots (bool)  
* save_fig: path to save figures if diagnostic_plots is True (str)  

Additionally, users can update paths to store CSV files in settings. The current paths are setup to save directly to the downloaded repository but can be changed for different systems if needed.  

The input class looks like this:  
```
class Inputs:
    """
    Class to hold user specified inputs to run examples.
    See README or documentation site for additional information.
    """

    # name of csv file to store calculations
    csv_name = 'example.csv'

    # name of instrument to use for calculation of RV model
    # choose either 'NEID' or 'HARPS-N'
    inst = 'NEID'

    # querying cadence in seconds
    cadence = 24 * 60 * 60

    # start date for calculationsx
    start_date = datetime.datetime(2021, 2, 10, 0, 0, 0)

    # end date for calculations
    end_date = datetime.datetime(2021, 2, 14, 0, 0, 0)

    # True if outputting diagnostic plots
    diagnostic_plots = True
    # path to save diagnostic figure or none
    save_fig = None
```  

# Outputted results

Our package produces a CSV or pickle file which includes calculation results of velocity components, model RV variations, and 
various solar observables. These results are stored as follows.

* date_obs: calculation time in UT (str)  
* date_jd: calculation time in JD (float)  
* rv_model: model RV variation [m/s]  
* v_quiet: quiet-Sun velocity [m/s]  
* v_disc: velocity of full solar disk [m/s]  
* v_phot: photometric velocity component [m/s]  
* v_conv: convective velocity component [m/s]  
* f_bright: filling factor due to bright regions [%]  
* f_spot: filling factor due to spots [%]  
* f: filling factor [%]  
* Bobs: unsigned magnetic flux [G]  
* vphot_bright: photometric velocity component due to bright regions [m/s]  
* vphot_spot: photometric velocity component due to spots [m/s]  
* f_small: filling factor due to small regions [%]  
* f_large: filling factor due to large regions [%]  
* f_network: filling factor due to network regions [%]  
* f_plage: filling factor due to plage regions [%]  
* quiet_flux: magnetic flux due to quiet-Sun regions [G]  
* ar_flux: magnetic flux due to active Sun regions [G]  
* conv_flux: magnetic flux due to large active regions [G]  
* pol_flux: polarized magnetic flux [G]  
* pol_conv_flux: polarized magnetic flux due to large active regions [G]  
* vconv_quiet: convective velocity component due to quiet-Sun regions [m/s]  
* vconv_large: convective velocity component due to large active regions [m/s]  
* vconv_quiet: convective velocity component due to small active regions [m/s]  

