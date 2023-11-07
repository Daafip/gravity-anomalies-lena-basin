# Mass anomalies obtained by gravimetry in the Lena Basin
_Conversion of stoke coefficients to mass anaomalies using GRACE data in the LENA basin_

## Overview of code: 
- `Project_grace_filtered.ipynb` contains the main bulk of code and should be understandable 
- `Data overview.ipynb` contains some helper functions to speed up data collection etc.
- Other notebooks can be found in `\old code` these are poorly commented and not worked out in detail 

### `.py` files
- Multi core processing speeds up the computation times, jupyter lab is less usefull for this, thus these are run in:
    > - `multicore_filtered.py` and `multicore_raw.py` which run in around 1/6th of the time from my testing (depening on your machine performance offcouse)  
- `calc_geoid_change_alt.py` is used to run the normalized Legendre function, this was provided by P. Ditmar (TU Delft)
- `geospatial_functions.py` is a combination of functions to plot background maps easily using contexitly

<!-- - `Forcing.txt` contains the real observational data
- `Callibration model.ipynb` contains the code running the function from top to bottom
- `HBVMod.py` contains the actual model
- `MC2_NSE.txt`,`MC2_NSE_log.txt`,`MC2_NSE_sqrt.txt` are three observational runs - MCS of 5000 runs takes a minute thus we store the result to improve experience
- `Weigfun.py` contains a weighting function to add lag to the model
- `Forward model.ipynb` old version of forward model - timestep based - running but not OOP
- `Forward model BMI.ipynb` Newer version: Object Oriented Programming, where the model state is a class instance. running model.update() advances the model timestep according to the Basic Model Interface (BMI). -->