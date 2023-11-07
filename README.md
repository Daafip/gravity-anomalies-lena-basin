# Mass anomalies obtained by gravimetry in the Lena Basin
_Conversion of stoke coefficients to mass anaomalies using GRACE data in the LENA basin_

# Explanation of repo
## Overview of code 
- `Project_grace_filtered.ipynb` contains the main bulk of code and should be understandable 
- `Data overview.ipynb` contains some helper functions to speed up data collection etc.
- Other notebooks can be found in `\old code` these are poorly commented and not worked out in detail 

### The `.py` files
- Multi core processing speeds up the computation times, jupyter lab is less usefull for this, thus these are run in:
    - `multicore_filtered.py` and `multicore_raw.py` which run in around 1/6th of the time from my testing (depening on your machine performance offcouse)  
- `calc_geoid_change_alt.py` is used to run the normalized Legendre function, this was provided by P. Ditmar (TU Delft)
- `geospatial_functions.py` is a combination of functions to plot background maps easily using contexitly

## Data
Contains an assortment of data, mainly the processed geospatial files. Download stokes coefficients yourself using `Data overview.ipynb` from the owners.

## Files
Contains report and presentation.
