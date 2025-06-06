# Condensed phase properties


This repo contains the functions and exemplary input data for water that were used to compute thermodynamic properties. In the data folder, all published data is available (thermodynamic data and diffusion constants).


## Features

All functions used to compute thermodynamic properties can be found in the `compute_properties.py` script. In the data folder, there is some examplary data from a NPT water simulation using TIP3 water. In the data/all_property_data folder, all raw data (thermodynamic properties and diffusion constants) can be found.

## Examplary scripts

The python script `compute_properties.py` uses the water report in the data folder and the property functions to compute and print thermodynamic properties. It also computes the mean and standard deviation from NBoot = 100 bootstrapping samples from the water data using the property bootstrapping functions in the `cp_props.py` script. The python script `get_available_properties.py` collects and prints the raw data (stored in the data folder) for a given species. To execute those scripts (once installed this repository), run

- `python compute_properties.py`

or

- `python get_available_properties.py water`

water can be replaced by methanol, acetone, nma, benzene or hexane.

## Other used packages

For the calculation of MSD and diffusion constants with respect to the center of mass, the NewAnalysis package was used (https://github.com/cbc-univie/mdy-newanalysis-package). Radial distribution functions were computed using MDAnalysis (https://github.com/MDAnalysis/mdanalysis).

