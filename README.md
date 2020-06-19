# Modifying emission and concentration data in response to the COVID-19 pandemic

This repository takes a set of data from input4mips 
(https://esgf-node.llnl.gov/search/input4mips/) and modification factors from 
https://github.com/Priestley-Centre/COVID19_emissions to output netCDF files 
suitable for use after the pandemic. The files modified are as follows:

* lat/longitude gridded aerosol/precursor emissions by grid location and sector, 
    * will be modified by country activity levels
    * Requires the following files in "input/aerosols":
    "{}-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc"
    for {} = BC, CO, NH3, NOx, NMVOC, OC, SO2
* hemispheric and global well-mixed GHG emissions
    * will be modified only by global concentration levels
    * requires the files in "input":
    "mole-fraction-of-{}-in-ar_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-1_gr1-GMNHSH_201501-250012.nc"
    for {} = carbon-dioxide, methane, nitrous-oxide.
 * lat/long/altitude gridded aviation emissions for NOx. 
 
 In all cases, the modification factors are taken from the Priestly github data, which 
 is assumed to be in a sister folder to this folder, entitled "COVID19_emissions_data". 