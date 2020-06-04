import netCDF4 as nc
import pandas as pd
import numpy as np
import reverse_geocoder as rg
from multiprocessing import freeze_support, Pool

def main():
    freeze_support()
    # ____________________________Define inputs_______________________________________
    input_folder = "../input/"
    output_folder = "../output/"

    input_co2_mole = "mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-0_gr1-GMNHSH_2015-2500.nc"
    input_co2_air = "CO2-em-AIR-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc"

    # Input for the gases
    input_nox = "NOx-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc"
    input_bc = "BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc"
    input_so2 = "SO2-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc"
    input_oc = "OC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc"
    input_co = "CO-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc"
    input_nh3 = "NH3-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc"
    input_nmvoc = "NMVOC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc"
    input_blip = "Robin_sectors_V3.csv"
    convert_country_code_file = "convertCountryCodes.csv"
    files_to_blip = [
        input_co, input_nh3, input_nmvoc,
        #input_oc, input_so2, input_nox, input_bc
    ]
    key_variables = [
        "CO_em_anthro", "NH3_em_anthro", "NMVIC_em_anthro",
        # "OC_em_anthro", "SO2_em_anthro", "NOx_em_anthro", "BC_em_anthro"
    ]
    scenario_string = "_blip.nc"

    assert len(files_to_blip) == len(key_variables) # check input

    # The set of sectors in our blip need to be converted into our sectors in the netCDF
    # case. This uses:
    # 0: Agriculture; 1: Energy; 2: Industrial; 3: Transportation; 4: Residential,
    # Commercial, Other; 5: Solvents production and application; 6: Waste;
    # 7: International Shipping

    sector_dict = {
        "surface-transport": 3, "residential": 4, "public/commercial": -4,
        "industry": 2, "international-shipping": 7, "international-aviation": -1,
        "domestic-aviation": -2, "power": 1
    }
    # We will manage sectors 6 elsewhere, no change to sector 0 (agri).
    sectors_to_use = [1, 2, 3, 4, 5, 7]
    # __________________________________________________________________________________

    # Collect and clean the data
    # We need to make the blip factors consistent with the netCDF data structure. This will
    # require example data, although the results should not depend which example is chosen.

    nox_0 = nc.Dataset(input_folder + input_nox, "r", format="NETCDF4")
    blip_factors = pd.read_csv(input_folder + input_blip)
    convert_countries = pd.read_csv(
        input_folder + convert_country_code_file, keep_default_na=False, na_values=['_']
    )

    blip_factors = blip_factors[~blip_factors["1"].isna()]

    # Perform the sector weighting
    blip_factors_multi = blip_factors.copy()
    blip_factors_multi.drop(
        ["Country", "Base(MtCO2/day)", "Unnamed: 0"], axis=1, inplace=True
    )
    blip_factors_multi["Sector"] = [sector_dict[sect] for sect in blip_factors_multi["Sector"]]
    blip_factors_multi.set_index(
        blip_factors_multi.columns[:2].to_list(), drop=True, inplace=True
    )

    # We want to average the two sets of sector 4 together in the right ratio
    all_countries = blip_factors_multi.index.get_level_values("ISO_A3").unique()
    for country in all_countries:
        if (country, 4) in blip_factors_multi.index and (country, -4) in blip_factors_multi.index:
            blip_factors_multi.loc[country, 4] = (
                blip_factors_multi.loc[country, 4].values *
                blip_factors_multi["Base%"][country, 4] +
                blip_factors_multi.loc[country, -4].values *
                blip_factors_multi["Base%"][country, -4]
                ) / (
                    blip_factors_multi["Base%"][country, 4] +
                    blip_factors_multi["Base%"][country, -4]
                )
            blip_factors_multi["Base%"][country, 4] = blip_factors_multi["Base%"][country, 4] + \
                blip_factors_multi["Base%"][country, -4]
            blip_factors_multi.drop((country, -4), inplace=True)
        elif (country, -4) in blip_factors_multi.index:
            blip_factors_multi.loc[country, 4] = blip_factors_multi.loc[country, -4]
        elif (country, 4) in blip_factors_multi.index:
            continue
        else:
            print("no data for {}".format(country))

    # Test that this produces the right answers
    example_factor = blip_factors[
        (blip_factors["ISO_A3"] == "GBR") &
        (blip_factors["Sector"].isin(["residential", "public/commercial"]))
    ][["Base%", "100"]]
    assert np.isclose(blip_factors_multi.loc["GBR", 4][100], sum(
        example_factor["Base%"] * example_factor["100"]) / sum(example_factor["Base%"])
    )

    # Also assume solvent tracks industry
    for country in all_countries:
        blip_factors_multi.loc[(country, 5), :] = blip_factors_multi.loc[(country, 4), :]

    # Derive country and date relation
    # We need to assign each lat/long a country. This is slightly complicated by the
    # country index being 2 letters in the inverse geocoder but 3 letters in our data.

    lat, lon = nox_0.variables["lat"][:], nox_0.variables["lon"][:]

    convert_countries_dict = {
        convert_countries["A2 (ISO)"][i]: convert_countries["A3 (UN)"][i]
        for i in convert_countries.index
    }
    coords = []
    lon_length = len(lon)
    for latperm in lat:
        coords = coords + list(zip([latperm] * lon_length, lon))

    # Find the countries of each coordinate:
    results = rg.search(coords)
    lat_countries_dict = {
        coords[i]: convert_countries_dict[results[i]["cc"]] for i in range(len(coords))
                          if results[i]["cc"] in convert_countries_dict.keys()
    }
    # The process will be faster if we map the other way and use the index rather than
    #  the coordinates
    country_coord_dict = {}
    for k, v in lat_countries_dict.items():
        country_coord_dict[v] = country_coord_dict.get(v, [])
        country_coord_dict[v].append(
            (np.where(lat.data == k[0])[0][0], np.where(lon.data == k[1])[0][0])
        )

    # Now we must relate the dates. blip_factors uses days from 2020-01-01, and has
    # values for every day. The netCDFs use days since 2015-01-01, which is 5 * 365 + 1
    # days later and monthly.

    date_dif = 5 * 365 + 1
    netCDF_times = nox_0.variables["time"][:]
    netCDF_tseries = pd.Series(netCDF_times)
    bliptimes = blip_factors_multi.columns[blip_factors_multi.columns != "Base%"]
    bliptimes = pd.Series(pd.to_numeric(bliptimes))
    time_dict = {}
    remaining_times = bliptimes.copy()
    mappable_times = netCDF_tseries[
        (netCDF_tseries > date_dif) & (netCDF_tseries < date_dif + max(bliptimes))
    ]
    for t in mappable_times.index[:-1]:
        closeTimes = [bliptime for bliptime in remaining_times if (
            0.5 * (mappable_times[t + 1] + mappable_times[t]) - date_dif > bliptime
        )]
        time_dict[mappable_times[t]] = closeTimes
        remaining_times = remaining_times[~remaining_times.isin(closeTimes)]
    time_dict[mappable_times.iloc[-1]] = list(remaining_times)
    blip_factors_av = pd.DataFrame(
        index=blip_factors_multi.index, columns=time_dict.keys()
    )
    for key, val in time_dict.items():
        blip_factors_av[key] = blip_factors_multi[
            list(str(v) for v in val)
        ].mean(axis=1)

    # Perform the emissions blip
    # We now have a mapping between times and locations and the emissions we want.

    nox_0.close()
    all_valid_countries = [c for c in all_countries if c in country_coord_dict.keys()]

    for fileind in range(len(files_to_blip)):
        file = files_to_blip[fileind]
        print("Working on file {}".format(file))
        data = copy_netcdf_file(file, input_folder, output_folder, scenario_string)
        output = data.variables[key_variables[fileind]]
        for country in all_valid_countries:
            print(country)
            for time in blip_factors_av.columns:
                timeind = np.where(data.variables["time"][:] == time)[0]
                for sector in sectors_to_use:
                    try:
                        mult_fact = blip_factors_av[time].loc[country, sector] + 1
                        if mult_fact != 1.0:  # This saves operations
                            for lati, longi in country_coord_dict[country]:
                                output[timeind, sector, lati, longi] *= mult_fact
                    except KeyError as e:
                        print("Key error for country {}".format(e))
                        break
        data.variables[key_variables[fileind]] = output
        data.close()


def copy_netcdf_file(filename, input_folder, output_folder, scenario_string):
    src = nc.Dataset(input_folder + filename)
    trg = nc.Dataset(output_folder + filename + scenario_string, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions, complevel=9)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]

    # Return the data
    src.close()
    return trg


if __name__ == "__main__":
    main()