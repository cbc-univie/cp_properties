import pandas as pd
import sys




######################## PRINT AVAILABLE DATA #################################################




def load_initial_properties(species: str):
    file_name_props = f"cp_props_{species}.csv"
    data_path_props = f"data/all_property_data/{species}/thermodynamic_properties/initial_npt_run/{file_name_props}"
    props = pd.read_csv(data_path_props, sep="&")
    file_name_timeseries = f"cp_props_time_series_{species}.csv"
    data_path_timeseries = f"data/all_property_data/{species}/thermodynamic_properties/initial_npt_run/{file_name_timeseries}"
    timeseries = pd.read_csv(data_path_timeseries, sep="&")
    return props, timeseries

def load_repetition_properties(species: str, method: str):
    file_name = f"{method}_props_per_run.csv"
    data_path = f"data/all_property_data/{species}/thermodynamic_properties/repetition_runs/{file_name}"
    return pd.read_csv(data_path, sep="&")
    
def load_msd_diffusion(species: str, method: str):
    file_name_diffusion = f"diffusion_{species}_{method}_NVT_NH.dat"
    data_path = f"data/all_property_data/{species}/diffusion/{file_name_diffusion}"
    diffusion = pd.read_csv(data_path, sep="&")
    file_name_msd = f"msd_{species}_{method}_NVT_NH.dat"
    data_path = f"data/all_property_data/{species}/diffusion/{file_name_msd}"
    msd = pd.read_csv(data_path, sep="&")
    return msd, diffusion




species = sys.argv[1]
if species not in ['water', 'methanol', 'acetone', 'nma', 'hexane', 'benzene']:
    print('Data available for water, methanol, acetone, nma, hexane, benzene. \n' \
    'Please use one of those!\n')
    exit()

props, timeseries = load_initial_properties(species)

print("\n")
print(f"Properties from initial NPT 1ns run for {species}")
print(props)
print("-"*50)
print("\n")

if species in ['water', 'acetone', 'hexane']:
    theories = ['mm', 'ani2x', 'mace_s', 'mace_m_fp32']
else:
    theories = ['mm', 'ani2x', 'mace_s']
for theory in theories:
    water_5_reps = load_repetition_properties(species, theory)
    print(f"Properties from repitition NPT runs for {species} with {theory}")
    print(water_5_reps)
    print("-"*50)
    print("\n")

if species!="nma":
    for theory in ['mm', 'ani2x', 'mace_s']:
        msd, diffusion = load_msd_diffusion(species, theory)
        print(f"Diffusion coefficient for {theory} {species}")
        print(diffusion.columns[0])
        print("-"*50)
else:
    print("No diffusion data for NMA available!")

print("\n")
