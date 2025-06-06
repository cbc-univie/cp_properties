import numpy as np
from numpy.typing import NDArray
from openff.units import unit, Quantity
from typing import cast, Callable
from collections import namedtuple
import pandas as pd
from openff.units import unit




##################### DEFINE PROPERTY FUNCTIONS ####################################


Constants = namedtuple('Constants', ['GAS_CONSTANT', 'BOLTZMANN_CONSTANT'])

CONSTANTS = Constants(
    GAS_CONSTANT=8.31446261815324 * unit.joule / unit.mole / unit.kelvin,
    BOLTZMANN_CONSTANT=1.380649e-23 * unit.joule / unit.kelvin,
)


def calc_heat_capacity_units(
    total_energy: NDArray[np.float64],
    number_particles: int,
    temp: float,
    molar_mass: float,
    printing: bool
) -> Quantity:
    """
    Compute the heat capacity.

    C_p = Variance(Total energy) / (number_particles * temp^2 * gas constant)
    the return value is in units cal/mole/kelvin
    for a correct unit transformation, the molar mass is required
    """
    pot_var = total_energy.var() * (unit.kilojoule / unit.mole) ** 2
    temp = temp * unit.kelvin

    val = cast(Quantity,
               pot_var
               / number_particles
               / temp**2
               / CONSTANTS.GAS_CONSTANT
               )
    val = val.to(unit.cal / unit.mole / unit.kelvin) / molar_mass
    if printing:
        print("heat capacity: ", val)
    return val


def calc_thermal_expansion(
        total_energy: NDArray[np.float64],
        volume: NDArray[np.float64],
        temp: float,
        printing: bool
) -> Quantity:
    """
    Compute the coefficient of thermal expansion.

    alpha = Cov(energy, vol) / (vol * temp^2 * gas constant)
    the return value is in units 1/Kelvin
    """
    cov_en_vol = cast(Quantity,
                      np.cov(total_energy, volume)[0][1]
                      * (unit.nanometer**3)
                      * unit.kJ / unit.mole
                      )

    T = temp * unit.kelvin
    volume = volume.mean() * (unit.nanometer**3)

    alpha = cov_en_vol / CONSTANTS.GAS_CONSTANT / T**2 / volume
    alpha_shift = cast(Quantity, alpha.to(1 / unit.kelvin))
    if printing:
        print("thermal expansion: ", alpha_shift)
    return alpha_shift


def calc_isothermal_compressibility(
        volume: NDArray[np.float64],
        temp: float,
        printing: bool
) -> Quantity:
    """
    Compute the isothermal compressibility.

    kappa = Variance(Box volume) / (k_B * temperature * volume)
    the return value is in units 1/bar
    """
    volume_var = volume.var() * (unit.nanometer**3) ** 2
    volume_mean = volume.mean() * (unit.nanometer**3)
    T = temp * unit.kelvin

    val = volume_var / CONSTANTS.BOLTZMANN_CONSTANT / T / volume_mean
    val = val.to(1 / unit.bar)
    if printing:
        print("thermal expansion: ", val)
    return val


def calc_heat_of_vaporization(
    pot_energy: NDArray[np.float64],
    pot_energy_mono: NDArray[np.float64],
    temp_traj: NDArray[np.float64],
    box_count: int,
    printing: bool
) -> Quantity:
    """
    Compute the heat of vaporization.

    Delta H_vap = mean_energy_gas - mean_energy_liquid + R*temperature
    the return value is in units kJ
    """
    pot_mean = pot_energy.mean() * unit.kilojoule / unit.mole / box_count
    pot_mono_mean = pot_energy_mono.mean() * unit.kilojoule / unit.mole
    temp_mean = temp_traj.mean() * unit.kelvin
    val = pot_mono_mean - pot_mean + CONSTANTS.GAS_CONSTANT * temp_mean
    if printing:
        print("heat of vaporozation: ", val)
    return val


def bootstrap_hov(
        liquid_pot: NDArray[np.float64],
        mono_pot: NDArray[np.float64],
        liquid_temp: NDArray[np.float64],
        box_count: int,
        Nboot: int,
        statfun: Callable[
            [NDArray[np.float64],
             NDArray[np.float64],
             NDArray[np.float64],
             int,
             bool], Quantity
             ]
) -> NDArray[np.float64]:
    """Calculate bootstrap statistics for a sample x."""
    liquid_pot = np.array(liquid_pot)
    liquid_temp = np.array(liquid_temp)
    mono_pot = np.array(mono_pot)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0, len(liquid_pot), len(liquid_pot))
        sample_pot = liquid_pot[index]
        sample_temp = liquid_temp[index]
        bastatistics = (
            statfun(sample_pot,
                    mono_pot,
                    sample_temp,
                    box_count,
                    False).magnitude
        )
        resampled_stat.append(bastatistics)

    return np.array(resampled_stat)


def bootstrap_hcap(
        liquid_total: NDArray[np.float64],
        box_count: int,
        liquid_temp: NDArray[np.float64],
        molar_mass: float,
        Nboot: int,
        statfun: Callable[
            [NDArray[np.float64],
             int,
             float,
             float,
             bool], Quantity
             ]
) -> NDArray[np.float64]:
    """Calculate bootstrap statistics for a sample x."""
    liquid_total = np.array(liquid_total)
    liquid_temp = np.array(liquid_temp)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0, len(liquid_total), len(liquid_total))
        sample_pot = liquid_total[index]
        sample_temp = liquid_temp[index]
        bastatistics = (
            statfun(sample_pot,
                    box_count,
                    sample_temp.mean(),
                    molar_mass,
                    False).magnitude
        )
        resampled_stat.append(bastatistics)

    return np.array(resampled_stat)


def bootstrap_texp(
        liquid_total: NDArray[np.float64],
        box_vol: NDArray[np.float64],
        liquid_temp: NDArray[np.float64],
        Nboot: int,
        statfun: Callable[
            [NDArray[np.float64],
             NDArray[np.float64],
             float,
             bool], Quantity
             ]
) -> NDArray[np.float64]:
    """Calculate bootstrap statistics for a sample x."""
    liquid_total = np.array(liquid_total)
    box_vol = np.array(box_vol)
    liquid_temp = np.array(liquid_temp)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0, len(liquid_total), len(liquid_total))
        sample_pot = liquid_total[index]
        sample_box_vol = box_vol[index]
        sample_temp = liquid_temp[index]
        bastatistics = (
            statfun(sample_pot,
                    sample_box_vol,
                    sample_temp.mean(),
                    False).magnitude
        )
        resampled_stat.append(bastatistics)

    return np.array(resampled_stat)


def bootstrap_icomp(
        box_vol: NDArray[np.float64],
        liquid_temp: NDArray[np.float64],
        Nboot: int,
        statfun: Callable[
            [NDArray[np.float64],
             float,
             bool], Quantity
             ]
) -> NDArray[np.float64]:
    """Calculate bootstrap statistics for a sample x."""
    box_vol = np.array(box_vol)
    liquid_temp = np.array(liquid_temp)

    resampled_stat = []
    for k in range(Nboot):
        index = np.random.randint(0, len(box_vol), len(box_vol))
        sample_box_vol = box_vol[index]
        sample_temp = liquid_temp[index]
        icomp = statfun(
            sample_box_vol,
            sample_temp.mean(),
            False
            ).magnitude
        resampled_stat.append(icomp)

    return np.array(resampled_stat)






############## COMPUTE PROPERTIES FROM WATER DATA #############################




def load_csv(filename: str, method: str) -> pd.DataFrame:
    data_path = f"data/water_traj/{method}/{filename}"
    return pd.read_csv(data_path, sep="\t")
    

def get_liquid_traj(method: str, start: int, end: int):
    csv_complete=pd.DataFrame()
    for i in range(start,end+1): 
        csv_part = load_csv(f"{method}_tip572_{i}_NPT.csv", method)
        csv_complete = pd.concat([csv_complete, csv_part])
    return csv_complete

def get_gas_traj(method: str):
    csv_complete=load_csv(f"gas_{method}.csv", method)
    return csv_complete



skip_size = 0.090909
box_count=572
molar_mass = 18.015 * unit.gram / unit.mole

theories = ['mm', 'ani2x', 'mace_s', 'mace_m', 'mace_m_fp32']

print("\nCOMPUTING CONDENSED PHASE PROPERTIES\n")
for theory in theories:
    liquid = get_liquid_traj(theory, 1, 11)
    gas = get_gas_traj(theory)

    skip_part_gas = int(round(gas["Potential Energy (kJ/mole)"].count()*skip_size,0))
    gas_cut = gas[skip_part_gas-1:-1] # skip the first 10%

    skip_part_liquid = int(round(liquid["Potential Energy (kJ/mole)"].count()*skip_size,0))
    liquid_cut = liquid[skip_part_liquid-1:-1] # skip the first 10%
    
    heat_capacity = calc_heat_capacity_units(
        liquid_cut["Total Energy (kJ/mole)"].to_numpy(), 
        box_count, 
        liquid_cut["Temperature (K)"].mean(), 
        molar_mass, 
        False
    )
    
    thermal_expansion = calc_thermal_expansion(
        liquid_cut["Total Energy (kJ/mole)"].to_numpy(), 
        liquid_cut["Box Volume (nm^3)"].to_numpy(), 
        liquid_cut["Temperature (K)"].mean(), 
        False
    )
    
    iso_comp = calc_isothermal_compressibility(
        liquid_cut["Box Volume (nm^3)"].to_numpy(), 
        liquid_cut["Temperature (K)"].mean(), 
        False
    )
    hov = calc_heat_of_vaporization (
        liquid_cut["Potential Energy (kJ/mole)"].to_numpy(), 
        gas_cut["Potential Energy (kJ/mole)"].to_numpy(), 
        liquid_cut["Temperature (K)"].to_numpy(), 
        box_count, 
        False
    )

    density = liquid_cut["Density (g/mL)"].mean()
    
    print(f"{'Property':35} {'Model':10} {'Value':>10} {'Unit':>25}")
    print("-" * 85)
    print(f"{'Heat capacity':35} {theory:<10} {round(heat_capacity.magnitude, 2):>10} {str(heat_capacity.units):>25}")
    print(f"{'Thermal expansion (*1e2)':35} {theory:<10} {round(thermal_expansion.magnitude*1e2, 2):>10} {str(thermal_expansion.units):>25}")
    print(f"{'Isothermal compressibility (*1e4)':35} {theory:<10} {round(iso_comp.magnitude*1e4, 2):>10} {str(iso_comp.units):>25}")
    print(f"{'Heat of vaporization':35} {theory:<10} {round(hov.magnitude, 2):>10} {str(hov.units):>25}")
    print(f"{'Density':35} {theory:<10} {round(density, 2):>10} {str((unit.gram / unit.milliliter)):>25}")
    print("-" * 85)




########################### COMPUTE BOOTSTRAPPING SAMPLES ##################################################





def print_results(boot_hov_mean: float,
                  boot_hov_std: float,
                  boot_hcap_mean: float,
                  boot_hcap_std: float,
                  boot_texp_mean: float,
                  boot_texp_std: float,
                  boot_icomp_mean: float,
                  boot_icomp_std: float):
    print(f"{'Heat of vaporization - mean:':45} {boot_hov_mean:.6f}")
    print(f"{'Heat of vaporization - std:':45} {boot_hov_std:.6f}\n")
    
    print(f"{'Heat capacity - mean:':45} {boot_hcap_mean:.6f}")
    print(f"{'Heat capacity - std:':45} {boot_hcap_std:.6f}\n")
    
    print(f"{'Thermal expansion coeff. - mean:':45} {boot_texp_mean:.6f}")
    print(f"{'Thermal expansion coeff. - std:':45} {boot_texp_std:.6f}\n")
    
    print(f"{'Isothermal compressibility - mean:':45} {boot_icomp_mean:.6f}")
    print(f"{'Isothermal compressibility - std:':45} {boot_icomp_std:.6f}\n")
    print('-'*60)

skip_size = 0.090909
box_count=572
molar_mass = 18.015 * unit.gram / unit.mole
Nboot=100


theories = ['mm', 'ani2x', 'mace_s', 'mace_m', 'mace_m_fp32']

print("\nCOMPUTING BOOTSTRAPPING VALUES FOR CONDENSED PHASE PROPERTIES\n")
for theory in theories:
    liquid = get_liquid_traj(theory, 1, 11)
    gas = get_gas_traj(theory)

    skip_part_gas = int(round(gas["Potential Energy (kJ/mole)"].count()*skip_size,0))
    gas_cut = gas[skip_part_gas-1:-1] # skip the first 10%

    skip_part_liquid = int(round(liquid["Potential Energy (kJ/mole)"].count()*skip_size,0))
    liquid_cut = liquid[skip_part_liquid-1:-1] # skip the first 10%


    boot_hov = bootstrap_hov(liquid_pot=liquid_cut["Potential Energy (kJ/mole)"],
                                               mono_pot=gas_cut["Potential Energy (kJ/mole)"],
                                               liquid_temp=liquid_cut["Temperature (K)"],
                                               box_count=box_count,
                                               Nboot=Nboot,
                                               statfun=calc_heat_of_vaporization
                                               )

    boot_hcap = bootstrap_hcap(liquid_total=liquid_cut["Total Energy (kJ/mole)"],
                                                 box_count=box_count,
                                                 liquid_temp=liquid_cut["Temperature (K)"],
                                                 molar_mass=molar_mass,
                                                 Nboot=Nboot,
                                                 statfun=calc_heat_capacity_units
                                                 )

    boot_texp = bootstrap_texp(liquid_total=liquid_cut["Total Energy (kJ/mole)"],
                                                 box_vol=liquid_cut['Box Volume (nm^3)'],
                                                 liquid_temp=liquid_cut["Temperature (K)"],
                                                 Nboot=Nboot,
                                                 statfun=calc_thermal_expansion
                                                 )

    boot_icomp= bootstrap_icomp(box_vol=liquid_cut['Box Volume (nm^3)'],
                                                  liquid_temp=liquid_cut["Temperature (K)"],
                                                  Nboot=Nboot,
                                                  statfun=calc_isothermal_compressibility
                                                  )
    
    print(f'Results for {theory}:')
    print('-'*60)
    print_results(boot_hov.mean(),
                  boot_hov.std(),
                  boot_hcap.mean(),
                  boot_hcap.std(),
                  boot_texp.mean(),
                  boot_texp.std(),
                  boot_icomp.mean(),
                  boot_icomp.std()
                  )
    
    


