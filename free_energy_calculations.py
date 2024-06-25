import numpy as np

def calculate_free_energy(T, m, Iabc, sigma, M, nu, Edft, Linearity, p):
    """
    Calculate the free energy of a molecule.

    Parameters:
    T (float): Temperature in Kelvin
    m (float): Mass of the molecule in kg
    Iabc (float or list): Moment of inertia (kg·m²) or list of moments of inertia for non-linear molecules
    sigma (int): Symmetry number
    M (int): Spin multiplicity
    nu (array): Array of vibrational frequencies in cm^-1
    Edft (float): Electronic energy from DFT calculations in eV
    Linearity (bool): True if the molecule is linear, False otherwise
    p (float): Pressure in Pa

    Returns:
    float: Total free energy in eV
    """
    kB = 8.617333262145e-5  # Boltzmann constant in eV/K   
    kB_joules = 1.380649E-23  # Boltzmann constant in J/K
    h_joules = 6.62607015e-34  # Planck's constant in J·s (same value as eV·s)
    h = 4.135667696e-15  # Planck's constant in eV·s
    
    # Convert vibrational frequencies from cm^-1 to Hz
    nu_hz = nu * 2.99792458e10  # Speed of light in cm/s
    
    # Vibrational free energy
    F_vib = np.sum((h * nu_hz) / 2) + np.sum(kB * T * np.log(1 - np.exp(-h * nu_hz / (kB * T))))
    
    # Translational free energy (Gibbs free energy version)
    F_trans = -kB * T * np.log((2 * np.pi * m * kB * T / h**2) ** (3/2) * kB * T / p)
    
    # Rotational free energy
    if Linearity:
        F_rot_joules = -kB_joules * T * np.log(8 * np.pi**2 * Iabc * kB_joules * T / h_joules**2)
    else:
        F_rot_joules = -kB_joules * T * np.log(8 * np.pi**2 * (2 * np.pi * kB_joules * T / h_joules**2)**(3/2)) - (kB_joules * T * np.log(Iabc) / 2)
    
    # Convert rotational free energy from joules to eV
    F_rot = F_rot_joules * 6.241509074e18  # Conversion factor from J to eV
    
    # Symmetry free energy
    F_symm = kB * T * np.log(sigma)
    
    # Spin free energy
    F_spin = -kB * T * np.log(M)
    
    # Total free energy
    F_total = F_trans + F_rot + F_vib + F_symm + F_spin + Edft
    
    return F_total

def calculate_reaction_free_energy(T, p, reactants, products):
    """
    Calculate the Gibbs free energy change for a reaction.

    Parameters:
    T (float): Temperature in Kelvin
    p (float): Pressure in Pa
    reactants (list): List of tuples with parameters for reactant molecules
    products (list): List of tuples with parameters for product molecules

    Returns:
    float: Gibbs free energy change for the reaction in eV
    """
    F_reactants = sum(calculate_free_energy(T, *params, p) for params in reactants)
    F_products = sum(calculate_free_energy(T, *params, p) for params in products)
    ΔG = F_products - F_reactants
    
    return ΔG

# Example reaction data
reactions = [
    ([  # Reactants for reaction 1 (MgSiO3 + CH')
        (1.6669934e-25, 8.9416055e-135, 2, 1, np.array([165.00471544, 300.31341353, 423.49748267, 467.42456698, 602.34483041, 688.2404795, 808.27235342, 853.26287239, 1275.33240978]), -19503.223762886, False),
        (2.1617958e-26, 1.9550563e-47, 1, 2, np.array([2830.57687704]), -1046.613540389, True)
    ], [
        (1.4000651e-25, 3.6772064e-135, 2, 1, np.array([167.33986661, 492.6030138, 533.09967848, 688.15369804, 819.76849106, 829.86371134, 943.0655396, 1054.27974773, 1821.47201843]), -12636.738916183, False),
        (4.8310788e-26, 3.8005954e-47, 1, 2, np.array([1984.22718546]), -7914.103305082, True)
    ])
]

# Example usage
T = 298.15  # Temperature in Kelvin
p = 101325  # Pressure in Pa

for reactants, products in reactions:
    ΔG = calculate_reaction_free_energy(T, p, reactants, products)
    print(f"ΔG for the reaction: {ΔG:.4f} eV")
