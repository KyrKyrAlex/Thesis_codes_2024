# Thesis_codes_2024
This is a repository that containts a python code that I created to evaluate the favourability of reactions based on the Gibbs free energy difference ΔG

## Introduction
To evaluate the favorability of reactions, we employed thermodynamics, specifically statistical mechanics. We calculated the Gibbs free energy change (ΔG) for the reactions concerning temperatures and pressures specified in our adopted protoplanetary disk model.

## Thermodynamics in Low-Density Environments
A point of debate is whether thermodynamics is suitable for an environment where densities are low and chemistry is kinetically dominated. While this can be considered an approximation in our study, it's important to note that this research serves as an initial step in evaluating potential interactions or reactions between magnesium silicate clusters and carbon-containing molecules. Protoplanetary disks are not typical regions in space; they exhibit higher densities and unusual phenomena such as accretion and turbulence. For these reasons, we believe statistical mechanics is a suitable initial approximation for evaluating the viability of the reactions in our study.

## Gibbs Free Energy Calculation
The main feature we calculate is the difference in Gibbs free energy between the reactants and the products of the reactions. When this difference is below zero, the reaction will occur spontaneously under those conditions. The Gibbs free energy difference can be written as:

$$
\Delta G = \Delta H - T\Delta S
$$

To calculate this difference, we need to compute the molecular partition function of the molecules. We are working under the ideal gas law and the canonical ensemble (N, V, T), where the number of molecules, the volume, and the temperature of our system remain constant. The partition function for a many-particle system is defined as:

$$
Q = \sum_{s} e^{-\beta E_s}
$$

where $\beta = (kT)^{-1}$, $k$ is the Boltzmann constant, $s$ is the number of microstates that the system possesses, and $E_s$ is the energy of each state. This creates an ensemble of the possible states our particles can access, which can then be translated into macroscopic properties. The connection of free energy with the partition function is given by:

$$
G(T,V,N) = -K_bT \ln Q(T,V,N) + pV
$$

The partition function can be calculated by defining its different contributions: translational, rotational, vibrational, electronic, and configurational degrees of freedom. This makes the equation look like:

$$
G = -K_b T \ln (q^{trans} q^{rot} q^{vib} q^{elec}) + pV
$$

## Translational Partition Function
The translational energy levels are very closely spaced, thus, at normal temperatures, large numbers of them are typically accessible. For a system having macroscopic dimensions, the translational partition function is:

$$
q^{trans} = \left(\frac{2\pi m k_B T}{h^2}\right)^\frac{3}{2} V
$$

Including this in the Gibbs free energy expression, the translational contribution to the Gibbs free energy is:

$$
G^{trans} = -k_B T \ln \left[\left(\frac{2\pi m k_B T}{h^2}\right)^\frac{3}{2} \frac{k_B T}{P}\right] - k_B T
$$

## Rotational Partition Function
Under the rigid rotor approximation, the rotational partition function is:

$$
q^{rot} = \sum_{J=0}^{\infty} (2J+1)e^{-\frac{J(J+1)B_O}{K_b T}}
$$

where $B_O$ is the rotational constant:

$$
B_O = \frac{\hbar^2}{2I}
$$

At high temperatures, the rotational free energy for non-linear molecules is:

$$
G^{rot} = -k_B T \ln \left(8\pi^2 \left(\frac{2\pi k_B T}{h^2}\right)^\frac{3}{2}\right) - \frac{1}{2}k_B T \ln \left(I_a I_b I_c\right)
$$

And for linear molecules:

$$
G^{rot} = -k_B T \ln \left(\frac{8\pi^2 k_B T I_a}{h^2}\right)
$$

The symmetry contribution is expressed as:

$$
G^{symm} = k_B T \ln \sigma
$$

## Vibrational Partition Function
Under the harmonic oscillator approximation, the vibrational partition function is:

$$
q^{vib} = \sum_{i=0}^{M}\sum_{n=0}^{\infty} e^{-\frac{(n+\frac{1}{2})\hbar \omega_i}{k_B T}}
$$

The vibrational free energy is:

$$
G^{vib} = \sum_{i=1}^{M} \left[\frac{\hbar \omega_i}{2} + k_B T \ln \left[1 - e^{-\frac{\hbar \omega_i}{k_B T}}\right]\right]
$$

## Electronic Partition Function
For the electronic partition function, considering the possible spin degeneracy:

$$
G^{elec} = E^{el.struct} - k_B T \ln M
$$

where $E^{el.struct}$ is the total energy of the electronic structure calculation and $M$ is the multiplicity of the molecule.

## Total Free Energy
Combining all these contributions, the total free energy of a molecule in isolation is:

$$
G = G^{trans} + G^{rot} + G^{symm} + G^{vib} + G^{elec}
$$

The $pV$ term is canceled by the $-K_b T$ from the translational contribution due to the ideal gas law and molecular formulas.

## Conclusion
This document outlines the theoretical foundation and computational methodology used to evaluate the favorability of reactions within protoplanetary disks using statistical mechanics. The approximations and considerations made provide a robust framework for understanding these complex environments.

