#!/usr/bin/env python3
from math import pi

# Constants
h = 6.63e-34            # Planck's Constant
pi = pi                 # π

# Engineering maths
epsilon0 = 8.85e-12     # Free space permittivity
mu0 = 4 * pi * 10**-7   # Free space permeability 4*pi*10^-7
kCoulomb = 9e9          # Coulomb's Constant 1/(4*pi*epsilon0)
c = 3e8                 # Speed of Light in air/vacuum 1/(mu0*epsilon0)^0.5

# Info and System
VT = 26e-3              # Thermal Voltage when T = 300 K (room temperature)
kBolzman = 1.38e-23

# Power and Energy
e = -1.6e-19            # The charge of an electron
rowAir = 1.2            # Air density at sea level
Cpmax = 16 / 27         # Maximum power coefficient
