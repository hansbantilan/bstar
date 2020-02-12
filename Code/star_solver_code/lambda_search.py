from Bstar.ComplexBosonStar import *
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
import matplotlib

# =====================
#  All imporntnat definitions
# =====================

# Physics defintions
phi0 = 0.40         # centeral phi
D = 5.0             # Dimension (total not only spacial)
# Solver definitions
Rstart = 3
Rend = 50.00
deltaR = 1
N = 100000
#Lambda = -0.2       # Cosmological constant
#e_pow_minus_delta_guess = 0.4999
Lambda = -5.5       # Cosmological constant
e_pow_minus_delta_guess = 0.1912

verbose = 2
eps = 1e-10  # Small epsilon to avoid r \neq 0

Lstart =-5.6 #first cosmological constant Lambda
Lend=-6.0    #target cosmological constant Lambda
LN=5         #number of steps in Lambda

# ====================================
#   Main routine
# ====================================

lambda_list=np.linspace(Lstart,Lend,LN)

for Lambda in lambda_list:

	pewpew = Complex_Boson_Star(e_pow_minus_delta_guess, phi0, D, Lambda, verbose)
	pewpew.print_parameters()

	e_pow_minus_delta_guess = pewpew.radial_walker(Rstart, Rend, deltaR, N, eps)

	r, sol = pewpew.get_solution()
	pewpew.normalise_edelta()
	path = pewpew.get_path()
	pewpew.plot_solution()
	pewpew.print_solution()
	
	print "For next iteration, now using\n e_pow_minus_delta_guess=",e_pow_minus_delta_guess,"\n Rstart=",Rstart,"\n Rend=",Rend
