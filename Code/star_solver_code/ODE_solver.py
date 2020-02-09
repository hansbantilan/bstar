import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
import matplotlib 
matplotlib.use("agg")
import os 
import matplotlib.pyplot as plt
from bstar.ComplexBosonStar import * 

#=====================
#  All imporntnat definitions 
#=====================

# Physics defintions 
phi0 = 0.40         # centeral phi
D = 5.0             # Dimension (total not only spacial)
Lambda = -0.2       # Cosmological constant  
# Solver definitions 
Rstart = 3
Rend = 50.00
deltaR = 1
N = 100000
edelta_guess = 0.4999

verbose = 1
eps = 1e-10 # Small epsilon to avoid r \neq 0 

# ====================================
#   Main routine 
# ====================================

pewpew = Complex_Boson_Star(edelta_guess,phi0,D,Lambda,verbose) 

pewpew.print_parameters()

alpha0 = pewpew.radial_walker(Rstart,Rend,deltaR,N,eps)


#=====================================
#   Output and plotting
#=====================================
r,sol = pewpew.get_solution()

omega, sol =  pewpew.normalise_edelta(sol) 

# ===============================
path = pewpew.get_path()
pewpew.plot_solution()

plt.savefig(path+"/overview.png")
np.savetxt(path+"/omega.dat",[omega]),
np.savetxt(path+"/rvals.dat",r),
np.savetxt(path+"/edelta.dat",1/sol[:, 0]),
np.savetxt(path+"/m.dat",sol[:, 1]),
np.savetxt(path+"/phi.dat",sol[:, 2]),

