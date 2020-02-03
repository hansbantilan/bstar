import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
import matplotlib 
matplotlib.use("agg")
import os 
import matplotlib.pyplot as plt
from Shooting_method import * 

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

verbose = 10
eps = 1e-10 # Small epsilon to avoid r \neq 0 

# ====================================
#   Main routine 
# ====================================

pewpew = Complex_Boson_Star(edelta_guess,phi0,D,Lambda,verbose) 

alpha0 = pewpew.radial_walker(Rstart,Rend,deltaR,N,eps)

#=====================================
#   Output and plotting
#=====================================
r = np.linspace(eps, Rend, N)  
sol = pewpew.get_solution()

omega, sol =  pewpew.normalise_edelta(sol) 

plt.plot(r, 1/sol[:, 0], 'b', label='edelta(t)')
plt.plot(r, 1+sol[:, 1], 'g', label='1+m(t)')
plt.plot(r, 1+sol[:, 2], 'r', label='1+phi(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylim([0.99,max(1+sol[:,2])*1.2])
plt.grid()

# ===============================
path = pewpew.get_path()

plt.savefig(path+"/overview.png")
np.savetxt(path+"/omega.dat",[omega]),
np.savetxt(path+"/rvals.dat",r),
np.savetxt(path+"/edelta.dat",1/sol[:, 0]),
np.savetxt(path+"/m.dat",sol[:, 1]),
np.savetxt(path+"/phi.dat",sol[:, 2]),

