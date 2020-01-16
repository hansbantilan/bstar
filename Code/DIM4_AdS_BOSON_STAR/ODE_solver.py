import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
import matplotlib 
matplotlib.use("agg")
import os 
import matplotlib.pyplot as plt

def eqns(y, r, D, Lambda):
    """ Differential equation for scalar fields 

    Parameters:
        y (list with reals): current status vector ( a(r), alpha(r), phi(r), pi(r) ) 
	r (real) : current position 

    Returns:
        dydr (list with reals): derivative for y in r 

    """
    D = float(D)
    edelta, m, phi, pi = y
    # Where edelta  = e^{-\delta}

    F  = ( 1 - 2*m /r**( D-3 ) - 2 * Lambda * r**2 / (( D - 2 )*( D - 1 )) )
   
   

    dedeltadr = r*(edelta*pi**2.0 + edelta**(-1) *phi**2/F**2)
    dmdr      = r**( D - 2 ) * 0.5 * ( F * pi**2 + phi**2 + edelta**(-2) * phi**2 / F) 
    dphidr    = pi 
    
    dFdr      =  (-4*Lambda*r)/((-2 + D)*(-1 + D)) - 2*(3 - D)*r**(2 - D)*m - 2*r**(3 - D)*dmdr

    dpidr     =  -(phi/(edelta**2*F**2)) + phi/F - (dedeltadr*pi)/edelta - (dFdr*pi)/F + (2*pi)/r - (D*pi)/r
    dydr = [dedeltadr,dmdr,dphidr,dpidr]
    return dydr

# solve

def shoot(edelta_guess,phi0,D,Lambda,r):
    """ Solved differential equation for shooting process.

    Parameters:
        edelta_guess (real): The lapse value guess at r = rmin 
        D (int) :   Dimension of the problem
        Lambda (real) : cosmological constant 
	phi0 (real) : The phi value at r = rmin

    Returns:
        phi_end (real):. The phi value at r = rmax    

    """
    
    # Define initial data vector 
    y0 = [edelta_guess,0,phi0,0]
    # Solve differential equaion 
    sol = spi.odeint(eqns, y0, r, args = (D,Lambda))
    phi_end = sol[-1,2]	
    
    return phi_end

def radial_walker(edelta_guess,phi0,D,Lambda,rstart,rend,deltaR,N,eps): 
    """ Performs shooting for multiple radii rmax shooting process.

    Parameters:
        edelta_guess (real) : alpha guess for rmin calculation 
        phi0 (real) : phi value at r = 0 `
        D (int) :   Dimension of the problem
        Lambda (real) : cosmological constant 
        rstart (real) : first rmax for which shooting is performed
	rend (real) : maximum rmax for which shooting is performed
	deltaR (real) : stelpsize
	N (real) : number of gridpoints 

    Returns:
        alpha0 (real):. alpha0 for rmax   
    """
    range_list = np.arange(rstart,rend,deltaR)
    alpha0 = edelta_guess

    for R in range_list:
        r = np.linspace(eps, R, N)

        fun = lambda x: shoot(x,phi0,D,Lambda,r)
        root = opi.root(fun,alpha0)
        alpha0 = root.x 

        print("step ",R)
        print("alpha0 ",alpha0)
    
    return alpha0[0]

def normalise_edelta(sol): 
    """ Extractsomega for edelta by the coordinate transformation  t -> omega t 

    Parameters:
        sol (real array) : were the sol[:,1] corresponds to edelta^(-1) and
                           and asymtotic value that does not go to 1 
    Returns:
        omega (real): frequency of scalar field 
        sol (real array) : sol array with fixed edelta
    """
    edelta = 1./sol[:,0]
    N = len(edelta)
    omega = edelta[N-1]
    edelta = edelta/omega
    sol[:,0] = 1./edelta
    return omega , sol 

def make_filesystem(phi0,D,Lambda):
    """ Creates Folder for current physics problem if they do not yet exist
        
    Parameters:
        phi0 (real) : phi value at r = 0 `
        D (int) :   Dimension of the problem
        Lambda (real) : cosmological constant 
    Returns:
        path (string) : path to current file

    """

    name_Lambda_D = "Lambda"+str(Lambda)+"D"+str(D)
    if not os.path.exists(name_Lambda_D):
        os.mkdir(name_Lambda_D)

    name_phi = "phi"+str(phi0)
    if not os.path.exists(name_Lambda_D + "/"+name_phi):
        os.mkdir(name_Lambda_D + "/"+name_phi)
    
    path = name_Lambda_D + "/"+name_phi

    return path 

#=====================
#  All imporntnat definitions 
#=====================

# Physics defintions 
phi0 = 0.04         # centeral phi
D = 5.0             # Dimension (total not only spacial)
Lambda = -0.1       # Cosmological constant  
# Solver definitions 
Rstart = 3
Rend = 50.00
deltaR = 1
N = 100000
edelta_guess = 0.75

eps = 1e-10 # Small epsilon to avoid r \neq 0 
 
# ====================================
#   Main routine 
# ====================================

alpha0 = radial_walker(edelta_guess,phi0,D,Lambda,Rstart,Rend,deltaR,N,eps)

#=====================================
#   Output and plotting
#=====================================

r = np.linspace(eps, Rend, N)
y0 = [ alpha0, 0 ,phi0,0]
sol = spi.odeint(eqns, y0, r,args = (D,Lambda))

omega,sol =  normalise_edelta(sol) 

plt.plot(r, 1/sol[:, 0], 'b', label='edelta(t)')
plt.plot(r, 1+sol[:, 1], 'g', label='1+m(t)')
plt.plot(r, 1+sol[:, 2], 'r', label='1+phi(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylim([0.99,max(1+sol[:,2])*1.2])
plt.grid()

# ===============================
path = make_filesystem(phi0,D,Lambda)

plt.savefig(path+"/overview.png")
np.savetxt(path+"/omega.dat",[omega]),
np.savetxt(path+"/rvals.dat",r),
np.savetxt(path+"/edelta.dat",1/sol[:, 0]),
np.savetxt(path+"/m.dat",sol[:, 1]),
np.savetxt(path+"/phi.dat",sol[:, 2]),

