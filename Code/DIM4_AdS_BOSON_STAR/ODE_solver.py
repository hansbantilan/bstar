import numpy as np
import scipy.integrate as spi
import scipy.optimize as opi
import matplotlib 
matplotlib.use("agg")
import os 
import matplotlib.pyplot as plt

def eqns(y, r):
    """ Differential equation for scalar fields 

    Parameters:
        y (list with reals): current status vector ( a(r), alpha(r), phi(r), pi(r) ) 
	r (real) : current position 

    Returns:
        dydr (list with reals): derivative for y in r 

    """
    D = 5.0 
    Lambda = -0.1
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

def shoot(alpha0_guess,phi0,r):
    """ Solved differential equation for shooting process.

    Parameters:
        alpha0_guess (real): The lapse value guess at r = rmin 
	phi0 (real) : The phi value at r = rmin

    Returns:
        phi_end (real):. The phi value at r = rmax    

    """
    
    # Define initial data vector 
    y0 = [alpha0_guess,0,phi0,0]
    # Solve differential equaion 
    sol = spi.odeint(eqns, y0, r)
    phi_end = sol[-1,2]	
    
    return phi_end

def radial_walker(alpha0_guess,phi0,rstart,rend,deltaR,N): 
    """ Performs shooting for multiple radii rmax shooting process.

    Parameters:
        alpha0_guess (real) : alpha guess for rmin calculation 
        phi0 (real) : phi value at r = 0 
        rstart (real) : first rmax for which shooting is performed
	rend (real) : maximum rmax for which shooting is performed
	deltaR (real) : stelpsize
	N (real) : number of gridpoints 

    Returns:
        alpha0 (real):. alpha0 for rmax   
    """

    eps = 1e-10 # distance from zero
    range_list = np.arange(rstart,rend,deltaR)
    alpha0 = alpha0_guess

    for R in range_list:
        r = np.linspace(eps, R, N)

        fun = lambda x: shoot(x,phi0,r)
        root = opi.root(fun,alpha0)
        alpha0 = root.x 

        print("step ",R)
        print("alpha0 ",alpha0)
    
    return alpha0[0]

# We want to evaluate the system on 30 linearly
# spaced times between t=0 and t=3.


phi0 = 0.1
# Resolution of diff eqn 
Rstart = 6
Rend = 30.00
deltaR = 0.5
N = 100000

alpha0 = radial_walker(1,phi0,Rstart,Rend,deltaR,N)

name = "phi"+str(phi0)
if not os.path.exists(name):
    os.mkdir(name)


r = np.linspace(1e-10, Rend, N)
y0 = [ alpha0, 0 ,phi0,0]
sol = spi.odeint(eqns, y0, r)

plt.plot(r, 1/sol[:, 0], 'b', label='edelta(t)')
plt.plot(r, 1+sol[:, 1], 'g', label='1+m(t)')
plt.plot(r, 1+sol[:, 2], 'r', label='1+phi(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylim([0.99,max(1+sol[:,2])*1.5])
plt.grid()

plt.savefig(name+"/overview.png")

np.savetxt(name+"/edelta.dat",1/sol[:, 0]),
np.savetxt(name+"/m.dat",sol[:, 1]),
np.savetxt(name+"/phi.dat",sol[:, 2]),

