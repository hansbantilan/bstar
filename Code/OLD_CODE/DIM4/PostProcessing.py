import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import sys

A_data = np.loadtxt("general_a202.dat")
C_data = np.loadtxt("general_c202.dat")
phi_data = np.loadtxt("general_phi202.dat")

dx = 0.35
deltax = phi_data[1,0]-phi_data[0,0]
red = int(dx/deltax)

phi_data = phi_data[::red]
A_data = A_data[::red]
C_data = C_data[::red]


N = len(A_data[0,:])
Len = len(A_data)
kappa = 8*np.pi;
deltaxA = (A_data[1,0]-A_data[0,0])
deltaxC = (C_data[1,0]-C_data[0,0])
deltax = phi_data[1,0]-phi_data[0,0]

if( (deltaxA != deltaxC) or (deltaxC != deltax)):
	sys.exit("Inconsistent Stepsize! Maybe mixed up files?")

print("DeltaX = ",  deltax)
print("Order of fouier = " , N)
print("Number of points = ", Len)
Omega = np.sqrt(C_data[Len-1,1]/A_data[Len-1,1])
print("Omega = ", Omega)

t = 0;#  np.pi/2.*1/Omega

#if( t != 0 ):
#	sys.exit("Not made for t not equal 0 ")


# ====================================
# Resummation 
# ====================================

r     = phi_data[:,0]
phi    = np.zeros(Len)
dphidt = np.zeros(Len)
A      = np.zeros(Len)
C      = np.zeros(Len)

for i in range(N-2):
	phi    +=   kappa**(-0.5) * phi_data[:,i+1] *               np.cos((i+1)*Omega*t)
	dphidt += - kappa**(-0.5) * phi_data[:,i+1] * (i+1)*Omega * np.sin((i+1)*Omega*t) 

for i in range(N-2):
	A  += A_data[:,i+1] * np.cos(i*Omega*t)
	C  += C_data[:,i+1] * np.cos(i*Omega*t)


alpha  = Omega*np.sqrt(A/C)

dphidr   = np.gradient(phi  )/deltax
dAdr     = np.gradient(A    )/deltax
dalphadr = np.gradient(alpha)/deltax


# ====================================
# ADM mass 
# ====================================

rr1 = 0.70*r[Len-1]
rr2 = r[Len-1]
print("Asymtotic values to calculate MADM")
print("r1 = ", rr1)
print("r2 = ", rr2)

N1 = int(rr1/deltax)
N2 = int(rr2/deltax)

A1 = A[N1]
A2 = A[N2]

MADM = -(A1*rr1**2*rr2**2 - A2*rr1**2*rr2**2)/(2*A1*A2*(rr1**2-rr2**2))

# Asymtotic value of grr

grrInf = - (A2*rr1**2 - A1*rr2**2)/(A1*A2*(rr2**2-rr1**2))

print("MADM = ", MADM)
print("Asymtotic grr ", grrInf) 


#=====================================
# Output
#=====================================

if(t == 0): 
	np.savetxt("grr001.csv", A)
	np.savetxt("Phi001.csv", phi)
	np.savetxt("alpha001.csv",alpha)


#=====================================
# Violation of Ham
#=====================================

Ttt = dphidt**2 - alpha**2*(-dphidr**2/(2.*A) + dphidt**2/(2.*alpha**2) - phi**2/2.)
#Ttr = dphidr*dphidt
Trr = dphidr**2 + A*(-dphidr**2/(2.*A) + dphidt**2/(2.*alpha**2) - phi**2/2.) 


Gtt   = (3*alpha**2*(-2*A + 2*A**2 + dAdr*r))/(2.*A**2*r**2)
#Gtr   = (3*dAdt)/(2.*A*r)
Grr   = (3*(-((-1 + A)*alpha) + dalphadr*r))/(alpha*r**2)


res = Gtt - 8 * np.pi * Ttt

rel_res = abs(Gtt/Ttt - 8*np.pi)*100


Ttt_max = max(Ttt)
Tttfun = interp1d(r, Ttt - 0.01*Ttt_max, kind="linear")
size = fsolve(Tttfun,max(r)/20)[0]	
print("99 % size of star = ",size)
#=====================================
# Plotting 
#=====================================

plt.figure(figsize=(14, 9), dpi=200)
plt.plot(r[1:],rel_res[1:],linewidth = 2.5,label = r'$G_tt-8 \pi T_tt$ at time '+str(t) )
plt.xlim(0,size)
plt.ylim(0,100)
plt.legend()
plt.savefig("rel_Ham_profile.png",bbox_inches = 'tight')
plt.close()

plt.figure(figsize=(14, 9), dpi=200)
plt.plot(r[1:],res[1:],linewidth = 5,label = r'$G_tt-8 \pi T_tt$ at time '+str(t) )
plt.legend()
plt.savefig("Ham_profile.png",bbox_inches = 'tight')
plt.close()


plt.figure(figsize=(14, 9), dpi=200)
plt.plot(r,Ttt,linewidth = 5,label = r'$\Phi$ at time '+str(t) )
plt.legend()
plt.savefig("Trr_profile.png",bbox_inches = 'tight')
plt.close()


plt.figure(figsize=(14, 9), dpi=200)
plt.plot(r,phi,linewidth = 5,label = r'$\Phi$ at time '+str(t) )
plt.xlabel(r'$r~[1/M_{pl}]$')
plt.ylabel(r'$r]$')
plt.legend()
plt.savefig("Phi_profile.png",bbox_inches = 'tight')
plt.close()

plt.figure(figsize=(14, 9), dpi=200)
plt.plot(r,A,linewidth = 5,label = r'$g_rr$ at time '+str(t) )
plt.plot(r[N1:N2],(grrInf-2*MADM/r[N1:N2]**2)**(-1) ,linewidth = 5,linestyle = '--', label = r'$g_rr$ of Schwarzschild at time '+str(t) )
plt.xlabel(r'$r~[1/M_{pl}]$')
plt.ylabel(r'$r]$')
plt.legend()
plt.savefig("grr_profile.png",bbox_inches = 'tight')
plt.close()


