# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 13:08:16 2020

@author: Gaétan
"""

# =============================================================================
# Import des bibliothèques
# =============================================================================

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib


#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('font', size=14)
#plt.rc('axes',titlesize=16)
#plt.rc('legend',fontsize=14)
#plt.rc('figure',titlesize=20)
# =============================================================================
# Data importation
# =============================================================================

filename = 'TraitementDataPDF/GENERAL_u_mom.11mar2020'
T = np.genfromtxt(filename,skip_header=2)
y = T[:,0]*1e-2 #mm
mean_U = T[:,1]
mean_U2 = T[:,2]
mean_U3 = T[:,3]
mean_U4 = T[:,4]
u_rms = T[:,5]
u_skew = T[:,6]
u_kurt = T[:,7]

# =============================================================================
# Mean velocity profile
# =============================================================================

# Log fitting
def cost_fit(x):
    mean_u_fit = x[0]*np.log(1e-3*y[2:-4]) + x[1]
    return np.sqrt(np.mean(np.square(mean_u_fit-mean_U[2:-4]))/np.mean(np.square(mean_U[2:-4])))
res = scipy.optimize.minimize(cost_fit,[1,1])
x = res.x
yfit = np.linspace(np.min(y),np.max(y),1000)
mean_u_fit = x[0]*np.log(yfit*1e-3) + x[1]


def cost_lin(alpha):
    mean_u_lin = alpha*y[:3]
    return np.sqrt(np.mean(np.square((mean_u_lin - mean_U[:3])))/np.mean(np.square(mean_U[:3])))

res_lin = scipy.optimize.minimize(cost_lin,[1.])
alpha = res_lin.x[0]


# Log-law fitting
k = 0.41
a = 5.2
ustar = x[0]*k
nu = ustar*np.exp(-k*(x[1]/ustar-a))

nu_theo = 1.5e-5
y_5 = nu*ustar*5*1e3 # y (mm) @ which y+ = 5
y_30 = nu_theo*ustar*30*1e3 # y (mm) @ wjich y+ = 30


#ustar = np.sqrt(alpha*nu_theo)

print('Log-law fitting with coefficients (%.3e,%.3e)'%(x[0],x[1]))
print('Normalised residual RMS error : %.3e' % cost_fit(x))

h = 30.

plt.figure()
plt.plot(yfit/h,mean_u_fit)
plt.plot(y/h,mean_U,linestyle='dotted', marker='o',markersize=5.)
plt.xlabel('Distance from wall $\\eta = y/h$ - log scale')
plt.ylabel('Mean velocity (m/s)')
plt.xscale('log')
#plt.yscale('log')
plt.vlines(y_30/h,np.min(mean_U),np.max(mean_U),linestyle='dashed')
plt.text(y_30/h*1.2, 4., '$y^{+} = 30$')
plt.vlines(0.3,np.min(mean_U),np.max(mean_U),linestyle='dashed')
plt.text(0.34, 3.5, '$\\eta = \\frac{y}{h} = 0.3$')
plt.tight_layout()
plt.savefig('Mean_Velocity_Profile.pgf')
#plt.close()
# A comparer avec log law cf. Pope p.274

# =============================================================================
# RMS profile
# =============================================================================
plt.figure()
plt.plot(y,u_rms,linestyle='dotted', marker='o',markersize=5.)
plt.xlabel('Distance from wall $y$ - log scale (mm)')
plt.ylabel('$u_{RMS}$ (m/s)') 
plt.xscale('log')
plt.tight_layout()
plt.savefig('U_RMS_Profile.pgf')
plt.close()
# =============================================================================
# Skewness profile
# =============================================================================
plt.figure()
plt.plot(y,u_skew,linestyle='dotted', marker='o',markersize=5.)
plt.xlabel('Distance from wall $y$ - log scale (mm)')
plt.ylabel('Skewness $S_u$')
plt.xscale('log')
plt.tight_layout()
plt.savefig('Skewness_Profile.pgf')
plt.close()