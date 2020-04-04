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
    mean_u_fit = x[0]*np.log10(y) + x[1]
    return np.sqrt(np.mean(np.square(mean_u_fit-mean_U))/np.mean(np.square(mean_U)))
res = scipy.optimize.minimize(cost_fit,[1,1])
x = res.x
yfit = np.linspace(np.min(y),np.max(y),1000)
mean_u_fit = x[0]*np.log10(yfit) + x[1]

print('Log-law fitting with coefficients (%.3e,%.3e)'%(x[0],x[1]))
print('Normalised residual RMS error : %.3e' % cost_fit(x))

plt.figure()
plt.plot(yfit,mean_u_fit)
plt.plot(y,mean_U,linestyle='dotted', marker='o',markersize=5.)
plt.xlabel('Distance from wall $y$ - log scale (mm)')
plt.ylabel('Mean velocity (m/s)')
plt.xscale('log')
plt.tight_layout()
plt.savefig('Mean_Velocity_Profile.pgf')

# A comparer avec log law cf. Pope p.274

# =============================================================================
# RMS profile
# =============================================================================
plt.figure()
plt.plot(y,u_rms,linestyle='dotted', marker='o',markersize=5.)
plt.xlabel('Distance from wall $y$ (mm)')
plt.ylabel('u RMS (m/s)') 
plt.xscale('log')
plt.tight_layout()
plt.savefig('U_RMS_Profile.pgf')

# =============================================================================
# Skewness profile
# =============================================================================
plt.figure()
plt.plot(y,u_skew,linestyle='dotted', marker='o',markersize=5.)
plt.xlabel('Distance from wall $y$ (mm)')
plt.ylabel('Skewness $S_u$')
plt.xscale('log')
plt.tight_layout()
plt.savefig('Skewness_Profile.pgf')
