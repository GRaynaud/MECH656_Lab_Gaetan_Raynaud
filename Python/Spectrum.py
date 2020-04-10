# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 10:20:34 2020

@author: Gaétan
"""

# =============================================================================
# Import des librairies
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import splrep,splev
import tikzplotlib

# =============================================================================
# Acquisition des données
# =============================================================================

filename = 'TraitementDataSpec/specrsltfilnam020.11mar2020.txt'

T = np.genfromtxt(filename,skip_header=1)

f = T[:,0]
Spec = T[:,1]

# =============================================================================
#  Conversion temporel --> spatial
# y = [  0.2 ,   5.5 ,  20.  ]
# mean_U = [2.3158,  4.1619, 5.0156]
# =============================================================================
nu = 1.5e-5 #m/s à vérifier
U1 = 2.3158 #m/s
K1 = 2*np.pi*f/U1
F11_K1 = Spec*U1/(2*np.pi)

# =============================================================================
# Selections des wavenumber LP
# =============================================================================
K_cut = 7e3
id_K_cut = np.argwhere(K1>K_cut)[0,0]
K1 = K1[:id_K_cut]
F11_K1 = F11_K1[:id_K_cut] 

dK = np.mean(np.diff(K1))
K = K1[:-2]
Klog = np.logspace(np.log10(np.min(K1)),np.log10(np.max(K1)),len(K))
#tck = splrep(K1,F11_K1)
#F11_Klog = splev(Klog,tck) 
# =============================================================================
# Filtrage de F11_K1
# =============================================================================

b,a = signal.butter(5,4.8e-3,fs=1./dK) # 0.03 marche pas mal
#F11_filt = np.exp(signal.filtfilt(b,a,np.log(F11_Klog)))
F11_filt = signal.filtfilt(b,a,F11_K1)

# =============================================================================
# Conversion F11_K1 --> E(K)
# =============================================================================
bE,aE = signal.butter(3,1.4e-3,fs=1./dK)
E = np.power(K,3.)*np.diff( (1./K1[:-1])*np.diff(F11_filt)/dK )/dK
E_filt = signal.filtfilt(bE,aE,E)

plt.figure();
plt.plot(K1,F11_K1,label='Raw F_11');
plt.plot(K1,F11_filt,label='Filt. F_11');
plt.plot(K,E_filt,label='Filt. E')
plt.legend()
plt.xscale('log');plt.yscale('log')


#tckE = splrep(K,E)
#E_log = splev(Klog,tckE)
#E_filt = np.exp(signal.filtfilt(bE,aE,np.log(E_log)))

#E = np.power(Klog[:-2],3.)*np.diff( (1./Klog[:-1])*np.diff(F11_filt)/np.diff(Klog))/np.diff(Klog[:-1])
#Efilt = np.exp(signal.filtfilt(b,a,np.log10(E))*np.log(10))
#plt.figure();plt.plot(Klog[:-1],E);plt.plot(K1,F11_K1);plt.plot(Klog,F11_filt);plt.xscale('log');plt.yscale('log')
# =============================================================================
# Calcul epsilon
# =============================================================================

#epsilon = 2*nu*np.sum( K1*K1*F11_K1)*np.mean(np.diff(K1))
epsilon = 2*nu*np.sum( (K**2)*E*dK)

mean_U2 = (1./3.)*np.sum(E)*dK  #--> A peu près ok

# =============================================================================
#  Non dimensionnalize
# =============================================================================

eta = np.power(nu**3/epsilon,0.25)

F11star = F11_K1*np.power(nu,-5./4.)*np.power(epsilon,-1./4.)
F = E_filt*np.power(nu,-5./4.)*np.power(epsilon,-1./4.)

Kstar = K*eta

#Kstarlog = np.logspace(np.log10(np.min(Kstar)),np.log10(np.max(Kstar)))
#F_Kolmogorov = 0.5*np.power(epsilon,2./3.)*np.power(Kstarlog,-5./3.)*np.power(nu,-5./4.)*np.power(epsilon,-1./4.)

plt.figure()
plt.plot(K1*eta,F11star)
plt.xlabel('Kstar')
plt.ylabel('F11 star')
plt.xscale('log')
plt.yscale('log')


plt.figure()
#plt.plot(Kstarlog,Kstarlog**(-5/3))
plt.plot(Kstar,F)
plt.xlabel('Kstar')
plt.ylabel('Fstar')
plt.xscale('log')
plt.yscale('log')
