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
#import tikzplotlib

# =============================================================================
# Acquisition des données
# =============================================================================

filename0020 = 'TraitementDataSpec/specrsltfilnam020.11mar2020.txt'
filename0550 = 'TraitementDataSpec/specrsltfilnam550.11mar2020.txt'
filename2000 = 'TraitementDataSpec/specrsltfilnam2000.11mar2020.txt'

T0020 = np.genfromtxt(filename0020,skip_header=1)
T0550 = np.genfromtxt(filename0550,skip_header=1)
T2000 = np.genfromtxt(filename2000,skip_header=1)

f0020 = T0020[:,0]
Spec0020 = T0020[:,1]

f0550 = T0550[:,0]
Spec0550 = T0550[:,1]

f2000 = T2000[:,0]
Spec2000 = T2000[:,1]

# =============================================================================
#  Conversion temporel --> spatial
# y = [  0.2 ,   5.5 ,  20.  ]
# mean_U = [2.3158,  4.1619, 5.0156]
# =============================================================================
nu = 1.5e-5 #m/s à vérifier
U0020 = 2.3158 #m/s
U0550 = 4.1619
U2000 = 5.0156


K1_0020 = 2*np.pi*f0020/U0020
F11_K1_0020 = Spec0020*U0020/(2*np.pi)

K1_0550 = 2*np.pi*f0550/U0550
F11_K1_0550 = Spec0550*U0550/(2*np.pi)

K1_2000 = 2*np.pi*f2000/U2000
F11_K1_2000 = Spec2000*U2000/(2*np.pi)

# =============================================================================
# Selections des wavenumber LP
# =============================================================================
K_cut_0020 = 7e3
id_K_cut_0020 = np.argwhere(K1_0020>K_cut_0020)[0,0]
K1_0020 = K1_0020[:id_K_cut_0020]
F11_K1_0020 = F11_K1_0020[:id_K_cut_0020] 

K_cut_0550 = 5e3
id_K_cut_0550 = np.argwhere(K1_0550>K_cut_0550)[0,0]
K1_0550 = K1_0550[:id_K_cut_0550]
F11_K1_0550 = F11_K1_0550[:id_K_cut_0550] 

K_cut_2000 = 3.5e3
id_K_cut_2000 = np.argwhere(K1_2000>K_cut_2000)[0,0]
K1_2000 = K1_2000[:id_K_cut_2000]
F11_K1_2000 = F11_K1_2000[:id_K_cut_2000] 


dK0020 = np.mean(np.diff(K1_0020))
K0020 = K1_0020[:-2]

dK0550 = np.mean(np.diff(K1_0550))
K0550 = K1_0550[:-2]

dK2000 = np.mean(np.diff(K1_2000))
K2000 = K1_2000[:-2]

#Klog = np.logspace(np.log10(np.min(K1)),np.log10(np.max(K1)),len(K))
#tck = splrep(K1,F11_K1)
#F11_Klog = splev(Klog,tck) 
# =============================================================================
# Filtrage de F11_K1
# =============================================================================

b0020,a0020 = signal.butter(5,4.8e-3,fs=1./dK0020) # 0.03 marche pas mal
F11_filt_0020 = signal.filtfilt(b0020,a0020,F11_K1_0020)
#F11_filt = np.exp(signal.filtfilt(b,a,np.log(F11_Klog)))

b0550,a0550 = signal.butter(5,4.8e-3,fs=1./dK0550) # 0.03 marche pas mal
F11_filt_0550 = signal.filtfilt(b0550,a0550,F11_K1_0550)

b2000,a2000 = signal.butter(5,4.8e-3,fs=1./dK2000) # 0.03 marche pas mal
F11_filt_2000 = signal.filtfilt(b2000,a2000,F11_K1_2000)



# =============================================================================
# Conversion F11_K1 --> E(K)
# =============================================================================
bE0020,aE0020 = signal.butter(3,1.4e-3,fs=1./dK0020)
E0020 = np.power(K0020,3.)*np.diff( (1./K1_0020[:-1])*np.diff(F11_filt_0020)/dK0020 )/dK0020
E_filt0020 = signal.filtfilt(bE0020,aE0020,E0020)

bE0550,aE0550 = signal.butter(3,1.4e-3,fs=1./dK0550)
E0550 = np.power(K0550,3.)*np.diff( (1./K1_0550[:-1])*np.diff(F11_filt_0550)/dK0550 )/dK0550
E_filt0550 = signal.filtfilt(bE0550,aE0550,E0550)

bE2000,aE2000 = signal.butter(3,1.4e-3,fs=1./dK2000)
E2000 = np.power(K2000,3.)*np.diff( (1./K1_2000[:-1])*np.diff(F11_filt_2000)/dK2000 )/dK2000
E_filt2000 = signal.filtfilt(bE2000,aE2000,E2000)




plt.figure();
plt.plot(K1_0020,F11_K1_0020,label='Raw $F_{11}$');
plt.plot(K1_0020,F11_filt_0020,label='Filt. $F_{11}$');
plt.plot(K0020,E_filt0020,label='Filt. $E$')
plt.xlabel('$K_{11}$ - log scale ($m^{-1}$)')
plt.legend()
plt.xscale('log');plt.yscale('log')
plt.tight_layout()
plt.savefig("Filetering_protocol.pgf")

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
epsilon0020 = 2*nu*np.sum( (K0020**2)*E0020*dK0020)
mean_U2_0020 = (1./3.)*np.sum(E0020)*dK0020  #--> A peu près ok

epsilon0550 = 2*nu*np.sum( (K0550**2)*E0550*dK0550)
mean_U2_0550 = (1./3.)*np.sum(E0550)*dK0550

epsilon2000 = 2*nu*np.sum( (K2000**2)*E2000*dK2000)
mean_U2_2000 = (1./3.)*np.sum(E2000)*dK2000

# =============================================================================
#  Non dimensionnalize
# =============================================================================

eta0020 = np.power(nu**3/epsilon0020,0.25)
F11star0020 = F11_K1_0020*np.power(nu,-5./4.)*np.power(epsilon0020,-1./4.)
F0020 = E_filt0020*np.power(nu,-5./4.)*np.power(epsilon0020,-1./4.)
Kstar0020 = K0020*eta0020

eta0550 = np.power(nu**3/epsilon0550,0.25)
F11star0550 = F11_K1_0550*np.power(nu,-5./4.)*np.power(epsilon0550,-1./4.)
F0550 = E_filt0550*np.power(nu,-5./4.)*np.power(epsilon0550,-1./4.)
Kstar0550 = K0550*eta0550

eta2000 = np.power(nu**3/epsilon2000,0.25)
F11star2000 = F11_K1_2000*np.power(nu,-5./4.)*np.power(epsilon2000,-1./4.)
F2000 = E_filt2000*np.power(nu,-5./4.)*np.power(epsilon2000,-1./4.)
Kstar2000 = K2000*eta2000

#Kstarlog = np.logspace(np.log10(np.min(Kstar)),np.log10(np.max(Kstar)))
#F_Kolmogorov = 0.5*np.power(epsilon,2./3.)*np.power(Kstarlog,-5./3.)*np.power(nu,-5./4.)*np.power(epsilon,-1./4.)

plt.figure()
plt.plot(K1_0020*eta0020,F11star0020, label='0.20')
plt.plot(K1_0550*eta0550,F11star0550, label='5.50')
plt.plot(K1_2000*eta2000,F11star2000, label='20.0')
plt.xlabel('$K^{*}$')
plt.ylabel('$F_{11}^{*}$')
plt.xscale('log')
plt.yscale('log')
plt.legend(title='$y$ (mm)')
plt.tight_layout()
plt.savefig('F_11_comparison.pgf')

#plt.figure()
#plt.plot(K1_0020*eta0020,F11_filt_0020*np.power(nu,-5./4.)*np.power(epsilon0020,-1./4.), label='0.20')
#plt.plot(K1_0550*eta0550,F11_filt_0550*np.power(nu,-5./4.)*np.power(epsilon0550,-1./4.), label='5.50')
#plt.plot(K1_2000*eta2000,F11_filt_2000*np.power(nu,-5./4.)*np.power(epsilon2000,-1./4.), label='20.00')
#plt.xlabel('Kstar')
#plt.ylabel('F11 star')
#plt.xscale('log')
#plt.yscale('log')
#plt.legend(title='$y$ (mm)')
#plt.tight_layout()

t=5.

plt.figure()
plt.plot(Kstar0020,F0020,label='0.20')
plt.plot(Kstar0550,F0550,label='5.50')
plt.plot(Kstar2000,F2000,label='20.00')
plt.plot([4e-2,4e-2*t,4e-2,4e-2],[18.9,18.9*t**(-5./3.),18.9*t**(-5./3.),18.9 ], color='black')
plt.text(4e-2,0.5, '$-5/3$ slope')
plt.xlabel('$K^{*}$')
plt.ylabel('$F^{*}$')
plt.xscale('log')
plt.yscale('log')
plt.legend(title='$y$ (mm)')
plt.tight_layout()
plt.savefig('F_comparison.pgf')