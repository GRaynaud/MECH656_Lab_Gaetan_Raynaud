# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 13:42:22 2020

@author: Gaétan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
#import tikzplotlib


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)
plt.rc('axes',titlesize=16)
plt.rc('legend',fontsize=14)
plt.rc('figure',titlesize=20)

def Custom_PDF(u,N_reg,f_LP=3e-3):
    '''
    Compute and return the Probability Density Function (PDF)
    of variable u on a N_reg equally spaced range of u values
    u : input vector of values to analyse
    N_reg : Number of points to evaluate PDF
    f_LP = 3e-3 : LowPass dimensionless frequency of Butterworth filter for PDF
    frequency min : 1. --> n variations at higer frequency thant 1./dx
    frequency max : 1./N_reg --> All variations are of wavelength > N_reg are cut
    '''
    
    # Step 1 : sort u
    index_sorted = np.argsort(u)
    u_sorted = u[index_sorted]

    #Step 2 : Construction of Cumulative Density Function x --> F
    N = len(u_sorted)
    tol = (u_sorted[-1]-u_sorted[0])*1e-5
    
    x = []
    F = []
    x.append(u_sorted[0])
    F.append(0.)
    
    for k in range(1,N):
        if u_sorted[k]>u_sorted[k-1]+tol: # Evite les redondances
            x.append(u_sorted[k]) # A simplifier par x = u_sorted
            F.append(1.*k/N)
        
    x = np.asarray(x)
    F = np.asarray(F)
    
    # Step 3 : interpolation of F on regular mesh x_reg
    x_reg = np.linspace(u_sorted[0],u_sorted[-1],N_reg)
    F_reg = np.interp(x_reg, x, F)

    
    # Step 4 : d/dx of CDF --> PDF
    dx = (x_reg[-1]-x_reg[0])/N_reg
    x_o = x_reg[:-1]
    f_o = np.diff(F_reg)/dx
    
    # Step 5 : Low pass filtering
    b,a = signal.butter(3,f_LP/dx,fs=1./dx)
    f_filtered = signal.filtfilt(b,a,f_o)
    return x_o, f_filtered


# =============================================================================
# Testing
# =============================================================================

filename_data = 'TraitementDataSpec/U_series_spec020'
U = np.genfromtxt(filename_data)
fs = 8000 #Hz
N_reg = 1000
x,f = Custom_PDF(U,N_reg,5e-2)



# =============================================================================
# Verification
# =============================================================================

IntPDF = np.sum(f)*np.mean(np.diff(x))
print('Integral PDF - 1. = %.2e' % (IntPDF-1.))


# =============================================================================
# Comparaison avec fonctions existantes
# =============================================================================
Nbins = 100
n_hist,u_hist,_ = plt.hist(U,bins=Nbins)
plt.close()
plt.plot(u_hist[:-1],n_hist*Nbins/(len(U)*(np.max(U)-np.min(U))),label='Numpy')
plt.plot(x,f,label='Custom')
plt.legend(title='Method')
plt.xlabel('Speed range $u$ (m/s)')
plt.ylabel('PDF (s/m)')
plt.savefig('Comparison_PDF_Custom_Hist.pgf')
#tikzplotlib.save("Comparison_PDF_Custom_Hist.tex")#, encoding="utf-8")

# =============================================================================
# Plot PDF @ 3 measurments points
# =============================================================================

listey = [0.20,5.50,20.]
listeystr = ['020','550','2000']

plt.figure()

for k in range(3):
    filename = 'TraitementDataSpec/U_series_spec'+listeystr[k]
    U = np.genfromtxt(filename)
    
    N_reg = 1000
    f_LP = 5e-2
    
    x,f = Custom_PDF(U,N_reg,f_LP)
    
    IntPDF = np.sum(f)*np.mean(np.diff(x))
    print('Integral PDF - 1. = %.2e' % (IntPDF-1.))
    
    plt.plot(x,f,label=str(listey[k]))

plt.legend(title='$y$ (mm)')
plt.xlabel('Speed range $u$ (m/s)')
plt.ylabel('PDF (s/m)')
plt.savefig('Dimensionnal_PDF.pgf')

#tikzplotlib.save("Dimensionnal_PDF.tex")#, encoding="utf-8")
