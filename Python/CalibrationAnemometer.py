# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 10:54:55 2020

@author: Gaétan
"""

# =============================================================================
# Import des bibliothèques
# =============================================================================

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib
# =============================================================================
# Data importation
# =============================================================================

filename = '2020-03-05-Calibration.txt'
T = np.genfromtxt(filename)
U = T[:,0]
Esquare = T[:,1]



# =============================================================================
# Fitting
# =============================================================================

x0 = [1.,1.,1.] #initialisation for [A,B,n]

def cost(x):
    A = x[0]
    B = x[1]
    n = x[2]
    err = Esquare - A - B*np.power(U,n) 
    return np.mean(np.square(err))


result = scipy.optimize.minimize(cost,x0)
x = result.x
A,B,n = x
# =============================================================================
# Output
# =============================================================================

print('Fitted values with normalised rms error %.3e' % (np.sqrt(cost(x))/np.sqrt(np.mean(np.square(Esquare)))))
print('A = %.3e - B = %.3e - n = %.3e'%(A,B,n))

Ufit = np.linspace(np.min(U),np.max(U),1000)



plt.figure()
plt.plot(U,Esquare,linestyle='none',marker='s',label='Provided data')
plt.plot(Ufit,A+B*np.power(Ufit,n),label='King\'s Law calibration')
plt.xlabel('Velocity U')
plt.ylabel('Square output voltage E^2')
plt.legend()
#plt.savefig('FittingAnemometer_toLatex.pgf')