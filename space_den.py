''' Volume test plot for space density calculations '''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import mass_distr_functs as mdf
import scipy
from scipy import integrate

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

theta = np.deg2rad(np.sqrt(116))
t1 = theta/2.0
alpha_c3 = np.deg2rad(61.4) # degrees
alpha_c6 = np.deg2rad(50.4) # degrees
z1 = 0.5 # kpc
z2 = 1.0
z3 = 1.5
z4 = 2.0

c3 = pd.read_csv('/home/bmr135/C3_ojh')
c6 = pd.read_csv('/home/bmr135/C6_ojh')

x2 = lambda x: (z1**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - x)**-3)
x3 = lambda xa: (z2**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - xa)**-3)
x4 = lambda xb: (z3**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - xb)**-3)
x5 = lambda xc: (z4**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - xc)**-3)

vol1 = integrate.quad(x2,0,theta)
vol2 = integrate.quad(x3,0,theta)
vol3 = integrate.quad(x4,0,theta)
vol4 = integrate.quad(x5,0,theta)

vola = vol2[0] - vol1[0]
volb = vol3[0] - vol2[0]
volc = vol4[0] - vol3[0]

print(vol1, vol2, vol3, vola, volb)
print(volc-volb,volb-vola, vola-vol1[0])

print((z1*(2*z1*np.tan(t1))**2)/3)
print((z2*(2*z2*np.tan(t1))**2)/3)
print((z3*(2*z3*np.tan(t1))**2)/3)
