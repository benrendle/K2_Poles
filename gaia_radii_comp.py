''' Quick-fire comparisons between Seismic and Gaia Radii '''

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sys

df = pd.read_csv('/home/bmr135/Desktop/C3_gaia')
df1 = pd.read_csv('/home/bmr135/Desktop/C6_gaia')
df2 = pd.read_csv('/home/bmr135/Desktop/C3_spec_gaia')
df3 = pd.read_csv('/home/bmr135/Desktop/C6_spec_gaia')

f1, ax = plt.subplots()
ax.plot([3,25],[3,25],color='orange',alpha=0.5,linestyle='--')
ax.scatter(df2['radius_val'],df2['radius'],label=r'Phot.')
ax.scatter(df2['radius_val'],df2['sRad'],label=r'Spec.')
x = np.linspace(3,25,100)
ax.plot(x,x*1.125,color='m',alpha=0.5,linestyle='--')
ax.set_xlabel(r'R$_{\rm{Gaia}}$')
ax.set_ylabel(r'R$_{\rm{Other}}$')
ax.set_xlim(3,25)
ax.set_ylim(3,25)
ax.legend()

f1, ax = plt.subplots()
ax.plot([3,25],[3,25],color='orange',alpha=0.5,linestyle='--')
ax.scatter(df2['radius'],df2['sRad'],label=r'Seismic')
ax.scatter(df2['radius'],df2['radius_val'],label=r'Gaia')
x = np.linspace(3,25,100)
ax.plot(x,x*0.875,color='m',alpha=0.5,linestyle='--')
ax.set_xlabel(r'R$_{\rm{Phot.}}$')
ax.set_ylabel(r'R$_{\rm{Gaia}}$')
ax.set_xlim(3,25)
ax.set_ylim(3,25)
ax.legend()
plt.show()
