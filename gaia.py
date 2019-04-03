''' Gaia DR2 cross match filtering '''

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import abj2016

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

# df = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/CoRoGEE_x_gaiadr2.csv')
df = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/GAP_Gaia_x_gaiadr2.csv')
df1 = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_x_gaiadr2.csv')
df2 = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/K2_Gaia_x_gaiadr2.csv')
G3 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/GAP3')
G6 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/GAP6')
# C3s = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Normal/C3_07092018')
#df1.rename(columns={'#Id':'EPIC'},inplace=True)
# C6s = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Normal/C6_07092018')
#df2.rename(columns={'#Id':'EPIC'},inplace=True)

#df1 = pd.merge(df1,df[['EPIC','parallax','parallax_error']],how='inner',on=['EPIC'])
#df1 = df1.drop_duplicates(subset=['EPIC']).reset_index(drop=True)

#df2 = pd.merge(df2,df[['EPIC','parallax','parallax_error']],how='inner',on=['EPIC'])
#df2 = df2.drop_duplicates(subset=['EPIC']).reset_index(drop=True)

''' Convert parallax to distance using the Astraatmadja and Bailer-Jones 2016 method - results in kpc '''
# df['dist_G'] = 0.0
# df['sig_dist_G'] = 0.0
df1['dist_ABJ'] = 0.0
df1['sig_dist_ABJ'] = 0.0
df1['dist_1P'] = 0.0
df1['sig_dist_1P'] = 0.0
# df2['dist_G'] = 0.0
# df2['sig_dist_G'] = 0.0
# for i in range(len(df)):
#     df['dist_G'].iloc[i] = abj2016.distpdf(df['parallax'].iloc[i], df['parallax_error'].iloc[i]).modedist   # Initiates the distpdf object
#     df['sig_dist_G'].iloc[i] = abj2016.distpdf(df['parallax'].iloc[i], df['parallax_error'].iloc[i]).diststd   # Initiates the distpdf object

for i in range(len(df1)):
    if df1['EPIC'].iloc[i] < 207000000:
        df1['parallax'].iloc[i] = df1['parallax'].iloc[i] + 0.01506
    else:
        df1['parallax'].iloc[i] = df1['parallax'].iloc[i] + 0.00149

    df1['dist_ABJ'].iloc[i] = abj2016.distpdf(df1['parallax'].iloc[i], df1['parallax_error'].iloc[i]).modedist   # Initiates the distpdf object
    df1['sig_dist_ABJ'].iloc[i] = abj2016.distpdf(df1['parallax'].iloc[i], df1['parallax_error'].iloc[i]).diststd   # Initiates the distpdf object
    df1['dist_1P'].iloc[i] = 1/df1['parallax'].iloc[i]
    df1['sig_dist_1P'].iloc[i] = (df1['parallax_error'].iloc[i]/df1['parallax'].iloc[i])*df1['dist_1P'].iloc[i]

#for i in range(len(df2)):
#    if df1['EPIC'].iloc[i] < 207000000:
#        df1['parallax'].iloc[i] = df1['parallax'].iloc[i] - 0.006
#    else: df1['parallax'].iloc[i] = df1['parallax'].iloc[i] - 0.017
#    df2['dist_G'].iloc[i] = abj2016.distpdf(df2['parallax'].iloc[i], df2['parallax_error'].iloc[i]).modedist   # Initiates the distpdf object
#    df2['sig_dist_G'].iloc[i] = abj2016.distpdf(df2['parallax'].iloc[i], df2['parallax_error'].iloc[i]).diststd   # Initiates the distpdf object

#df1.rename(columns={'EPIC':'#Id'},inplace=True)
#df2.rename(columns={'EPIC':'#Id'},inplace=True)

# df.to_csv('/home/ben/K2_Poles/Mass_Distr_In/CoRoGEE_x_gaiadr2.csv',index=False)
df1.to_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_x_gaiadr2.csv',index=False)
#df2.to_csv('/home/ben/K2_Poles/Mass_Distr_In/Gaia/K2_Gaia_x_gaiadr2.csv',index=False)
sys.exit()

C3s['sig_dist'] = (C3s['dist_68U']-C3s['dist_68L'])/2
C6s['sig_dist'] = (C6s['dist_68U']-C6s['dist_68L'])/2

df = df.sort_values(['dist_G']).groupby(['CoRoT_ID']).head(1).reset_index(drop=True)
df1 = df1.sort_values(['dist_G']).groupby(['2MASS']).head(1).reset_index(drop=True)
df2 = df2.sort_values(['dist_G']).groupby(['2MASS']).head(1).reset_index(drop=True)

''' Merge EPIC ID to K2 seismic sample '''
df2_G3 = pd.merge(df2,G3[['EPIC','2MASS']],how='inner',on=['2MASS'])
df2_G6 = pd.merge(df2,G6[['EPIC','2MASS']],how='inner',on=['2MASS'])

C3 = pd.merge(df2_G3,C3s,how='inner',on=['EPIC'])
C3['ddist'] = C3['dist']*1e-3 - C3['dist_G']
C3 = C3.dropna(subset=['parallax'])
C6 = pd.merge(df2_G6,C6s,how='inner',on=['EPIC'])
C6['ddist'] = C6['dist']*1e-3 - C6['dist_G']
C6 = C6.dropna(subset=['parallax'])
# print(len(C3),len(C6))

# C3.to_csv('/home/ben/K2_Poles/Mass_Distr_In/C3_Gaia_Sample',index=False)
# C6.to_csv('/home/ben/K2_Poles/Mass_Distr_In/C6_Gaia_Sample',index=False)

''' Stars with sig_parallax > 10% '''
C3_pe = C3[(C3['parallax_error']/C3['parallax'])>0.1]
C6_pe = C6[(C6['parallax_error']/C6['parallax'])>0.1]

''' Stars with duplicity flag '''
C3d = C3[C3['duplicated_source'] == 1]
C6d = C6[C6['duplicated_source'] == 1]

''' Compare Distances and Residuals '''
# fig, (ax,ax1) = plt.subplots(1,2,sharey=True)
fig, ((ax,ax1),(rax,rax1)) = plt.subplots(2,2,sharex=True,sharey='row',gridspec_kw={"height_ratios" : [5,1]})
ax.errorbar(C3['dist']*1e-3,C3['dist_G'],xerr=C3['sig_dist']*1e-3,yerr=C3['sig_dist_G'],fmt='o',ecolor='k',alpha=0.5)
# ax.scatter(C3['dist']*1e-3,C3['dist_G'])
ax.scatter(C3_pe['dist']*1e-3,C3_pe['dist_G'],color='g')
# ax.scatter(C3d['dist']*1e-3,C3d['dist_G'],alpha=1.0,color='r')
ax.set_ylabel(r'Gaia Dist. [kpc]')
ax.plot([0,7.5],[0,7.5],linestyle='--',color='k')
ax.set_xlim(0,7.5)
ax.set_title(r'C3')
ax1.plot([0,7.5],[0,7.5],linestyle='--',color='k')
ax1.errorbar(C6['dist']*1e-3,C6['dist_G'],xerr=C6['sig_dist']*1e-3,yerr=C6['sig_dist_G'],fmt='o',ecolor='k',alpha=0.5)
# ax1.scatter(C6['dist']*1e-3,C6['dist_G'])
ax1.scatter(C6_pe['dist']*1e-3,C6_pe['dist_G'],color='g')
# ax1.scatter(C6d['dist']*1e-3,C6d['dist_G'],color='r',alpha=1.0)
ax1.set_xlim(0,7.5)
ax1.set_title(r'C6')

rax.scatter(C3['dist']*1e-3,100*C3['ddist']/(C3['dist']*1e-3),alpha=0.5)
rax.plot([0,7.5],[0,0])
rax.set_xlabel(r'Seismic Dist. [kpc]')
rax.set_ylabel(r'Residual [%]')
rax.set_ylim(-100,100)
rax1.scatter(C6['dist']*1e-3,100*C6['ddist']/(C6['dist']*1e-3),alpha=0.5)
rax1.plot([0,7.5],[0,0])
rax1.set_ylim(-100,100)
rax1.set_xlabel(r'Seismic Dist. [kpc]')

plt.tight_layout()
plt.show()
sys.exit()
''' Plot Uncertainties '''
fig, (ax,ax1) = plt.subplots(1,2,sharey=True)
ax.scatter(C3['dist']*1e-3,C3['sig_dist']*1e-3,alpha=0.5,label=r'$\sigma_{seismic}$')
ax.scatter(C3_pe['dist']*1e-3,C3_pe['sig_dist_G'],color='g',alpha=0.5,label=r'$\sigma_{Gaia}$')
ax.set_xlabel(r'Seismic Dist. [kpc]')
ax.set_ylabel(r'$\sigma$ Dist. [kpc]')
# ax.plot([0,7.5],[0,7.5],linestyle='--',color='k')
ax.set_title(r'C3')
ax1.scatter(C6['dist']*1e-3,C6['sig_dist']*1e-3,alpha=0.5,label=r'$\sigma_{seismic}$')
ax1.scatter(C6_pe['dist']*1e-3,C6_pe['sig_dist_G'],color='g',alpha=0.5,label=r'$\sigma_{Gaia}$')
ax1.set_xlabel(r'$\sigma$ Seismic Dist. [kpc]')
# ax1.plot([0,7.5],[0,7.5],linestyle='--',color='k')
# ax1.set_xlim(0,7.5)
ax1.set_title(r'C6')
ax1.legend()
plt.tight_layout()
plt.show()
