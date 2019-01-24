# Comparison of Yvonne's, Savita's and Benoit's C3/C6 numax and dnu values
# Plots styled with the help of GRD
# Last Modified: 27/04/17, Ben Rendle

import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.mlab as mlab
import pandas as pd
import sys
from pandas import DataFrame, read_csv
import numpy as np
from scipy import stats
from scipy.stats import norm
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde
import scipy.odr.odrpack as odrpack
from sklearn.metrics import r2_score
# from pyqt_fit import kde
import colormaps
import K2_plot as k2p
import K2_properties as prop
import K2_constants as const
import K2_data as dat
from numbers import Number
import seaborn as sns
import time

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
plt.rcParams["font.family"] = "serif"

''' Dropbox Path '''
ext_DB = '/home/bmr135/' # Work
# ext_DB = '/home/ben/'   # Laptop
''' GA directory '''
ext_GA = '/media/bmr135/SAMSUNG/' # Hard-Drive

def hist_orig(df,df1,df2,bins,n):
    '''
    df = original data
    df1 = reduced data set
    '''
    fig, axes = plt.subplots(3,2)
    ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
    hist = [0,0,0,0,0]
    hist0 = [0,0,0,0,0]
    b = [0,0,0,0,0]
    d = [0,0,0,0,0]
    ratio = [0,0,0,0,0]

    hist0[0], d[0], c = ax0.hist(df['mass'],bins=bins[0],histtype='step',linewidth=2,label=r'Full')
    hist[0], b[0], c = ax0.hist(df1['mass'],bins=bins[0],histtype='step',linewidth=2,label=r'EPIC Rad.')
    hist[0], b[0], c = ax0.hist(df2['mass'],bins=bins[0],histtype='step',linewidth=2,label=r'Gaia Rad.')
    ax0.set_xlabel(r'Mass [M$_{\odot}$]')
    ax0.legend()
    if n == 6:
        ax0.set_xlim(0.5,2.5)

    hist0[1], d[1], c = ax1.hist(df['Radius'],bins=bins[1],histtype='step',linewidth=2)
    hist[1], b[1], c = ax1.hist(df1['Radius'],bins=bins[1],histtype='step',linewidth=2)
    hist[1], b[1], c = ax1.hist(df2['Rgaia'],bins=bins[1],histtype='step',linewidth=2)
    ax1.set_xlabel(r'Radius [R$_{\odot}$]')
    if n == 6:
        ax1.set_xlim(0,30)

    hist0[2], d[2], c = ax2.hist(df['Teff'],bins=bins[2],histtype='step',linewidth=2)
    hist[2], b[2], c = ax2.hist(df1['Teff'],bins=bins[2],histtype='step',linewidth=2)
    hist[2], b[2], c = ax2.hist(df2['Teff'],bins=bins[2],histtype='step',linewidth=2)
    ax2.set_xlabel(r'T$_{\rm{eff}}$ [K]')
    if n == 6:
        ax2.set_xlim(4000,5250)

    hist0[3], d[3], c = ax3.hist(df['[Fe/H]'],bins=bins[3],histtype='step',linewidth=2)
    hist[3], b[3], c = ax3.hist(df1['[Fe/H]'],bins=bins[3],histtype='step',linewidth=2)
    hist[3], b[3], c = ax3.hist(df2['[Fe/H]'],bins=bins[3],histtype='step',linewidth=2)
    ax3.set_xlabel(r'[Fe/H] [dex]')
    if n == 6:
        ax3.set_xlim(-1.2,0.2)

    hist0[4], d[4], c = ax4.hist(df['logg'],bins=bins[4],histtype='step',linewidth=2)
    hist[4], b[4], c = ax4.hist(df1['logg'],bins=bins[4],histtype='step',linewidth=2)
    hist[4], b[4], c = ax4.hist(df2['glogg'],bins=bins[4],histtype='step',linewidth=2)
    ax4.set_xlabel(r'log$_{10}$(g)')
    if n == 6:
        ax4.set_xlim(1,4)

    ax5.axis('off')
    ratio[0] = (hist[0]/hist0[0])*100
    ratio[1] = (hist[1]/hist0[1])*100
    ratio[2] = (hist[2]/hist0[2])*100
    ratio[3] = (hist[3]/hist0[3])*100
    ratio[4] = (hist[4]/hist0[4])*100
    # ax0.set_title(r'Cut Applied: '+cut)
    plt.tight_layout()
    fig.savefig('Det_prob_cut_props.pdf', bbox_inches='tight')

''' Read in simulated and real data '''
# besa3, besa6 = dat.BESANCON()
# print("Besancon")
TRILEGAL_C3, TRILEGAL_C6 = dat.TRILEGAL()
C3 = dat.C3_cat()
C6 = dat.C6_cat()
GAP3, GAP6 = dat.K2_GAP()
GAP3 = GAP3.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
GAP6 = GAP6.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
Yvonne_C3, Yvonne_C6, Yvonne_EC6, Yvonne_EC3 = dat.Yvonne()
Savita_C3, Savita_C6, Savita_EC3, Savita_EC6 = dat.Savita()
Benoit_C3, Benoit_C6, Everest_C3, Everest_C6 = dat.Benoit()
RAVE3, RAVE6 = dat.RAVE()
GES3 = dat.Gaia_ESO()
APO3, APO6 = dat.APOGEE()
LAMOST3, LAMOST6 = dat.LAMOST()
oc = dat.occurrence()
print( "All files in")

''' Add K2P2 occurrence to GAP Target Lists '''
GAP6 = dat.n_epics(GAP6,oc)

''' Preparing GAP for probability detections '''
''' Merge data with GAP target lists '''
YC3 = pd.merge(Yvonne_C3,GAP3,how='inner',on=['EPIC'])
SC3 = pd.merge(Savita_C3,GAP3,how='inner',on=['EPIC'])
BC3 = pd.merge(Benoit_C3,GAP3,how='inner',on=['EPIC'])

''' Complete asteroseismic lists '''
camp3_0 = pd.concat([YC3,SC3,BC3],ignore_index=True)
camp3_0 = camp3_0.drop_duplicates(subset=['EPIC'])
camp3_0 = camp3_0.reset_index(drop=True)
camp3_0 = camp3_0.fillna(value='NaN',method=None)


GAP3 = pd.merge(GAP3,camp3_0[['EPIC','Bnumax','nmx','Snumax','BDnu','dnu','SDnu']],how='inner',on=['EPIC'])
GAP3 = prop.single_seismo(GAP3,['Bnumax','nmx','Snumax'],'NUMAX')
GAP3 = prop.single_seismo(GAP3,['BDnu','dnu','SDnu'],'DNU')
GAP3['NUMAX']=pd.to_numeric(GAP3['NUMAX'])
GAP3 = GAP3[GAP3['NUMAX'] < 280]
GAP3['foma'] = GAP3['Radius']**-1.85 * (GAP3['Teff']/5777.0)**0.92 * 3090 # numax for C3 GAP (frequency of maximum amplitude) <- coefficients from Bill
GAP3['Lumo'] = GAP3['Radius']**2 * (GAP3['Teff']/const.solar_Teff)**4
GAP3['fomag'] = GAP3['Rgaia']**-2 * (GAP3['Teff']/5777.0)**0.92 * 3090 # numax for C3 GAP (frequency of maximum amplitude)
GAP3['Lumog'] = GAP3['Rgaia']**2 * (GAP3['Teff']/const.solar_Teff)**4
GAP3_v2 = GAP3[GAP3['foma'] < 280]
GAP3_v2 = GAP3_v2[GAP3_v2['foma'] > 10]
GAP3_v2 = GAP3_v2[GAP3_v2['imag'] > 0.0]
GAP3_v2 = GAP3_v2.reset_index(drop=True)
GAP3_v2 = prop.det_prob_GAP(GAP3_v2,'NUMAX',3090,135.1)

GAP3_v3 = GAP3[GAP3['fomag'] < 280]
GAP3_v3 = GAP3_v3[GAP3_v3['fomag'] > 10]
GAP3_v3 = GAP3_v3[GAP3_v3['imag'] > 0.0]
GAP3_v3 = GAP3_v3.reset_index(drop=True)
GAP3_v3 = prop.det_prob_GAP_gaia(GAP3_v3,'NUMAX',3090,135.1)


# GAP3_v2.to_csv(ext_GA+'GA/K2Poles/GAP3_det_prob_gaia',index=False)

GAP6['foma'] = GAP6['Radius']**-1.85 * (GAP6['Teff']/5777.0)**0.92 * 3090 # numax for C6 GAP (frequency of maximum amplitude)
GAP6['Lumo'] = GAP6['Radius']**2 * (GAP6['Teff']/const.solar_Teff)**4
GAP6['fomag'] = GAP6['Rgaia']**-1.85 * (GAP6['Teff']/5777.0)**0.92 * 3090 # numax for C6 GAP (frequency of maximum amplitude)
GAP6['Lumog'] = GAP6['Rgaia']**2 * (GAP6['Teff']/const.solar_Teff)**4
GAP6_v2 = GAP6[GAP6['foma'] < 280]
GAP6_v2 = GAP6_v2[GAP6_v2['foma'] > 10]
GAP6_v2 = GAP6_v2[GAP6_v2['imag'] > 0.0]
GAP6_v2 = GAP6_v2.reset_index(drop=True)
GAP6_v2 = prop.det_prob_GAP(GAP6_v2,'foma',3090,135.1)
GAP6_v3 = GAP6[GAP6['fomag'] < 280]
GAP6_v3 = GAP6_v3[GAP6_v3['fomag'] > 10]
GAP6_v3 = GAP6_v3[GAP6_v3['imag'] > 0.0]
GAP6_v3 = GAP6_v3.reset_index(drop=True)
GAP6_v3 = prop.det_prob_GAP_gaia(GAP6_v3,'fomag',3090,135.1)

# GAP3_v2 = GAP3_v2[GAP3_v2['prob_s'] >= 0.95]
# GAP3_v3 = GAP3_v3[GAP3_v3['prob_s_gaia'] >= 0.95]
# GAP6_v2 = GAP6_v2[GAP6_v2['prob_s'] >= 0.95]
# GAP6_v3 = GAP6_v3[GAP6_v3['prob_s_gaia'] >= 0.95]

# GAP6_v2.to_csv(ext_GA+'GA/K2Poles/GAP6_det_prob_gaia',index=False)
K2_camp = pd.concat([GAP3,GAP6],ignore_index=True)
K2_camp = K2_camp.reset_index(drop=True)
K2_camp_v2 = pd.concat([GAP3_v2,GAP6_v2],ignore_index=True)
K2_camp_v2 = K2_camp_v2.reset_index(drop=True)
K2_camp_v3 = pd.concat([GAP3_v3,GAP6_v3],ignore_index=True)
K2_camp_v3 = K2_camp_v3.reset_index(drop=True)
K2_camp_v3[['EPIC','Rgaia','radius_val','Kabs','glogg']].to_csv(ext_DB+'K2_Poles/Mass_Distr_In/K2_det_prob_gaia',index=False)




# fig,ax = plt.subplots()
# a = np.where(K2_camp_v2['prob_s'] >= 0.95)
# b = np.where(K2_camp_v3['prob_s_gaia'] >= 0.95)
# ax.scatter(K2_camp['JK'],K2_camp['Kabs'],alpha=0.5,label=r'GAP')
# ax.scatter(K2_camp_v2['JK'].iloc[a],K2_camp_v2['Kabs'].iloc[a],alpha=0.5,label=r'R$_{\rm{catalogue}}$')
# ax.scatter(K2_camp_v3['JK'].iloc[b],K2_camp_v3['Kabs'].iloc[b],alpha=0.5,label=r'R$_{\rm{Gaia}}$')
# ax.set_xlabel(r'J - K',fontsize=15)
# ax.set_ylabel(r'K$_{\rm{abs}}$',fontsize=15)
# ax.set_xlim(0.475,1.325)
# ax.invert_yaxis()
# ax.legend(loc=4)
# # plt.show()
# fig.savefig('Det_prob_cut_HRD.pdf', bbox_inches='tight')
# sys.exit()

# bins = [np.linspace(min(K2_camp['mass']),max(K2_camp['mass']),50), \
#         np.linspace(min(K2_camp['Radius']),25,50), \
#         np.linspace(min(K2_camp['Teff']),max(K2_camp['Teff']),50), \
#         np.linspace(min(K2_camp['[Fe/H]']),max(K2_camp['[Fe/H]']),50), \
#         np.linspace(min(K2_camp['logg']),max(K2_camp['logg']),50)]
# hist_orig(K2_camp,K2_camp_v2,K2_camp_v3,bins,0)

# cols = ['EPIC','2MASS','RA','Dec']#,'Teff','[Fe/H]','logg']
# GAP3_v2.to_csv('/home/ben/Desktop/C3_GAP_Gaia',columns=cols,index=False)
# GAP6_v2.to_csv('/home/ben/Desktop/C6_GAP_Gaia',columns=cols,index=False)
# sys.exit()

# GAP3 = GAP3_v3
# GAP6 = GAP6_v3

# ''' Merge data with GAP target lists '''
# YC3 = pd.merge(Yvonne_C3,GAP3,how='inner',on=['EPIC'])
# YC6 = pd.merge(Yvonne_C6,GAP6,how='inner',on=['EPIC'])
# SC3 = pd.merge(Savita_C3,GAP3,how='inner',on=['EPIC'])
# SC6 = pd.merge(Savita_C6,GAP6,how='inner',on=['EPIC'])
# BC3 = pd.merge(Benoit_C3,GAP3,how='inner',on=['EPIC'])
# BC6 = pd.merge(Benoit_C6,GAP6,how='inner',on=['EPIC'])
# EC6 = pd.merge(Everest_C6,GAP6,how='inner',on=['EPIC'])
# YEC3 = pd.merge(Yvonne_EC3,GAP3,how='inner',on=['EPIC'])
# YEC6 = pd.merge(Yvonne_EC6,GAP6,how='inner',on=['EPIC'])
# SEC3 = pd.merge(Savita_EC3,GAP3,how='inner',on=['EPIC'])
# SEC6 = pd.merge(Savita_EC6,GAP6,how='inner',on=['EPIC'])
# EC3 = pd.merge(Everest_C3,GAP3,how='inner',on=['EPIC'])
# GG3 = pd.merge(GES3,GAP3,how='inner',on=['EPIC'])
# AC3 = pd.merge(APO3,GAP3,how='inner',on=['EPIC'])
# AC6 = pd.merge(APO6,GAP6,how='inner',on=['EPIC'])
# LC3 = pd.merge(LAMOST3,GAP3,how='inner',on=['EPIC'])
# LC6 = pd.merge(LAMOST6,GAP6,how='inner',on=['EPIC'])


# ''' Complete asteroseismic lists '''
# camp3_0 = pd.concat([YC3,SC3,BC3],ignore_index=True)
# camp3_0 = camp3_0.drop_duplicates(subset=['EPIC'])
# camp3_0 = camp3_0.reset_index(drop=True)
# camp3_0 = camp3_0.fillna(value='NaN',method=None)
#
# camp6_0 = pd.concat([YC6,SC6,BC6],ignore_index=True)
# camp6_0 = camp6_0.drop_duplicates(subset=['EPIC'])
# camp6_0 = camp6_0.reset_index(drop=True)
# camp6_0 = camp6_0.fillna(value='NaN',method=None)
#
# C3R = pd.merge(camp3_0[['EPIC']],GAP3_v2,how='inner',on=['EPIC'])
# C3Rg = pd.merge(camp3_0[['EPIC']],GAP3_v3,how='inner',on=['EPIC'])
# C6R = pd.merge(camp6_0[['EPIC']],GAP6_v2,how='inner',on=['EPIC'])
# C6Rg = pd.merge(camp6_0[['EPIC']],GAP6_v3,how='inner',on=['EPIC'])
#
# print(len(GAP3_v3),len(C3Rg))
# print(len(GAP3_v2),len(C3R))
# print(len(GAP6_v3),len(C6Rg))
# print(len(GAP6_v2),len(C6R))
#
# C3Rg = pd.merge(C3Rg,camp3_0[['EPIC','Bnumax','nmx','Snumax','BDnu','dnu','SDnu']],how='inner',on=['EPIC'])
# C3Rg = prop.single_seismo(C3Rg,['Bnumax','nmx','Snumax'],'NUMAX')
# C3Rg = prop.single_seismo(C3Rg,['BDnu','dnu','SDnu'],'DNU')
# C3Rg = C3Rg[C3Rg['NUMAX'] < 280]
# C3Rg['sRad'] = (C3Rg['NUMAX']/3090) * (C3Rg['DNU']/135.1)**-2  * (C3Rg['Teff']/5777)**0.5
cthree = pd.merge(GAP3_v2,GAP3_v3[['EPIC','Rgaia']],how='inner',on=['EPIC'])
# x = np.linspace(3,23,100)
# plt.figure()
# plt.plot(x,x,'r',alpha=0.5,linestyle='--')
# plt.scatter(cthree['Radius'],cthree['Rgaia_y'])
# plt.xlabel(r'R(EPIC)')
# plt.ylabel(r'R(Gaia)')
# plt.xlim(3,23)
# plt.ylim(3,23)

# nx, (ax,ax1) = plt.subplots(1,2)
# y = np.linspace(10,280,271)
# ax.scatter(C3Rg['NUMAX'],C3Rg['fomag'])
# ax.plot(y,y,'r',alpha=0.5,linestyle='--')
# ax.set_xlabel(r'$\nu_{\rm{max},true}$')
# ax.set_ylabel(r'$\nu_{\rm{max},scaling}$')
# ax.set_xlim(10,280)
# ax.set_ylim(10,280)
# ax1.hist(GAP3_v3['fomag'],bins=np.linspace(10,280,50),alpha=0.5,label=r'GAP$_{Gaia}$, predicted')
# ax1.hist(C3Rg['fomag'],bins=np.linspace(10,280,50),alpha=0.5,label=r'GAP$_{Gaia}$, actual')
# ax1.set_xlabel(r'$\nu_{\rm{max},scaling}$, C3')
# ax1.legend()
# ax1.set_xlim(10,280)

''' Detection Probability Plots '''
prob, ax = plt.subplots(1)
x = np.linspace(0,1,21)
# print(len(C3R))
a = len(GAP3_v2)
b = len(GAP3_v3)

for i in x:
    GAP3v2 = GAP3_v2[GAP3_v2['prob_s'] >= i]
    GAP3v3 = GAP3_v3[GAP3_v3['prob_s_gaia'] >= i]
    # GAP6v2 = GAP6_v2[GAP6_v2['prob_s'] >= i]
    # GAP6v3 = GAP6_v3[GAP6_v3['prob_s_gaia'] >= i]
    C3R = pd.merge(camp3_0[['EPIC']],GAP3v2,how='inner',on=['EPIC'])
    C3Rg = pd.merge(camp3_0[['EPIC']],GAP3v3,how='inner',on=['EPIC'])
    # C6R = pd.merge(camp6_0[['EPIC']],GAP6v2,how='inner',on=['EPIC'])
    # C6Rg = pd.merge(camp6_0[['EPIC']],GAP6v3,how='inner',on=['EPIC'])
    ax.scatter(i,len(GAP3v3),color='b')
    ax.scatter(i,len(C3Rg),color='b',marker='D')
    ax.scatter(i,len(GAP3v2),color='orange')
    ax.scatter(i,len(C3R),color='orange',marker='D')
    # ax.scatter(i,len(GAP6v3),color='r')
    # ax.scatter(i,len(C6Rg),color='r',marker='D')
    # ax.scatter(i,len(GAP6v2),color='m')
    # ax.scatter(i,len(C6R),color='m',marker='D')
ax.axhline(y=a, color='orange', linestyle='--',alpha=0.5)
ax.axhline(y=b, color='blue', linestyle='--',alpha=0.5)
ax.set_xlabel(r'Detection Probability Threshold',fontsize=15)
ax.set_ylabel(r'Number of stars',fontsize=15)
# ax.legend(labels=[r'GAP$_{Gaia}$ C3',r'Actual$_{Gaia}$ C3' \
#                 ,r'GAP$_{Gaia}$ C6',r'Actual$_{Gaia}$ C6'],ncol=2)#,r'GAP$_{EPIC}$ C6',r'Actual$_{EPIC}$ C6'],ncol=2)
        # ,r'GAP$_{EPIC}$ C3',r'Actual$_{EPIC}$ C3'
plt.show()
sys.exit()
#
# fig,((ax,ax1),(ax4,ax5),(ax2,ax3),(ax6,ax7)) = plt.subplots(4,2,figsize=(8,10))
# ax.hist(GAP3_v2['Radius'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP, predicted')
# ax.hist(C3R['Radius'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP, actual')
# ax.set_xlabel(r'R [R$_{\odot}$]',fontsize=15)
# ax.set_ylabel(r'C3',fontsize=15)
# ax.set_xlim(2.5,19.)
# ax.set_title(r'With R$_{\rm{EPIC}}$',fontsize=15)
# ax.legend()
#
# ax4.hist(GAP3_v2['Hmag'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP C3, EPIC')
# ax4.hist(C3R['Hmag'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP C3, actual')
# ax4.set_xlabel(r'H',fontsize=15)
# ax4.set_ylabel(r'C3',fontsize=15)
# ax4.set_xlim(7.,12.)
#
#
# ax2.hist(GAP6_v2['Radius'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP C6, EPIC')
# ax2.hist(C6R['Radius'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP C6, actual')
# ax2.set_xlabel(r'R [R$_{\odot}$]',fontsize=15)
# ax2.set_ylabel(r'C6',fontsize=15)
# ax2.set_xlim(2.5,19.)
# ax2.set_ylabel(r'C6')
#
# ax6.hist(GAP6_v2['Vcut'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP, predicted')
# ax6.hist(C6R['Vcut'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP, actual')
# ax6.set_xlabel(r'V$_{cut}$',fontsize=15)
# ax6.set_ylabel(r'C6',fontsize=15)
# ax6.set_xlim(9.,15.)
#
#
# ax1.hist(GAP3_v3['Radius'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP, predicted',color='r')
# ax1.hist(C3Rg['Radius'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP, actual',color='k')
# ax1.set_xlabel(r'R [R$_{\odot}$]',fontsize=15)
# ax1.set_xlim(2.5,19.)
# ax1.set_title(r'With R$_{\rm{Gaia}}$',fontsize=15)
# ax1.legend()
#
# ax5.hist(GAP3_v3['Hmag'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP C3, Gaia',color='r')
# ax5.hist(C3Rg['Hmag'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP C3, actual',color='k')
# ax5.set_xlabel(r'H',fontsize=15)
# ax5.set_xlim(7.,12.)
#
# ax3.hist(GAP6_v3['Radius'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP C6, Gaia',color='r')
# ax3.hist(C6Rg['Radius'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP C6, actual',color='k')
# ax3.set_xlabel(r'R [R$_{\odot}$]',fontsize=15)
# ax3.set_xlim(2.5,19.)
#
# ax7.hist(GAP6_v3['Vcut'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP, predicted',color='r')
# ax7.hist(C6Rg['Vcut'],bins=np.linspace(0,20,50),alpha=0.5,label=r'GAP, actual',color='k')
# ax7.set_xlabel(r'V$_{cut}$',fontsize=15)
# ax7.set_xlim(9.,15.)
# plt.tight_layout()

# plt.show()
# fig.savefig('det_bias.pdf',bbox_inches='tight')
# sys.exit()

''' Data processing and parameter calculation '''
seismo_list = [YC3,YC6,SC3,SC6,BC3,BC6,EC6,YEC6,EC3,YEC3,SEC3,SEC6]
seismo_name = ['YC3','YC6','SC3','SC6','BC3','BC6','EC6','YEC6','EC3','YEC3','SEC3','SEC6']
seismo3_list = [YC3,SC3,BC3,EC3]
seismo3_name = ['YC3','SC3','BC3','EC3']
seismo6_list = [YC6,SC6,BC6,EC6]
seismo6_name = ['YC6','SC6','BC6','EC6']
numax = ['nmx','nmx','Snumax','Snumax','Bnumax','Bnumax','Enumax','nmx','Enumax','nmx','Snumax','Snumax']
Numax = [const.solar_Numax_y,const.solar_Numax_y,const.solar_Numax_s,const.solar_Numax_s, \
        const.solar_Numax,const.solar_Numax,const.solar_Numax,const.solar_Numax_y, \
        const.solar_Numax,const.solar_Numax_y,const.solar_Numax_s,const.solar_Numax_s]
dnu = ['dnu','dnu','SDnu','SDnu','BDnu','BDnu','EDnu','dnu','EDnu','dnu','SDnu','SDnu']
Dnu = [const.solar_Dnu_y,const.solar_Dnu_y,const.solar_Dnu_s,const.solar_Dnu_s, \
      const.solar_Dnu,const.solar_Dnu,const.solar_Dnu,const.solar_Dnu_y, \
      const.solar_Dnu,const.solar_Dnu_y,const.solar_Dnu_s,const.solar_Dnu_s]
sel_numax = ['nmx','Snumax','Bnumax','Enumax','nmx','Snumax','nmx','Snumax','Bnumax','Enumax','nmx','Snumax']
sel_list = [YC3,    SC3,      BC3,      EC3,   YEC3,  SEC3,   YC6,   SC6,      BC6,     EC6,   YEC6,  SEC6]

YC3,YC6,SC3,SC6,BC3,BC6,EC6,YEC6,EC3,YEC3,SEC3,SEC6 = prop.galactic_coords(seismo_list)
C3,C6,GAP3,GAP6 = prop.galactic_coords([C3,C6,GAP3,GAP6])
YC3,YC6,SC3,SC6,BC3,BC6,EC6,YEC6,EC3,YEC3,SEC3,SEC6 = prop.lmrl_comps(seismo_list,numax,dnu,Numax,Dnu,1)
YC3,SC3,BC3,EC3,YEC3,SEC3,YC6,SC6,BC6,EC6,YEC6,SEC6 = prop.selection_function(sel_list,sel_numax)
print('GAP selection funciton implemented')

''' Add detection flags to data/save out values for comparisons '''
BC3,YC3,SC3 = prop.individ(BC3,YC3,SC3)
BC6,YC6,SC6 = prop.individ(BC6,YC6,SC6)
# YEC3,EC3,SEC3 = prop.individ(YEC3,EC3,SEC3)
# YEC6,EC6,SEC6 = prop.individ(YEC6,EC6,SEC6)
# sys.exit()

''' TRILEGAL selection cuts '''
# TRILEGAL_C3 = prop.det_prob(TRILEGAL_C3,'numax',3090.0,135.1)
# TRILEGAL_C6 = prop.det_prob(TRILEGAL_C6,'numax',3090.0,135.1)
# TRILEGAL_C3 = TRILEGAL_C3[ (TRILEGAL_C3['numax'] > 10) & (TRILEGAL_C3['numax'] < 280) & \
# (TRILEGAL_C3['Hmag'] > 7) & (TRILEGAL_C3['Hmag'] < 12) & (TRILEGAL_C3['JK'] > 0.5) & \
# (TRILEGAL_C3['prob_s'] > 0.95)]
# TRILEGAL_C6 = TRILEGAL_C6[ (TRILEGAL_C6['numax'] > 10) & (TRILEGAL_C6['numax'] < 280) & \
# (TRILEGAL_C6['Vcut'] > 9) & (TRILEGAL_C6['Vcut'] < 15) & (TRILEGAL_C6['JK'] > 0.5) & \
# (TRILEGAL_C6['prob_s'] > 0.95)]

# C3orig = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3_full')
# C6orig = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6_full')
# bins = [np.linspace(min(C3orig['mass']),max(C3orig['mass']),50), \
#         np.linspace(min(C3orig['Radius']),max(C3orig['Radius']),50), \
#         np.linspace(min(C3orig['Teff']),max(C3orig['Teff']),50), \
#         np.linspace(min(C3orig['[Fe/H]']),max(C3orig['[Fe/H]']),50), \
#         np.linspace(min(C3orig['logg']),max(C3orig['logg']),50)]
# hist_orig(TRILEGAL_C3,TRILEGAL_C3_cut,'All',bins,'/media/ben/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3/',6)
# bins = [np.linspace(min(C6orig['mass']),max(C6orig['mass']),50), \
#         np.linspace(min(C6orig['Radius']),max(C6orig['Radius']),50), \
#         np.linspace(min(C6orig['Teff']),max(C6orig['Teff']),50), \
#         np.linspace(min(C6orig['[Fe/H]']),max(C6orig['[Fe/H]']),50), \
#         np.linspace(min(C6orig['logg']),max(C6orig['logg']),50)]
# hist_orig(TRILEGAL_C6,TRILEGAL_C6_cut,'All',bins,'/media/ben/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6/',6)

''' Save out TRILEGAL simulations after selection function has been applied
    for use elsewhere '''

# TRI3 = TRILEGAL_C3[TRILEGAL_C3['logAge'] > np.log10(9e9)]
# TRI6 = TRILEGAL_C6[TRILEGAL_C6['logAge'] > np.log10(9e9)]
# TRI3.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/TRILEGAL_C3_old_stars',index=False)
# TRI6.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/TRILEGAL_C6_old_stars',index=False)
# TRILEGAL_C3.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/TRILEGAL_C3_self',index=False)
# TRILEGAL_C6.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/TRILEGAL_C6_self',index=False)
# print('Trilegal saved out')

# plt.figure()
# # plt.scatter(X[numax],X['KepMag'],c=X['prob_s'],cmap=colormaps.parula,alpha=0.5)
# plt.scatter(YC3['nmx'],YC3['KepMag'])
# x = np.where(YC3['det_prob_flag'] == 1)
# plt.scatter(YC3['nmx'].iloc[x],YC3['KepMag'].iloc[x],marker='.')
# plt.gca().invert_yaxis()
# # cbar = plt.colorbar()
# plt.show()
# sys.exit()

''' Yvonne Detects and Benoit doesn't
    - Merge, concatenate, delete duplicates '''
# YB3 = pd.merge(YC3,BC3[['EPIC']],how='inner',on=['EPIC'])
# YB6 = pd.merge(YC6,BC6[['EPIC']],how='inner',on=['EPIC'])
# Y3 = pd.concat([YC3,YB3]).drop_duplicates(subset=['EPIC'],keep=False).reset_index(drop=True)
# Y6 = pd.concat([YC6,YB6]).drop_duplicates(subset=['EPIC'],keep=False).reset_index(drop=True)
# print(len(YB3),len(YB6))
# print(len(Y3),len(Y6))
# print(len(YC3),len(YC6))
#
# Y3.to_csv(ext_GA+'GA/C3_Yvonne_det',index=False,columns=['EPIC','nmx','nmx_err','dnu','dnu_err'])
# Y6.to_csv(ext_GA+'GA/C6_Yvonne_det',index=False,columns=['EPIC','nmx','nmx_err','dnu','dnu_err'])
# sys.exit()



''' Flag any Super-Nyquist stars '''
# seismo_list = [YC3,YC6,SC3,SC6,BC3,BC6,EC6,YEC6,EC3,YEC3,SEC3,SEC6]
# seismo_name = ['YC3','YC6','SC3','SC6','BC3','BC6','EC6','YEC6','EC3','YEC3','SEC3','SEC6']
# dnu = ['YDnu','YDnu','SDnu','SDnu','BDnu','BDnu','EDnu','dnu','EDnu','dnu','SDnu','SDnu']
# # for i in seismo_list:
# #     print( i['Nseismo'][0])
# YC3,YC6,SC3,SC6,BC3,BC6,EC6,YEC6,EC3,YEC3,SEC3,SEC6 = prop.nyq_refl(seismo_list,dnu,30,seismo_name)


''' Merging of GES data with single asteroseismic dets '''
YG3,SG3,BG3,EG3 = dat.GES_merge(seismo3_list,GES3,seismo3_name)
GES = pd.concat([BG3,SG3,YG3],ignore_index=True)
GES = GES.drop_duplicates(subset=['EPIC'])
GES = GES.fillna(value='NaN',method=None)
GES = GES.reset_index(drop=True)
# GES.to_csv(ext_GA+'GA/K2Poles/Gaia_ESO/GES_full.csv',index=False,na_rep='Inf')
# print( "Gaia-ESO saved out")

''' Merging of APOGEE data with single asteroseismic dets '''
YA3,SA3,BA3,EA3 = dat.APO_merge(seismo3_list,APO3,seismo3_name)
AP3 = pd.concat([BA3,SA3,YA3],ignore_index=True)
AP3 = AP3.drop_duplicates(subset=['EPIC'])
AP3 = AP3.fillna(value='NaN',method=None)
AP3 = AP3.reset_index(drop=True)
# AP3.to_csv(ext_GA+'GA/K2Poles/APO_LAMOST/APOGEE_full_C3.csv',index=False,na_rep='Inf')

YA6,SA6,BA6,EA6 = dat.APO_merge(seismo6_list,APO6,seismo6_name)
AP6 = pd.concat([BA6,SA6,YA6],ignore_index=True)
AP6 = AP6.drop_duplicates(subset=['EPIC'])
AP6 = AP6.fillna(value='NaN',method=None)
AP6 = AP6.reset_index(drop=True)
# AP6.to_csv(ext_GA+'GA/K2Poles/APO_LAMOST/APOGEE_full_C6.csv',index=False,na_rep='Inf')
# print( "APOGEE saved out")

''' Merging of LAMOST data with single asteroseismic dets '''
YL3 = pd.merge(YC3,LAMOST3,how='inner',on=['EPIC'])
BL3 = pd.merge(BC3,LAMOST3,how='inner',on=['EPIC'])
SL3 = pd.merge(SC3,LAMOST3,how='inner',on=['EPIC'])
LST = pd.concat([BL3,SL3,YL3],ignore_index=True)
LST = LST.drop_duplicates(subset=['EPIC'])
LST = LST.fillna(value='NaN',method=None)
LST = LST.reset_index(drop=True)
# LST.to_csv(ext_GA+'GA/K2Poles/APO_LAMOST/LAMOST_full_C3.csv',index=False,na_rep='Inf')
# print( "LAMOST C3 saved out ", len(LST))

YL6 = pd.merge(YC6,LAMOST6,how='inner',on=['EPIC'])
BL6 = pd.merge(BC6,LAMOST6,how='inner',on=['EPIC'])
SL6 = pd.merge(SC6,LAMOST6,how='inner',on=['EPIC'])
LAMOST = pd.concat([BL6,SL6,YL6],ignore_index=True)
LAMOST = LAMOST.drop_duplicates(subset=['EPIC'],keep='first')
LAMOST = LAMOST.fillna(value='NaN',method=None)
LAMOST = LAMOST.reset_index(drop=True)
# LAMOST.to_csv(ext_GA+'GA/K2Poles/APO_LAMOST/LAMOST_full_C6.csv',index=False,na_rep='Inf')
# print( "LAMOST C6 saved out ", len(LAMOST))

''' Merging RAVE data with individual asteroseismic data sets for full RAVE
    spectro-seismic data list '''
RAVE_list = [RAVE3,RAVE6,RAVE3,RAVE6,RAVE3,RAVE6]
YR3,YR6,BR3,BR6,SR3,SR6 = dat.RAVE_merge([YC3,YC6,BC3,BC6,SC3,SC6], \
                                RAVE_list,['YR3','YR6','BR3','BR6','SR3','SR6'])

RAVE3 = pd.concat([YR3,BR3,SR3],ignore_index=True)
RAVE3 = RAVE3.drop_duplicates(subset=['EPIC'])
RAVE3 = RAVE3.fillna(value='NaN',method=None)
RAVE3 = RAVE3.reset_index(drop=True)
# RAVE3.to_csv(ext_GA+'GA/K2Poles/RAVE_C3.csv',index=False,na_rep='Inf')
RAVE6 = pd.concat([YR6,BR6,SR6],ignore_index=True)
RAVE6 = RAVE6.drop_duplicates(subset=['EPIC'])
RAVE6 = RAVE6.fillna(value='NaN',method=None)
RAVE6 = RAVE6.reset_index(drop=True)
# RAVE6.to_csv(ext_GA+'GA/K2Poles/RAVE_C6.csv',index=False,na_rep='Inf')
RAVE3 = RAVE3[RAVE3['[Fe/H]_RAVE'] > -900]
RAVE6 = RAVE6[RAVE6['[Fe/H]_RAVE'] > -900]
# print(RAVE3.columns.values)
RAVE3['TEFF'] = RAVE3['Teff_RAVE']
RAVE6['TEFF'] = RAVE6['Teff_RAVE']
# sys.exit()

''' Complete asteroseismic lists '''
camp3_0 = pd.concat([YC3,SC3,BC3],ignore_index=True)
camp3_0 = camp3_0.drop_duplicates(subset=['EPIC'])
camp3_0 = camp3_0.reset_index(drop=True)
camp3_0 = camp3_0.fillna(value='NaN',method=None)

camp6_0 = pd.concat([YC6,SC6,BC6],ignore_index=True)
camp6_0 = camp6_0.drop_duplicates(subset=['EPIC'])
camp6_0 = camp6_0.reset_index(drop=True)
camp6_0 = camp6_0.fillna(value='NaN',method=None)

cols = ['EPIC','RA','Dec','Teff','[Fe/H]','slogg','logg','radius','radius_val']
K2_camp = pd.concat([GAP3,GAP6],ignore_index=True)
K2_camp = K2_camp.reset_index(drop=True)

K2_camp_v2 = pd.concat([camp3_0,camp6_0],ignore_index=True)
K2_camp_v3 = pd.merge(K2_camp_v3,K2_camp_v2[['EPIC']],how='inner',on=['EPIC'])

fig,ax = plt.subplots()
a = np.where(K2_camp_v2['prob_s'] >= 0.95)
b = np.where(K2_camp_v3['prob_s_gaia'] >= 0.95)
ax.scatter(K2_camp_v2['JK'].iloc[a],K2_camp_v2['Kabs'].iloc[a],alpha=0.5,label=r'R$_{\rm{catalogue}}$')
ax.scatter(K2_camp_v3['JK'].iloc[b],K2_camp_v3['Kabs'].iloc[b],alpha=0.5,label=r'R$_{\rm{Gaia}}$')
ax.set_xlabel(r'J - K',fontsize=15)
ax.set_ylabel(r'K$_{\rm{abs}}$',fontsize=15)
ax.set_xlim(0.475,1.325)
ax.invert_yaxis()
ax.legend()
# plt.show()
f = plt.figure()
plt.hist(K2_camp_v3['Rgaia'],bins=np.linspace(0,20,50),histtype='step',label=r'R$_{Gaia}$',normed=True,linewidth=2)
plt.hist(K2_camp_v2['radius'],bins=np.linspace(0,20,50),histtype='step',label=r'R$_{seismo}$',normed=True,linewidth=2)
plt.xlim(3.5,18)
plt.xlabel(r'Radius [R$_{\odot}$]',fontsize=15)
plt.legend(prop={'size':10})
# plt.show()
# sys.exit()

fig,ax = plt.subplots()
a = np.where(K2_camp_v2['prob_s'] >= 0.95)
b = np.where(K2_camp_v3['prob_s_gaia'] >= 0.95)
ax.scatter(K2_camp['JK'],K2_camp['Kabs'],alpha=0.5,label=r'GAP')
ax.scatter(K2_camp_v2['JK'].iloc[a],K2_camp_v2['Kabs'].iloc[a],alpha=0.5,label=r'R$_{\rm{catalogue}}$')
ax.scatter(K2_camp_v3['JK'].iloc[b],K2_camp_v3['Kabs'].iloc[b],alpha=0.5,label=r'R$_{\rm{Gaia}}$')
ax.set_xlabel(r'J - K',fontsize=15)
ax.set_ylabel(r'K$_{\rm{abs}}$',fontsize=15)
ax.invert_yaxis()
ax.legend()
# plt.show()
# fig.savefig('Det_prob_cut_HRD.pdf', bbox_inches='tight')
# sys.exit()

# K2_camp.to_csv('/home/ben/Desktop/GAP_Gaia',columns=cols,index=False)
# camp3_0.to_csv('/home/bmr135/Desktop/C3_gaia',columns=cols,index=False)
# camp6_0.to_csv('/home/bmr135/Desktop/C6_gaia',columns=cols,index=False)

# K2_v3 = pd.merge(K2_camp_v3,camp3_0[['EPIC']],how='inner',on=['EPIC'])
# K2_v6 = pd.merge(K2_camp_v3,camp6_0[['EPIC']],how='inner',on=['EPIC'])
#
# print('Seismic cross with GAP',len(K2_v3),len(K2_v6))
print('GAP lengths:',len(GAP3_v3),len(GAP6_v3))
print('Full seismic:',len(camp3_0),len(camp6_0))
# sys.exit()

''' Complete spectroscopic lists '''
spec3_0 = pd.concat([AP3,RAVE3,GES],ignore_index=True)
spec3_0 = spec3_0.drop_duplicates(subset=['EPIC'])
spec3_0 = spec3_0.reset_index(drop=True)
spec3_0 = spec3_0.fillna(value='NaN',method=None)

spec6_0 = pd.concat([AP6,RAVE6,LAMOST],ignore_index=True)
spec6_0 = spec6_0.drop_duplicates(subset=['EPIC'])
spec6_0 = spec6_0.reset_index(drop=True)
spec6_0 = spec6_0.fillna(value='NaN',method=None)

# spec3_0 = prop.single_seismo(spec3_0,['Bnumax','nmx','Snumax'],'NUMAX')
# spec3_0 = prop.single_seismo(spec3_0,['BDnu','dnu','SDnu'],'DNU')
# spec6_0 = prop.single_seismo(spec6_0,['Bnumax','nmx','Snumax'],'NUMAX')
# spec6_0 = prop.single_seismo(spec6_0,['BDnu','dnu','SDnu'],'DNU')
# spec3_0['sRad'] = (spec3_0['NUMAX']/3090.) * (spec3_0['DNU']/135.1)**-2 * (spec3_0['TEFF']/const.solar_Teff)**0.5
# spec6_0['sRad'] = (spec3_0['NUMAX']/3090.) * (spec3_0['DNU']/135.1)**-2 * (spec3_0['TEFF']/const.solar_Teff)**0.5

cols = ['EPIC','RA','Dec','Teff','[Fe/H]','slogg','logg','radius','radius_val','TEFF','LOGG','sRad']

# spec3_0.to_csv('/home/bmr135/Desktop/C3_spec_gaia',columns=cols,index=False)
# spec6_0.to_csv('/home/bmr135/Desktop/C6_spec_gaia',columns=cols,index=False)
# print(len(spec3_0),len(spec6_0))
# sys.exit()

# YC3.to_csv('/media/ben/SAMSUNG1/GA/K2Poles/YC3_TL',index=False)
# BC3.to_csv('/media/ben/SAMSUNG1/GA/K2Poles/BC3_TL',index=False)
# SC3.to_csv('/media/ben/SAMSUNG1/GA/K2Poles/SC3_TL',index=False)

''' Number of pipeline detections cut '''
YC3 = YC3[YC3['Nseismo'] >= 2]
BC3 = BC3[BC3['Nseismo'] >= 2]
SC3 = SC3[SC3['Nseismo'] >= 2]
YC6 = YC6[YC6['Nseismo'] >= 2]
BC6 = BC6[BC6['Nseismo'] >= 2]
SC6 = SC6[SC6['Nseismo'] >= 2]

''' Merging of multiple fields for comparison -> K2P2 '''
# Yvonne + Savita C3
cols_to_use = YC3.columns.difference(SC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YS_C3 = pd.merge(SC3,YC3[cols_to_use],how='inner',on=['EPIC'])
# Yvonne + Savita C6
cols_to_use = YC6.columns.difference(SC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YS_C6 = pd.merge(SC6,YC6[cols_to_use],how='inner',on=['EPIC'])
# Yvonne + Benoit C3
cols_to_use = YC3.columns.difference(BC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YB_C3 = pd.merge(BC3,YC3[cols_to_use],how='inner',on=['EPIC'])
# Yvonne + Benoit C6
cols_to_use = YC6.columns.difference(BC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YB_C6 = pd.merge(BC6,YC6[cols_to_use],how='inner',on=['EPIC'])
# Savita + Benoit C3
cols_to_use = SC3.columns.difference(BC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
BS_C3 = pd.merge(BC3,SC3[cols_to_use],how='inner',on=['EPIC'])
# Savita + Benoit C6
cols_to_use = SC6.columns.difference(BC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
BS_C6 = pd.merge(BC6,SC6[cols_to_use],how='inner',on=['EPIC'])

'''  Merging of multiple fields for comparison -> K2P2/EVEREST cross-match '''
cols_to_use = EC6.columns.difference(BC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
EB_C6 = pd.merge(BC6,EC6[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = EC3.columns.difference(YEC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
EE_C3 = pd.merge(YEC3,EC3[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = EC6.columns.difference(YEC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
EE_C6 = pd.merge(YEC6,EC6[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = EC3.columns.difference(BC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
EB_C3 = pd.merge(BC3,EC3[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = YEC3.columns.difference(YC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YY_C3 = pd.merge(YC3,YEC3[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = YEC6.columns.difference(YC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YY_C6 = pd.merge(YC6,YEC6[cols_to_use],how='inner',on=['EPIC'])

''' Merging of multiple fields for comparison -> EVERST '''
cols_to_use = SEC3.columns.difference(SC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
SS_C3 = pd.merge(SC3,SEC3[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = SEC6.columns.difference(SC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
SS_C6 = pd.merge(SC6,SEC6[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = SEC3.columns.difference(YEC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YSE_C3 = pd.merge(YEC3,SEC3[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = SEC6.columns.difference(YEC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YSE_C6 = pd.merge(YEC6,SEC6[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = SEC3.columns.difference(EC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
SE_C3 = pd.merge(EC3,SEC3[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = SEC6.columns.difference(YEC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
SE_C6 = pd.merge(EC6,SEC6[cols_to_use],how='inner',on=['EPIC'])

print( "Merging Complete")


''' Refining merged data sets so that only values of numax/dnu within 2 sigma
    of each other are maintained '''
merged = [YS_C3,YS_C6,YB_C3,YB_C6,BS_C3,BS_C6,EB_C3,EB_C6,EE_C3,EE_C6,YY_C3,YY_C6,SS_C3,SS_C6, \
            YSE_C3,YSE_C6,SE_C3,SE_C6]
nmx1 = ['nmx','nmx','nmx','nmx','Bnumax','Bnumax','Enumax','Enumax','Enumax', \
            'Enumax','nmx','nmx','Snumax','Snumax','nmx','nmx','Snumax','Snumax']
err_nmx1 = ['nmx_err','nmx_err','nmx_err','nmx_err','e_Bnumax','e_Bnumax','e_Enumax', \
            'e_Enumax','e_Enumax','e_Enumax','nmx_err','nmx_err','e_Snumax','e_Snumax', \
            'nmx_err','nmx_err','e_Snumax','e_Snumax']
nmx2 = ['Snumax','Snumax','Bnumax','Bnumax','Snumax','Snumax','Bnumax','Bnumax','nmx','nmx', \
        'nmx','nmx','Snumax','Snumax','Snumax','Snumax','Enumax','Enumax']
err_nmx2 = ['e_Snumax','e_Snumax','e_Bnumax','e_Bnumax','e_Snumax','e_Snumax','e_Bnumax', \
            'e_Bnumax','nmx_err','nmx_err','nmx_err','nmx_err','e_Snumax','e_Snumax', \
            'e_Snumax','e_Snumax','e_Enumax','e_Enumax']
dnu1 = ['dnu','dnu','dnu','dnu','BDnu','BDnu','EDnu','EDnu','EDnu','EDnu','dnu', \
        'dnu','SDnu','SDnu','dnu','dnu','SDnu','SDnu']
err_dnu1 = ['dnu_err','dnu_err','dnu_err','dnu_err','e_BDnu','e_BDnu','e_EDnu','e_EDnu', \
            'e_EDnu','e_EDnu','dnu_err','dnu_err','e_SDnu','e_SDnu','dnu_err','dnu_err', \
            'e_SDnu','e_SDnu']
dnu2 = ['SDnu','SDnu','BDnu','BDnu','SDnu','SDnu','BDnu','BDnu','dnu','dnu','dnu','dnu', \
        'SDnu','SDnu','dnu','dnu','SDnu','SDnu']
err_dnu2 = ['e_SDnu','e_SDnu','e_BDnu','e_BDnu','e_SDnu','e_SDnu','e_BDnu','e_BDnu', \
            'dnu_err','dnu_err','dnu_err','dnu_err','e_SDnu','e_SDnu','e_SDnu','e_SDnu', \
            'e_EDnu','e_EDnu']

''' Sigma Clip '''
for i in range(0,len(merged),1):
    merged[i] = prop.sigma_clip(merged[i],nmx1[i],nmx2[i],dnu1[i],dnu2[i],err_nmx1[i],err_nmx2[i],err_dnu1[i],err_dnu2[i],3)

''' Output of sigma clip with flag used instead of a straight cut. Maintains full
    merged datasets for usage in comparison tests. '''
YS_C3,YS_C6,YB_C3,YB_C6,BS_C3,BS_C6,EB_C3,EB_C6,EE_C3,\
EE_C6,YY_C3,YY_C6,SS_C3,SS_C6,YSE_C3,YSE_C6,SE_C3,SE_C6 = merged

''' Implementation of the sigma clip to reduce the samples '''
for i in range(0,len(merged),1):
    merged[i] = merged[i][(merged[i]['sig_clip_flag'] == 1)]

YS_C3,YS_C6,YB_C3,YB_C6,BS_C3,BS_C6,EB_C3,EB_C6,EE_C3,\
EE_C6,YY_C3,YY_C6,SS_C3,SS_C6,YSE_C3,YSE_C6,SE_C3,SE_C6 = merged

print( "Sigma clip complete")

# plt.figure()
# plt.hist(camp3_0['[Fe/H]'],bins=50,normed=True,label=r'K2 C3')
# plt.hist(camp6_0['[Fe/H]'],bins=50,normed=True,label=r'K2 C6')
# plt.hist(RAVE3['[Fe/H]_RAVE'],bins=50,normed=True,label=r'RAVE 3')
# plt.hist(RAVE6['[Fe/H]_RAVE'],bins=50,normed=True,label=r'RAVE 6')
# plt.xlabel(r'[Fe/H]')
# yt = plt.gca()
# yt.axes.yaxis.set_ticks([])
# yt.axes.yaxis.set_ticklabels([])
# plt.legend()
# plt.show()

# sys.exit()


APORAVE = pd.merge(AP6,RAVE6,how='inner',on=['EPIC'])
# print( len(APORAVE))

''' Concatenated lists of stars in each campaign -> K2P2.
    Merging code similar to this and flag generation still required '''
camp3 = pd.concat([YB_C3,YS_C3,BS_C3],ignore_index=True)
camp3 = camp3.drop_duplicates(subset=['EPIC'])
camp3 = camp3.reset_index(drop=True)
camp3 = camp3.fillna(value='NaN',method=None)

camp3 = prop.single_seismo(camp3,['Bnumax','nmx','Snumax'],'NUMAX')
camp3 = prop.single_seismo(camp3,['e_Bnumax','nmx_err','e_Snumax'],'NUMAX_err')
camp3 = prop.single_seismo(camp3,['BDnu','dnu','SDnu'],'DNU')
camp3 = prop.single_seismo(camp3,['e_BDnu','dnu_err','e_SDnu'],'DNU_err')
# camp3 = prop.met_filter(camp3)
camp6 = pd.concat([YB_C6,YS_C6,BS_C6],ignore_index=True)
camp6 = camp6.drop_duplicates(subset=['EPIC'])
camp6 = camp6.reset_index(drop=True)
camp6 = camp6.fillna(value='NaN',method=None)

camp6 = prop.single_seismo(camp6,['Bnumax','nmx','Snumax'],'NUMAX')
camp6 = prop.single_seismo(camp6,['e_Bnumax','nmx_err','e_Snumax'],'NUMAX_err')
camp6 = prop.single_seismo(camp6,['BDnu','dnu','SDnu'],'DNU')
camp6 = prop.single_seismo(camp6,['e_BDnu','dnu_err','e_SDnu'],'DNU_err')
# camp6 = prop.met_filter(camp6)

print('C3/C6 lengths:',len(camp3),len(camp6))

Luca = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_Gaia_BC_full.csv')
a = Luca[Luca['#Id'] < 208000000]
b = Luca[Luca['#Id'] > 208000000]
print(len(a),len(b))
Luca = Luca.rename(columns={'#Id':'EPIC'})


a = pd.merge(Luca,camp3,how='inner',on=['EPIC'])
b = pd.merge(Luca,camp6,how='inner',on=['EPIC'])
print('C3/C6 Luca cross:',len(a),len(b))
# sys.exit()

''' Data Flag for EPIC parametric values - which [Fe/H] to use? '''
# spectro_EPICS_3 = camp3[camp3['stpropflag'] != 'rpm']
# spectro_EPICS_6 = camp6[camp6['stpropflag'] != 'rpm']
# spectro_EPICS_3 = spectro_EPICS_3[spectro_EPICS_3['stpropflag'] != 'rav']
# spectro_EPICS_6 = spectro_EPICS_6[spectro_EPICS_6['stpropflag'] != 'rav']
# print(len(spectro_EPICS_6['stpropflag']))

# GAP_camp3 = pd.merge(GAP3_v2,camp3,how='inner',on=['EPIC'])
# GAP_camp6 = pd.merge(GAP6_v2,camp6,how='inner',on=['EPIC'])
# print(len(GAP_camp3),len(GAP_camp6))


''' Merging of RAVE and full asteroseismic data '''
cols_to_use = camp3.columns.difference(RAVE3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
RC3 = pd.merge(RAVE3,camp3[cols_to_use],how='inner',on=['EPIC'])
RC3['slogg_spec'] = np.log10(const.solar_g * (RC3['NUMAX'].astype('float64')/3090) * np.sqrt(RC3['Teff_RAVE'].astype('float64')/const.solar_Teff))
cols_to_use = camp6.columns.difference(RAVE6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
RC6 = pd.merge(RAVE6,camp6[cols_to_use],how='inner',on=['EPIC'])
RC6['slogg_spec'] = np.log10(const.solar_g * (RC6['NUMAX'].astype('float64')/3090.0) * np.sqrt(RC6['Teff_RAVE'].astype('float64')/const.solar_Teff))
print("RAVE: ", len(RC3), len(RC6))

# sys.exit()

''' Merging of GES data with multiple asteroseismic dets '''
cols_to_use = camp3.columns.difference(GES3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
GES = pd.merge(GES3,camp3[cols_to_use],how='inner',on=['EPIC'])
GES['slogg_spec'] = np.log10(const.solar_g * (GES['NUMAX'].astype('float64')/3090.) * np.sqrt(GES['TEFF'].astype('float64')/const.solar_Teff))

''' Alteration of alpha abundances in Gaia-ESO '''
# GES['ALPHA'] = GES['ALPHA'] + 0.1
# GES.to_csv(ext_GA+'GA/K2Poles/matlab_in/GES_p0.1_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# GES['ALPHA'] = GES['ALPHA'] + 0.15
# GES.to_csv(ext_GA+'GA/K2Poles/matlab_in/GES_p0.25_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# GES['ALPHA'] = GES['ALPHA'] - 0.35
# GES.to_csv(ext_GA+'GA/K2Poles/matlab_in/GES_m0.1_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# GES['ALPHA'] = GES['ALPHA'] -0.15
# GES.to_csv(ext_GA+'GA/K2Poles/matlab_in/GES_0.m25_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
print("Gaia-ESO: ", len(GES))

''' Stars in C3 with no spectra --> Using for WHT-ISIS proposal '''
# C3_nospec = pd.DataFrame()
# C3_nospec = pd.concat([camp3,RC3,GES],ignore_index=True)
# C3_nospec = C3_nospec.drop_duplicates(subset=['EPIC'],keep=False)
# # C6_nospec = C6_nospec.dropna()
# df = C3_nospec[['EPIC','RA','Dec','Vmag','Bmag','rmag','Jmag','Hmag','Kmag','JK','BV']]
# df = df[df.Vmag != 'NaN']
# # df = df[df.values != 'NaN']
# df = df.reset_index(drop=True)
# # df.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/C3_TL',index=False)
# plt.figure()
# # hist, bins = np.histogram(df['Vmag'],bins=[9,10,11,12,13,14,15])
# # print(hist,bins)
# plt.hist(df['Vmag'],bins=[9,10,11,12,13,14,15])
# plt.xlabel(r'V',fontsize=15)
# plt.show()
# # sys.exit()

''' Merging of LAMOST data with multiple asteroseismic dets '''
cols_to_use = camp3.columns.difference(LAMOST3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
L3 = pd.merge(LAMOST3,camp3[cols_to_use],how='inner',on=['EPIC'])
L3['slogg_spec'] = np.log10(const.solar_g * (L3['NUMAX'].astype('float64')/3090.) * np.sqrt(L3['teff_L'].astype('float64')/const.solar_Teff))
cols_to_use = camp6.columns.difference(LAMOST6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
L6 = pd.merge(LAMOST6,camp6[cols_to_use],how='inner',on=['EPIC'])
L6['slogg_spec'] = np.log10(const.solar_g * (L6['NUMAX'].astype('float64')/3090.) * np.sqrt(L6['teff_L'].astype('float64')/const.solar_Teff))
print("LAMOST: ", len(LAMOST3), len(LAMOST6))

''' Merging of APOGEE data with multiple asteroseismic dets '''
cols_to_use = camp3.columns.difference(APO3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
AP3 = pd.merge(APO3,camp3[cols_to_use],how='inner',on=['EPIC'])
AP3['slogg_spec'] = np.log10(const.solar_g * (AP3['NUMAX'].astype('float64')/3090.) * np.sqrt(AP3['TEFF'].astype('float64')/const.solar_Teff))
# print(AP3.columns.values)
# sys.exit()

cols_to_use = camp6.columns.difference(APO6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
AP6 = pd.merge(APO6,camp6[cols_to_use],how='inner',on=['EPIC'])
AP6['slogg_spec'] = np.log10(const.solar_g * (AP6['NUMAX'].astype('float64')/3090.) * np.sqrt(AP6['TEFF'].astype('float64')/const.solar_Teff))
print("APOGEE: ", len(AP3), len(AP6))

''' Complete spectroscopic lists '''
spec3 = pd.concat([AP3,RC3,GES],ignore_index=True)
spec3 = spec3.drop_duplicates(subset=['EPIC'])
spec3 = spec3.reset_index(drop=True)
spec3 = spec3.fillna(value='NaN',method=None)

spec6 = pd.concat([AP6,RC6,L6],ignore_index=True)
spec6 = spec6.drop_duplicates(subset=['EPIC'])
spec6 = spec6.reset_index(drop=True)
spec6 = spec6.fillna(value='NaN',method=None)

print('C3/C6 Spec. lengths:',len(spec3),len(spec6))

''' Spectroscopic Sample Overlaps '''
print(len(pd.merge(AP3,GES,how='inner',on=['EPIC'])))
print(len(pd.merge(AP3,RC3,how='inner',on=['EPIC'])))
print(len(pd.merge(AP3,L3,how='inner',on=['EPIC'])))
print(len(pd.merge(RC3,L3,how='inner',on=['EPIC'])))
print(len(pd.merge(GES,L3,how='inner',on=['EPIC'])))
print(len(pd.merge((pd.merge(GES,RC3,how='inner',on=['EPIC'])),AP3,how='inner',on=['EPIC'])))

print(len(pd.merge(AP6,RC6,how='inner',on=['EPIC'])))
print(len(pd.merge(L6,AP6,how='inner',on=['EPIC'])))
print(len(pd.merge(L6,RC6,how='inner',on=['EPIC'])))

''' Merge of LAMOST and RAVE '''
# cols_to_use = RC6.columns.difference(L6.columns)
# cols_to_use = cols_to_use.union(['EPIC'])
# LR6 = pd.merge(L6,RC6[cols_to_use],how='inner',on=['EPIC'])
# print('LAMOST/RAVE C6 saved out')

# print(RC3.columns.values)
# sys.exit()
AP6 = prop.alt_spec_params(AP6,0,'RAVE',['TEFF','LOGG','FE_H'],['TEFF_ERR','LOGG_ERR','FE_H_ERR'])
AP3 = prop.alt_spec_params(AP3,1,'GES',['TEFF','LOGG','FE_H'],['TEFF_ERR','LOGG_ERR','FE_H_ERR'])
RC3 = prop.alt_spec_params(RC3,1,'GES',['Teff_RAVE','logg_RAVE','[Fe/H]_RAVE'],['sig_Teff','sig_logg','sig_feh'])

# sys.exit()
''' Save out combined data sets to be processed for use with PARAM '''
# camp3.to_csv(ext_GA+'GA/K2Poles/matlab_in/C3_'+time.strftime("%d%m%Y")+'.csv',index=False)
# camp6.to_csv(ext_GA+'GA/K2Poles/matlab_in/C6_'+time.strftime("%d%m%Y")+'.csv',index=False)
# RC3.to_csv(ext_GA+'GA/K2Poles/matlab_in/RC3_GES_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# RC6.to_csv(ext_GA+'GA/K2Poles/matlab_in/RC6_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# GES.to_csv(ext_GA+'GA/K2Poles/matlab_in/GES_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# L3.to_csv(ext_GA+'GA/K2Poles/matlab_in/LAMOST3_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# L6.to_csv(ext_GA+'GA/K2Poles/matlab_in/LAMOST6_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# AP3.to_csv(ext_GA+'GA/K2Poles/matlab_in/AP3_GES_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# AP6.to_csv(ext_GA+'GA/K2Poles/matlab_in/AP6_RAVE_'+time.strftime("%d%m%Y")+'.csv',index=False,na_rep='Inf')
# sys.exit()

''' Check length of data sets '''
# print( len(BC3), len(YC3), len(SC3), len(BC6), len(YC6), len(SC6))
# print( len(BR3), len(YR3), len(SR3), len(BR6), len(YR6), len(SR6))
# print( len(YEC3), len(YC6), len(EE_C6))
# print( C3.columns.values)

print( "Processing Complete")

''' PLOTTING '''

''' Gaia-ESO log-g iteration calculations '''
# SG3['LOGG_Seis'] = np.log10(const.solar_g*(SG3['Snumax']/const.solar_Numax_s)*np.sqrt(SG3['TEFF']/const.solar_Teff))
# BG3['LOGG_Seis'] = np.log10(const.solar_g*(BG3['Bnumax']/const.solar_Numax)*np.sqrt(BG3['TEFF']/const.solar_Teff))
# GES = pd.concat([BG3,YG3,SG3],ignore_index=True)
# GES = GES.drop_duplicates(subset=['EPIC'])
# GES = GES.reset_index(drop=True)
# # print( GES.columns.values)
# # print( GES['EPIC'])
#
# ges = pd.DataFrame()
# ges['EPIC'] = GES['EPIC']
# ges['LOGG_Seis'] = GES['LOGG_Seis']
# ges['TEFF2'] = GES['TEFF']
# ges = ges.drop_duplicates(subset=['EPIC'])
# # ges['SPEC_FEH_ispec'] = GES['SPEC_FEH_ispec']
# # print( GES['TEFF'])
#
# for i in range(0,len(ges['EPIC'])):
#     # if ges['SPEC_TEFF_ispec'][i] == GES['SPEC_TEFF_ispec'][i]:
#     #     print( "Y ", ges['EPIC'][i])
#     if ges['TEFF2'][i] != GES['TEFF'][i]:
#         print( i , ges['EPIC'][i])
# # print( len(ges))
#
# out_file = pd.read_csv(ext_DB+'Dropbox/GES-K2/Diane_Feuillet/epinarbo_ite2_new.txt')#,delim_whitespace=True)
# print( len(out_file))
# # out_file.drop(['SPEC_TEFF_ispec','SPEC_FEH_ispec','EPIC','LOGG_Seis'],inplace=True,axis=1)
# out_file['EPIC'] = out_file['OBJECT'].map(lambda x: x.split('_')[-1])
# out_file['EPIC'] = out_file['EPIC'].convert_objects(convert_numeric=True)
# out = pd.merge(out_file,ges,how='inner',on=['EPIC'])
# out = out.reset_index(drop=True)
# # print( len(out))
# # print( out['TEFF_y'])
#
# for i in range(0,len(ges['EPIC'])):
#     # if out_file['SPEC_TEFF_ispec'][i] == out_file['SPEC_TEFF_ispec'][i]:
#     #     print( "Y ", out_file['EPIC'][i])
#     for j in range(0,len(out['EPIC'])):
#         if ges['EPIC'][i] == out['EPIC'][j]:
#             print( "Teff ", i, j, ges['EPIC'][i], out['EPIC'][j], ges['TEFF2'][i], out['TEFF'][j])
#         # if ges['SPEC_FEH_ispec'][i] == out['SPEC_FEH_ispec'][j]:
#         #     print( "FeH", i, j, ges['EPIC'][i], out['EPIC'][j], ges['SPEC_FEH_ispec'][i], out['SPEC_FEH_ispec'][j])
#
# out = out.fillna(value='NaN',method=None)
# out.to_csv(ext_GA+'Dropbox/GES-K2/Diane_Feuillet/seismo_2.txt',index=False,)#,sep='\t')
# in1 = pd.read_csv(ext_GA+'Dropbox/GES-K2/Diane_Feuillet/epinarbo_ite2_new.txt')#,delimiter='\t')
# # print( len(in1))
# in1['EPIC'] = in1['OBJECT'].map(lambda x: x.split('_')[-1])
# in1['EPIC'] = in1['EPIC'].convert_objects(convert_numeric=True)
# for i in range(0,len(in1['EPIC'])):
#     # if out_file['SPEC_TEFF_ispec'][i] == out_file['SPEC_TEFF_ispec'][i]:
#     #     print( "Y ", out_file['EPIC'][i])
#     for j in range(0,len(out['EPIC'])):
#         if (in1['EPIC'][i] == out['EPIC'][j]): #& (in1['TEFF'][i] == out['TEFF'][j]):
#             print( "Teff ", i, j, in1['EPIC'][i], out['EPIC'][j], in1['TEFF'][i] - out['TEFF'][j], in1['ALPHA_FE'][i] - out['ALPHA_FE'][j])
#         # if in1['EPIC'][i] == out['EPIC'][j]:
#         #     print( "FeH", i, j, in1['EPIC'][i], out['EPIC'][j], in1['SPEC_FEH_ispec'][i] - out['SPEC_FEH_ispec'][j])
#
#
# print( out_file['ID'])

''' CMD/HRD comp plots '''
# C3 CMD/HRDs
# plt.figure()
# plt.subplot(3,2,1)
# k2p.plot_CMD(C3,GAP3,YC3,lab=[r'C3',r'K2 GAP',r'Yvonne'],a='Ynumax')
# plt.subplot(3,2,2)
# k2p.plot_HRD(YC3,TRILEGAL_C3,lab=[r'Yvonne',r'TRILEGAL'])
# plt.subplot(3,2,3)
# k2p.plot_CMD(C3,GAP3,SC3,lab=[r'C3',r'K2 GAP',r'Savita'],a='Snumax')
# plt.subplot(3,2,4)
# k2p.plot_HRD(SC3,TRILEGAL_C3,lab=[r'Savita',r'TRILEGAL'])
# plt.subplot(3,2,5)
# k2p.plot_CMD(C3,GAP3,BC3,lab=[r'C3',r'K2 GAP',r'Benoit'],a='Bnumax')
# plt.subplot(3,2,6)
# k2p.plot_HRD(BC3,TRILEGAL_C3,lab=[r'Benoit',r'TRILEGAL'])
# plt.tight_layout()
#
# # C6 CMD/HRDs
# plt.figure()
# plt.subplot(3,2,1)
# k2p.plot_CMD(C6,GAP6,YC6,lab=[r'C6',r'K2 GAP',r'Yvonne'],a='Ynumax')
# plt.subplot(3,2,2)
# k2p.plot_HRD(YC6,TRILEGAL_C6,lab=[r'Yvonne',r'TRILEGAL'])
# plt.subplot(3,2,3)
# k2p.plot_CMD(C6,GAP6,SC6,lab=[r'C3',r'K2 GAP',r'Savita'],a='Snumax')
# plt.subplot(3,2,4)
# k2p.plot_HRD(SC6,TRILEGAL_C6,lab=[r'Savita',r'TRILEGAL'])
# plt.subplot(3,2,5)
# k2p.plot_CMD(C6,GAP6,BC6,lab=[r'C3',r'K2 GAP',r'Benoit'],a='Bnumax')
# plt.subplot(3,2,6)
# k2p.plot_HRD(BC6,TRILEGAL_C6,lab=[r'Benoit',r'TRILEGAL'])
# plt.tight_layout()

''' N-sigma scatter plots '''
# k2p.plt_comp_scat(YB_C6,a='YDnu',b='BDnu',a1='e_YDnu',b1='e_BDnu', \
# lab=[r'(Yvonne - Benoit)/$\sigma_{Y,B}$',r'Yvonne $\Delta\nu$ [$\mu$Hz]'], \
# p=0.23,u=5,d=4)

''' Comparative histogram '''
# k2p.comp_histo(8.5,15.5,71,dat=[YEC6,EC6,EE_C6],a1='Vmag',a2='Vmag_y', \
#              lab=[r'Y-Everest',r'B-Everest',r'Overlap',r'C6: V',r'Frequency'])

''' Regression score between seismic values from different pipelines '''
# print( r2_score(YB_C6['BDnu'],YB_C6['YDnu']))

''' Report Plots: 1 to 1 scatters '''
# mpar, cpar, empar, ecpar = prop.least_squares_2err(EB_C6)
# x = np.arange(0,280,0.1)
#
# plt.figure()
# plt.subplot(1,2,1)
# # plt.scatter(EB_C6['Bnumax'],EB_C6['Enumax'])
# plt.errorbar(EB_C6['Bnumax'],EB_C6['Enumax'], xerr=[EB_C6['e_Bnumax'],EB_C6['e_Bnumax']], yerr=[EB_C6['e_Enumax'],EB_C6['e_Enumax']], fmt='o', alpha=0.5)
# plt.xlabel(r'KASOC $\nu_{\rm{max}}$')
# plt.ylabel(r'K2P2 $\nu_{\rm{max}}$')
# # plt.plot([0,250],[0,250],color='r')
# plt.xlim(0,250)
# plt.ylim(0,250)
# plt.plot([0,max(x)],[0,max(x)],c='r')
# plt.subplot(1,2,2)
# plt.scatter(EB_C6['BDnu'],EB_C6['EDnu'])
# plt.xlabel(r'K2P2 $\Delta\nu$')
# plt.ylabel(r'Everest $\Delta\nu$')
# plt.plot([0,21],[0,21],color='r')
# plt.xlim(0,21)
# plt.ylim(0,21)
# plt.plot(x,(x*mpar[0])+cpar[0],color='r',linestyle='--')
# plt.fill_between(x,(x*(mpar[0]-empar[0]))+(cpar[0]-ecpar[0]),(x*(mpar[0]+empar[0]))+(cpar[0]+ecpar[0]), \
#                 facecolor='g',edgecolor='green', alpha=0.5)
# plt.tight_layout()
#
# plt.figure()
# plt.subplot(1,2,1)
# plt.scatter(EE_C6['nmx'],EE_C6['Enumax'])
# plt.scatter(Z9['nmx'],Z9['Enumax'],color='grey',alpha=0.50)
# plt.xlabel(r'Elsworth $\nu_{\rm{max}}$')
# plt.ylabel(r'Mosser $\nu_{\rm{max}}$')
# plt.plot([0,250],[0,250],color='r')
# plt.xlim(0,250)
# plt.ylim(0,250)
# plt.subplot(1,2,2)
# plt.scatter(EE_C6['dnu'],EE_C6['EDnu'])
# plt.scatter(Z9['dnu'],Z9['EDnu'],color='grey',alpha=0.50)
# plt.xlabel(r'Elsworth $\Delta\nu$')
# plt.ylabel(r'Mosser $\Delta\nu$')
# plt.plot([0,21],[0,21],color='r')
# plt.xlim(0,21)
# plt.ylim(0,21)
# plt.tight_layout()
#
# plt.figure()
# plt.subplot(1,2,1)
# plt.scatter(EB_C3['Bnumax'],EB_C3['Enumax'])
# plt.scatter(Z6['Bnumax'],Z6['Enumax'],color='grey',alpha=0.50)
# plt.xlabel(r'K2P2 $\nu_{\rm{max}}$')
# plt.ylabel(r'Everest $\nu_{\rm{max}}$')
# plt.plot([0,250],[0,250],color='r')
# plt.xlim(0,250)
# plt.ylim(0,250)
# plt.subplot(1,2,2)
# plt.scatter(EB_C3['BDnu'],EB_C3['EDnu'])
# plt.scatter(Z6['BDnu'],Z6['EDnu'],color='grey',alpha=0.50)
# plt.xlabel(r'K2P2 $\Delta\nu$')
# plt.ylabel(r'Everest $\Delta\nu$')
# plt.plot([0,21],[0,21],color='r')
# plt.xlim(0,21)
# plt.ylim(0,21)
# plt.tight_layout()
#
# plt.figure()
# plt.subplot(1,2,1)
# plt.scatter(YB_C3['Bnumax'],YB_C3['Ynumax'])
# plt.scatter(Z2['Bnumax'],Z2['Ynumax'],color='grey',alpha=0.50)
# plt.xlabel(r'Mosser $\nu_{\rm{max}}$')
# plt.ylabel(r'Elsworth $\nu_{\rm{max}}$')
# plt.plot([0,250],[0,250],color='r')
# plt.xlim(0,250)
# plt.ylim(0,250)
# plt.subplot(1,2,2)
# plt.scatter(YB_C3['BDnu'],YB_C3['YDnu'])
# plt.scatter(Z2['BDnu'],Z2['YDnu'],color='grey',alpha=0.50)
# plt.xlabel(r'Mosser $\Delta\nu$')
# plt.ylabel(r'Elsworth $\Delta\nu$')
# plt.plot([0,21],[0,21],color='r')
# plt.xlim(0,21)
# plt.ylim(0,21)
# plt.tight_layout()
#
# plt.figure()
# plt.subplot(1,2,1)
# plt.scatter(BS_C3['Bnumax'],BS_C3['Snumax'])
# plt.scatter(Z4['Bnumax'],Z4['Snumax'],color='grey',alpha=0.50)
# plt.xlabel(r'Mosser $\nu_{\rm{max}}$')
# plt.ylabel(r'Mathur $\nu_{\rm{max}}$')
# plt.plot([0,250],[0,250],color='r')
# plt.xlim(0,250)
# plt.ylim(0,250)
# plt.subplot(1,2,2)
# plt.scatter(BS_C3['BDnu'],BS_C3['SDnu'])
# plt.scatter(Z4['BDnu'],Z4['SDnu'],color='grey',alpha=0.50)
# plt.xlabel(r'Mosser $\Delta\nu$')
# plt.ylabel(r'Mathur $\Delta\nu$')
# plt.plot([0,21],[0,21],color='r')
# plt.xlim(0,21)
# plt.ylim(0,21)
# plt.tight_layout()
#
# plt.figure()
# plt.subplot(1,2,1)
# plt.scatter(YS_C3['Ynumax'],YS_C3['Snumax'])
# plt.scatter(Z0['Ynumax'],Z0['Snumax'],color='grey',alpha=0.50)
# plt.xlabel(r'Elsworth $\nu_{\rm{max}}$')
# plt.ylabel(r'Mathur $\nu_{\rm{max}}$')
# plt.plot([0,250],[0,250],color='r')
# plt.xlim(0,250)
# plt.ylim(0,250)
# plt.subplot(1,2,2)
# plt.scatter(YS_C3['YDnu'],YS_C3['SDnu'])
# plt.scatter(Z0['YDnu'],Z0['SDnu'],color='grey',alpha=0.50)
# plt.xlabel(r'Elsworth $\Delta\nu$')
# plt.ylabel(r'Mathur $\Delta\nu$')
# plt.plot([0,21],[0,21],color='r')
# plt.xlim(0,21)
# plt.ylim(0,21)
# plt.tight_layout()
#
# plt.figure()
# plt.subplot(1,2,1)
# plt.scatter(YB_C6['Bnumax'],YB_C6['Ynumax'])
# plt.scatter(Z3['Bnumax'],Z3['Ynumax'],color='grey',alpha=0.50)
# plt.xlabel(r'Mosser $\nu_{\rm{max}}$')
# plt.ylabel(r'Elsworth $\nu_{\rm{max}}$')
# plt.plot([0,250],[0,250],color='r')
# plt.xlim(0,250)
# plt.ylim(0,250)
# plt.subplot(1,2,2)
# plt.scatter(YB_C6['BDnu'],YB_C6['YDnu'])
# plt.scatter(Z3['BDnu'],Z3['YDnu'],color='grey',alpha=0.50)
# plt.xlabel(r'Mosser $\Delta\nu$')
# plt.ylabel(r'Elsworth $\Delta\nu$')
# plt.plot([0,21],[0,21],color='r')
# plt.xlim(0,21)
# plt.ylim(0,21)
# plt.tight_layout()
#
# plt.figure()
# plt.subplot(1,2,1)
# plt.scatter(BS_C6['Bnumax'],BS_C6['Snumax'])
# plt.scatter(Z5['Bnumax'],Z5['Snumax'],color='grey',alpha=0.50)
# plt.xlabel(r'Mosser $\nu_{\rm{max}}$')
# plt.ylabel(r'Mathur $\nu_{\rm{max}}$')
# plt.plot([0,250],[0,250],color='r')
# plt.xlim(0,250)
# plt.ylim(0,250)
# plt.subplot(1,2,2)
# plt.scatter(BS_C6['BDnu'],BS_C6['SDnu'])
# plt.scatter(Z5['BDnu'],Z5['SDnu'],color='grey',alpha=0.50)
# plt.xlabel(r'Mosser $\Delta\nu$')
# plt.ylabel(r'Mathur $\Delta\nu$')
# plt.plot([0,21],[0,21],color='r')
# plt.xlim(0,21)
# plt.ylim(0,21)
# plt.tight_layout()
#
# plt.figure()
# plt.subplot(1,2,1)
# plt.scatter(YS_C6['Ynumax'],YS_C6['Snumax'])
# plt.scatter(Z1['Ynumax'],Z1['Snumax'],color='grey',alpha=0.50)
# plt.xlabel(r'Elsworth $\nu_{\rm{max}}$')
# plt.ylabel(r'Mathur $\nu_{\rm{max}}$')
# plt.plot([0,250],[0,250],color='r')
# plt.xlim(0,250)
# plt.ylim(0,250)
# plt.subplot(1,2,2)
# plt.scatter(YS_C6['YDnu'],YS_C6['SDnu'])
# plt.scatter(Z1['YDnu'],Z1['SDnu'],color='grey',alpha=0.50)
# plt.xlabel(r'Elsworth $\Delta\nu$')
# plt.ylabel(r'Mathur $\Delta\nu$')
# plt.plot([0,21],[0,21],color='r')
# plt.xlim(0,21)
# plt.ylim(0,21)
# plt.tight_layout()

''' Light Curve Comp Plots: N - sigma plots '''
# Everest vs KASOC - Benoit
# k2p.plt_comp_scat(EB_C6,'Bnumax','Enumax','e_Bnumax','e_Enumax',[r'$(\nu_{\rm{max},EV}-\nu_{\rm{max},K2P2})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.47,4,4)
# k2p.plt_comp_scat(EB_C3,'Bnumax','Enumax','e_Bnumax','e_Enumax',[r'$(\nu_{\rm{max},EV}-\nu_{\rm{max},K2P2})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.42,4,4)
# k2p.plt_comp_scat(EB_C6,'BDnu','EDnu','e_BDnu','e_EDnu',[r'$(\Delta\nu_{EV}-\Delta\nu_{K2P2})/\sigma_{comp}$',r'$\Delta\nu$'],0.25,5,5)
# k2p.plt_comp_scat(EB_C3,'BDnu','EDnu','e_BDnu','e_EDnu',[r'$(\Delta\nu_{EV}-\Delta\nu_{K2P2})/\sigma_{comp}$',r'$\Delta\nu$'],0.24,5,4)
#
# Everest vs KASOC - Yvonne
# k2p.plt_comp_scat(YY_C6,'Ynumax','nmx','e_Ynumax','nmx_err',[r'$(\nu_{\rm{max},EV}-\nu_{\rm{max},K2P2})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.3,4,4)
# k2p.plt_comp_scat(YY_C3,'Ynumax','nmx','e_Ynumax','nmx_err',[r'$(\nu_{\rm{max},EV}-\nu_{\rm{max},K2P2})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.32,8,7)
# k2p.plt_comp_scat(YY_C6,'YDnu','dnu','e_YDnu','dnu_err',[r'$(\Delta\nu_{EV}-\Delta\nu_{K2P2})/\sigma_{comp}$',r'$\Delta\nu$'],0.53,4,3)
# k2p.plt_comp_scat(YY_C3,'YDnu','dnu','e_YDnu','dnu_err',[r'$(\Delta\nu_{EV}-\Delta\nu_{K2P2})/\sigma_{comp}$',r'$\Delta\nu$'],0.35,7,6)
#
# Benoit vs Yvonne EVEREST
# k2p.plt_comp_scat(EE_C6,'Enumax','nmx','e_Enumax','nmx_err',[r'$(\nu_{\rm{max},EVEl}-\nu_{\rm{max},EVMo})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.30,5,4)
# k2p.plt_comp_scat(EE_C3,'Enumax','nmx','e_Enumax','nmx_err',[r'$(\nu_{\rm{max},EVEl}-\nu_{\rm{max},EVMo})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.35,6,5)
# k2p.plt_comp_scat(EE_C6,'EDnu','dnu','e_EDnu','dnu_err',[r'$(\Delta\nu_{EVEl}-\Delta\nu_{EVMo})/\sigma_{comp}$',r'$\Delta\nu$'],0.35,5,3)
# k2p.plt_comp_scat(EE_C3,'EDnu','dnu','e_EDnu','dnu_err',[r'$(\Delta\nu_{EVEl}-\Delta\nu_{EVMo})/\sigma_{comp}$',r'$\Delta\nu$'],0.37,5,4)
#
# # Benoit vs Yvonne KASOC
# k2p.plt_comp_scat(YB_C6,'Bnumax','Ynumax','e_Bnumax','e_Ynumax',[r'$(\nu_{\rm{max},El}-\nu_{\rm{max},Mo})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.34,5.5,4)
# k2p.plt_comp_scat(YB_C3,'Bnumax','Ynumax','e_Bnumax','e_Ynumax',[r'$(\nu_{\rm{max},El}-\nu_{\rm{max},Mo})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.33,5,4)
# k2p.plt_comp_scat(YB_C6,'BDnu','YDnu','e_BDnu','e_YDnu',[r'$(\Delta\nu_{El}-\Delta\nu_{Mo}/\sigma_{comp}$',r'$\Delta\nu$'],0.3,5,5)
# k2p.plt_comp_scat(YB_C3,'BDnu','YDnu','e_BDnu','e_YDnu',[r'$(\Delta\nu_{El}-\Delta\nu_{Mo}/\sigma_{comp}$',r'$\Delta\nu$'],0.33,5,4)
#
# # Benoit vs Savita KASOC
# k2p.plt_comp_scat(BS_C6,'Bnumax','Snumax','e_Bnumax','e_Snumax',[r'$(\nu_{\rm{max},Ma}-\nu_{\rm{max},Mo})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.7,5,4)
# k2p.plt_comp_scat(BS_C3,'Bnumax','Snumax','e_Bnumax','e_Snumax',[r'$(\nu_{\rm{max},Ma}-\nu_{\rm{max},Mo})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.6,5,4)
# k2p.plt_comp_scat(BS_C6,'BDnu','SDnu','e_BDnu','e_SDnu',[r'$(\Delta\nu_{Ma}-\Delta\nu_{Mo})/\sigma_{comp}$',r'$\Delta\nu$'],0.18,6,5)
# k2p.plt_comp_scat(BS_C3,'BDnu','SDnu','e_BDnu','e_SDnu',[r'$(\Delta\nu_{Ma}-\Delta\nu_{Mo})/\sigma_{comp}$',r'$\Delta\nu$'],0.18,6,5)
#
# # Yvonne vs Savita KASOC
# k2p.plt_comp_scat(YS_C6,'Ynumax','Snumax','e_Ynumax','e_Snumax',[r'$(\nu_{\rm{max},Ma}-\nu_{\rm{max},El})/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.42,4,4)
# k2p.plt_comp_scat(YS_C3,'Ynumax','Snumax','e_Ynumax','e_Snumax',[r'$(\nu_{\rm{max},Ma}-\nu_{\rm{max},E}l)/\sigma_{comp}$',r'$\nu_{\rm{max}}$'],0.58,8,7)
# k2p.plt_comp_scat(YS_C6,'YDnu','SDnu','e_YDnu','e_SDnu',[r'$(\Delta\nu_{Ma}-\Delta\nu_{El})/\sigma_{comp}$',r'$\Delta\nu$'],0.35,7,5)
# k2p.plt_comp_scat(YS_C3,'YDnu','SDnu','e_YDnu','e_SDnu',[r'$(\Delta\nu_{Ma}-\Delta\nu_{El})/\sigma_{comp}$',r'$\Delta\nu$'],0.35,8,7)

''' HRDs of simulated and real data '''
# plt.figure()
# k2p.plot_simHRD(YC3,besa3,lab=[r'Besancon C3',r'C3 Seismic'],b=0)
# k2p.plot_HRD(SC3,lab=[],b=0)
# k2p.plot_HRD(BC3,lab=[],b=0)
#
# plt.figure()
# k2p.plot_simHRD(YC6,besa6,lab=[r'Besancon 6',r'C6 Seismic'],b=0)
# k2p.plot_HRD(SC6,lab=[],b=0)
# k2p.plot_HRD(BC6,lab=[],b=0)
# k2p.plot_HRD(EC6,lab=[],b=0)

''' Mass vs Vertical Height '''
# print( YC3.columns.values)
# plt.figure()
# plt.scatter(besa3['IniMass'],besa3['zgal'],alpha=0.4)
# plt.scatter(YC3['mass'],(YC3['Distance'])*np.sin(YC3['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(BC3['mass'],(BC3['Distance'])*np.sin(BC3['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(SC3['mass'],(SC3['Distance'])*np.sin(SC3['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(YR3['mass'],(YR3['distance_1'])*np.sin(YR3['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(SR3['mass'],(SR3['distance_1'])*np.sin(SR3['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(BR3['mass'],(BR3['distance_1'])*np.sin(BR3['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
#
# plt.figure()
# plt.scatter(besa6['IniMass'],besa6['zgal'],alpha=0.4)
# plt.scatter(YC6['mass'],(YC6['Distance'])*np.sin(YC6['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(BC6['mass'],(BC6['Distance'])*np.sin(BC6['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(SC6['mass'],(SC6['Distance'])*np.sin(SC6['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(YR6['mass'],(YR6['distance_1'])*np.sin(YR6['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(SR6['mass'],(SR6['distance_1'])*np.sin(SR6['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(BR6['mass'],(BR6['distance_1'])*np.sin(BR6['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)
# plt.scatter(ER6['mass'],(ER6['distance_1'])*np.sin(ER6['Glat']*np.pi/180)*1e-3,color='r',alpha=0.4)

''' Histograms of numax in 10muHz bins, histos for figs 17-24 '''
# plt.figure()
# for i in range(10,160,10):
    # plt.subplot(3,5,(i/10))
# A3 = besa3[(besa3['numax'] > i) & (besa3['numax'] <= i+10)]
# B3 = BC3[(BC3['Bnumax'] > i) & (BC3['Bnumax'] <= i+10)]
# A6 = besa6[(besa6['numax'] > i) & (besa6['numax'] <= i+10)]
# B6 = BC6[(BC6['Bnumax'] > i) & (BC6['Bnumax'] <= i+10)]
# GAP3 = GAP3[ (GAP3['Hmag'] > 7) & (GAP3['Hmag'] < 12) & (GAP3['JK'] > 0.5) ]
# # print( len(GAP6), len(besa6))
# a = np.linspace(7,12,50)
# plt.figure()
# # plt.hist(GAP3['Hmag'],bins=a,alpha=0.5,label=r'GAP')#,normed=True)
# # plt.hist(YC3['Hmag'],bins=a,color='r')#,normed=True)#,histtype='step')
# # plt.hist(SC3['Hmag'],bins=a,color='g')#,normed=True)#,histtype='step')
# plt.hist(besa3['Hmag'],bins=a,alpha=0.5,label=r'Sim')#,normed=True)#,histtype='step')
# plt.hist(BC3['Hmag'],bins=a,color='y',alpha=0.5,label=r'Mosser')#,normed=True)#,histtype='step')
# plt.xlabel(r'H')# ($%s < \nu_{\rm{max}} < %s \nu\rm{Hz}$)' %(i,i+10))
# plt.xlim(7,12)
# plt.legend(loc=2)
# #
# # GAP6 = GAP6[ (GAP6['Vcut'] > 9) & (GAP6['Vcut'] < 15) & (GAP6['JK'] > 0.5) ]
# a = np.linspace(8.5,15,50)
# plt.figure()
# # plt.hist(GAP6['Vcut'],bins=a,alpha=0.5,label=r'GAP')
# # plt.hist(YC6['Vcut'],bins=a,color='r')#,normed=True)#,histtype='step')
# # plt.hist(SC6['Vcut'],bins=a,color='g')#,normed=True)#,histtype='step')
# plt.hist(besa6['Vcut'],bins=a,alpha=0.5,label=r'Sim')#,histtype='step')
# plt.hist(BC6['Vcut'],bins=a,color='y',alpha=0.5,label=r'Mosser')#,normed=True)#,histtype='step')
# plt.xlabel(r'V')# ($%s < \nu_{\rm{max}} < %s \nu\rm{Hz}$)' %(i,i+10))
# plt.xlim(9,15)
#     # if i < 1:
# plt.legend(loc=2)
# #
# print( stats.ttest_ind(BC3['Hmag'],besa3['Hmag'],equal_var=False))
# print( stats.ttest_ind(BC3['Hmag'],besa3['Hmag'],equal_var=False))

''' High/low numax V-band distribution, KDE '''
# plt.figure()
# a = [besa3,YC3,SC3,BC3]#,EC6]
# nmax = ['numax','Ynumax','Snumax','Bnumax']#,'Enumax']
# colours = ['m','g', 'b', 'r','k']
# lab = [r'besa6 high $\nu_{\rm{max}}$',r'YC6',r'SC6',r'BC6',r'EC6']
# lab2 = [r'besa6 low $\nu_{\rm{max}}$',r'YC6',r'SC6',r'BC6',r'EC6']
# i=0
# while i < len(a):
#     data = a[i]
#     x = min(data[nmax[i]]) + (max(data[nmax[i]])-min(data[nmax[i]]))/2.0
#     b1 = data[data[nmax[i]] >= x]
#     k2p.plot_KDE(b1,param='Vmag',label=lab[i],colour=colours[i],linestyle='-')
#     i+=1
# i=0
# while i < len(a):
#     data = a[i]
#     x = min(data[nmax[i]]) + (max(data[nmax[i]])-min(data[nmax[i]]))/2.0
#     b2 = data[data[nmax[i]] < x]
#     k2p.plot_KDE(b2,param='Vmag',label=lab2[i],colour=colours[i],linestyle='--')
#     i+=1
# plt.xlabel(r'V',fontsize=20)
# plt.legend(prop={'size':15},loc=2)
# plt.tight_layout()

''' numax vs V-mag scatter, det. prob. colourbar '''
# sim = besa6[besa6['prob_s'] > 0.95]
# det = EC6[EC6['prob_s'] > 0.95]
# ndet = EC6[EC6['prob_s'] <= 0.95]
#
# plt.figure()
# plt.scatter(det['Enumax'],det['Vmag'],c='r',label=r'EC6 det')
# plt.scatter(ndet['Enumax'],ndet['Vmag'],c='lime',label=r'EC6 ndet')
# plt.scatter(sim['numax'],sim['Vmag'],c=sim['prob_s'],alpha=0.25,cmap=colormaps.parula,label=r'Besan')
# cbar = plt.colorbar()
# cbar.set_label(r'Detection Probability', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.xlabel(r'$\nu_{\rm{max}}$',fontsize=20)
# plt.ylabel(r'V',fontsize=20)
# plt.ylim(max(besa3['Vmag']),min(besa3['Vmag']))
# plt.legend()

''' 2D hists for completeness comps '''
# TRILEGAL_C6 = TRILEGAL_C6[TRILEGAL_C6.Radius < 8.0]
# camp6 = camp6[camp6.radius < 8.0]
#
# plt.figure()
# plt.subplot(2,2,1)
# hist1, xb1, yb1, im1 = plt.hist2d(BC3['JK'],BC3['Hmag'],bins=49,cmap=colormaps.parula)#,normed=True)
# cbar = plt.colorbar()
# cbar.set_label(r'Number', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'H',fontsize=20, labelpad=20)
# plt.ylim(max(BC3['Hmag']),min(BC3['Hmag']))
# plt.xlim(min(BC3['JK']),max(BC3['JK']))
# plt.xlabel(r'$\nu_{\rm{max}}$',fontsize=20, labelpad=10)
# plt.title(r'Benoit C3',fontsize=20)
# plt.tick_params(labelsize=15)
# # plt.tight_layout()
#
# plt.subplot(2,2,2)
# hist2, xb2, yb2, im2 = plt.hist2d(camp3['JK'],camp3['Hmag'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
# cbar = plt.colorbar()
# cbar.set_label(r'Number', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'H',fontsize=20, labelpad=20)
# plt.ylim(max(yb1),min(yb1))
# plt.xlim(min(xb1),max(xb1))
# # plt.ylim(max(TRILEGAL_C6['Vmag']),min(TRILEGAL_C6['Vmag']))
# # plt.xlim(min(TRILEGAL_C6['JK']),max(TRILEGAL_C6['JK']))
# plt.xlabel(r'J-K',fontsize=20, labelpad=10)
# plt.title(r'C3 Multi det',fontsize=20)
# plt.tick_params(labelsize=15)
# # plt.tight_layout()
#
# hist = hist2-hist1
# # hist = np.zeros((49,49))
# # print(hist.item(0))
# # h = 49**2
# # hist = [[(1 - hist2[i][j]/hist1[i][j])*100 for j in np.arange(49)] for i in np.arange(49)]
# # print(np.divide(hist1,hist2))
# # np.savetxt('/home/bmr165/GA/K2Poles/hist.txt',hist)
# # plt.figure()
# plt.subplot(2,2,3)
# plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
# cbar = plt.colorbar()
# cbar.set_label(r'$N_{\rm{obs}} - N_{\rm{sim}}$', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'H',fontsize=20, labelpad=20)
# plt.xlabel(r'J-K',fontsize=20, labelpad=10)
# plt.title(r'Multi det - Benoit (C3)',fontsize=20)
# plt.tick_params(labelsize=15)
#
# plt.tight_layout()
#
# # a = np.sum(hist2)
# # b = np.sum(hist1)
# # c = np.float(b/a)*100
# # print( a, b, c)
#
# C3_over = pd.concat([C3_1,GAP3_1],ignore_index=True)
# C3_over = C3_over.drop_duplicates(subset=['EPIC'],keep=False)
# C3_over = C3_over.reset_index(drop=True)
# C6_over = pd.concat([C6_1,GAP6_1],ignore_index=True)
# C6_over = C6_over.drop_duplicates(subset=['EPIC'],keep=False)
# C6_over = C6_over.reset_index(drop=True)
#
# plt.figure()
# plt.scatter(C6['Glon'],C6['Glat'],alpha=0.25,label=r'C6 Field')
# plt.scatter(GAP6['Glon'],GAP6['Glat'],alpha=0.25,color='g',label=r'GAP6 Field')
# plt.scatter(C6_over['Glon'],C6_over['Glat'],color='r',label=r'Non-GAP')
# plt.xlabel(r'l')
# plt.ylabel(r'b')
# plt.legend()
# plt.tight_layout()

''' Sky projection plot '''
plt.figure()
# # plt.hist(C6_nospec['Grad']/1000,bins=75)
# v9 = C6_nospec[C6_nospec['Vmag'] < 10]
# v10 = C6_nospec[(C6_nospec['Vmag'] >= 10) & (C6_nospec['Vmag'] < 11)]
# v11 = C6_nospec[(C6_nospec['Vmag'] >= 11) & (C6_nospec['Vmag'] < 12)]
# v12 = C6_nospec[(C6_nospec['Vmag'] >= 12) & (C6_nospec['Vmag'] < 13)]
# v13 = C6_nospec[(C6_nospec['Vmag'] >= 13) & (C6_nospec['Vmag'] < 14)]
# v14 = C6_nospec[(C6_nospec['Vmag'] >= 14) & (C6_nospec['Vmag'] <= 15)]
# # plt.scatter(C6_nospec['RA'],C6_nospec['Dec'])
# plt.scatter(v9['RA'],v9['Dec'],label=r'V = 9-10')
# plt.scatter(v10['RA'],v10['Dec'],label=r'V = 10-11',color='r')
# plt.scatter(v11['RA'],v11['Dec'],label=r'V = 11-12',color='g')
# plt.scatter(v12['RA'],v12['Dec'],label=r'V = 12-13',color='m')
# plt.scatter(v13['RA'],v13['Dec'],label=r'V = 13-14',color='c')
# plt.scatter(v14['RA'],v14['Dec'],label=r'V = 14-15',color='k')
plt.scatter(TRILEGAL_C3['RA'],TRILEGAL_C3['Dec'],label=r'TRILEGAL C3',alpha=0.5)
plt.scatter(camp3_0['RA'],camp3_0['Dec'],label=r'K2 C3',alpha=0.5)
plt.xlabel(r'RA')
plt.ylabel(r'DEC')
plt.legend()

plt.show()
