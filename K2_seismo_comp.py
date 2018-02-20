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

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

''' Dropbox Path '''
ext_DB = '/home/bmr135/' # Work
# ext_DB = '/home/ben/'   # Laptop
''' GA directory '''
ext_GA = '/home/bmr135/' # Work
# ext_GA = '/media/ben/SAMSUNG/' # Hard-Drive

''' Read in simulated and real data '''
besa3, besa6 = dat.BESANCON()
print("Besancon")
TRILEGAL_C3, TRILEGAL_C6 = dat.TRILEGAL()
print("TRILEGAL")
C3 = dat.C3_cat()
print( "C3")
C6 = dat.C6_cat()
print("C6")
GAP3, GAP6 = dat.K2_GAP()
print( "GAP")
Yvonne_C3, Yvonne_C6, Yvonne_EC6, Yvonne_EC3 = dat.Yvonne()
print( "Yvonne")
Savita_C3, Savita_C6, Savita_EC3, Savita_EC6 = dat.Savita()
print( "Savita")
Benoit_C3, Benoit_C6, Everest_C3, Everest_C6 = dat.Benoit()
print( "Benoit")
RAVE3, RAVE6 = dat.RAVE()
print( "RAVE")
GES3 = dat.Gaia_ESO()
print( "Gaia-ESO")
APO6 = dat.APOGEE()
print( "APOGEE ")
LAMOST3, LAMOST6 = dat.LAMOST()
print("LAMOST")
oc = dat.occurrence()
print( "Occurrence in")

''' Add K2P2 occurrence to GAP Target Lists '''
GAP6 = dat.n_epics(GAP6,oc)

''' Preparing GAP for propbability detections '''
GAP3['foma'] = GAP3['mass'] * GAP3['Radius']**-2 * (GAP3['Teff']/5777)**-0.5 * 3090 # numax for C3 GAP (frequency of maximum amplitude)
GAP3['Lumo'] = GAP3['Radius']**2 * (GAP3['Teff']/const.solar_Teff)**4
GAP3_v2 = GAP3[GAP3['foma'] < 280]
GAP3_v2 = GAP3_v2[GAP3_v2['foma'] > 10]
GAP3_v2 = GAP3_v2[GAP3_v2['imag'] > 0.0]
GAP3_v2 = GAP3_v2.reset_index(drop=True)
GAP3_v2 = prop.det_prob_GAP(GAP3_v2,'foma',3090,135.1)
GAP3_v2 = GAP3_v2[GAP3_v2['prob_s'] >= 0.95]
GAP3_v2.to_csv('/home/bmr135/GA/K2Poles/GAP3')

GAP6['foma'] = GAP6['mass'] * GAP6['Radius']**-2 * (GAP6['Teff']/5777)**-0.5 * 3090 # numax for C6 GAP (frequency of maximum amplitude)
GAP6['Lumo'] = GAP6['Radius']**2 * (GAP6['Teff']/const.solar_Teff)**4
GAP6_v2 = GAP6[GAP6['foma'] < 280]
GAP6_v2 = GAP6_v2[GAP6_v2['foma'] > 10]
GAP6_v2 = GAP6_v2[GAP6_v2['imag'] > 0.0]
GAP6_v2 = GAP6_v2.reset_index(drop=True)
GAP6_v2 = prop.det_prob_GAP(GAP6_v2,'foma',3090,135.1)
GAP6_v2 = GAP6_v2[GAP6_v2['prob_s'] >= 0.95]
GAP6_v2.to_csv('/home/bmr135/GA/K2Poles/GAP6')

''' Merge data with GAP target lists '''
YC3 = pd.merge(Yvonne_C3,GAP3,how='inner',on=['EPIC'])
YC6 = pd.merge(Yvonne_C6,GAP6,how='inner',on=['EPIC'])
SC3 = pd.merge(Savita_C3,GAP3,how='inner',on=['EPIC'])
SC6 = pd.merge(Savita_C6,GAP6,how='inner',on=['EPIC'])
BC3 = pd.merge(Benoit_C3,GAP3,how='inner',on=['EPIC'])
BC6 = pd.merge(Benoit_C6,GAP6,how='inner',on=['EPIC'])
EC6 = pd.merge(Everest_C6,GAP6,how='inner',on=['EPIC'])
YEC3 = pd.merge(Yvonne_EC3,GAP3,how='inner',on=['EPIC'])
YEC6 = pd.merge(Yvonne_EC6,GAP6,how='inner',on=['EPIC'])
SEC3 = pd.merge(Savita_EC3,GAP3,how='inner',on=['EPIC'])
SEC6 = pd.merge(Savita_EC6,GAP6,how='inner',on=['EPIC'])
EC3 = pd.merge(Everest_C3,GAP3,how='inner',on=['EPIC'])
GG3 = pd.merge(GES3,GAP3,how='inner',on=['EPIC'])
AC6 = pd.merge(APO6,GAP6,how='inner',on=['EPIC'])
# LC3 = pd.merge(LAMOST3,GAP3,how='inner',on=['EPIC'])
# LC6 = pd.merge(LAMOST6,GAP6,how='inner',on=['EPIC'])

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
sel_numax = ['nmx','Snumax','Bnumax','Enumax','nmx','Snumax','numax','nmx','Snumax','Bnumax','Enumax','numax','nmx','Snumax']
sel_list = [YC3,SC3,BC3,EC3,YEC3,SEC3,besa3,YC6,SC6,BC6,EC6,besa6,YEC6,SEC6]

YC3,YC6,SC3,SC6,BC3,BC6,EC6,YEC6,EC3,YEC3,SEC3,SEC6 = prop.galactic_coords(seismo_list)
C3,C6,GAP3,GAP6 = prop.galactic_coords([C3,C6,GAP3,GAP6])
YC3,YC6,SC3,SC6,BC3,BC6,EC6,YEC6,EC3,YEC3,SEC3,SEC6 = prop.lmrl_comps(seismo_list,numax,dnu,Numax,Dnu,1)

plt.figure()
GAP6 = GAP6[GAP6['[Fe/H]'] > -10]
plt.hist(GAP6['[Fe/H]'],bins=50)
plt.show()

def hist_orig(df,df1,cut,bins,ext,n):
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

    hist0[0], d[0], c = ax0.hist(df['Mass'],bins=bins[0],histtype='step')
    hist[0], b[0], c = ax0.hist(df1['Mass'],bins=bins[0])
    ax0.set_xlabel(r'Mass [M$_{\odot}$]')
    if n == 6:
        ax0.set_xlim(0.5,2.5)

    hist0[1], d[1], c = ax1.hist(df['radius'],bins=bins[1],histtype='step')
    hist[1], b[1], c = ax1.hist(df1['radius'],bins=bins[1])
    ax1.set_xlabel(r'Radius [R$_{\odot}$]')
    if n == 6:
        ax1.set_xlim(0,30)

    hist0[2], d[2], c = ax2.hist(df['Teff'],bins=bins[2],histtype='step')
    hist[2], b[2], c = ax2.hist(df1['Teff'],bins=bins[2])
    ax2.set_xlabel(r'T$_{\rm{eff}}$ [K]')
    if n == 6:
        ax2.set_xlim(4000,5250)

    hist0[3], d[3], c = ax3.hist(df['M_H'],bins=bins[3],histtype='step')
    hist[3], b[3], c = ax3.hist(df1['M_H'],bins=bins[3])
    ax3.set_xlabel(r'[Fe/H] [dex]')
    if n == 6:
        ax3.set_xlim(-1.2,0.2)

    hist0[4], d[4], c = ax4.hist(df['logg'],bins=bins[4],histtype='step')
    hist[4], b[4], c = ax4.hist(df1['logg'],bins=bins[4])
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
    plt.savefig(ext+'_'+'TC6.png')

    if n == 6:
        fig, axes = plt.subplots(3,2)
        ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()

        ax0.scatter(bins[0][:-1],ratio[0])
        ax0.set_xlabel(r'Mass [M$_{\odot}$]')
        ax0.set_xlim(0.5,2.5)

        ax1.scatter(bins[1][:-1],ratio[1])
        ax1.set_xlabel(r'Radius [R$_{\odot}$]')
        ax1.set_xlim(0,30)

        ax2.scatter(bins[2][:-1],ratio[2])
        ax2.set_xlabel(r'T$_{\rm{eff}}$ [K]')
        ax2.set_xlim(4000,5250)

        ax3.scatter(bins[3][:-1],ratio[3])
        ax3.set_xlabel(r'[Fe/H] [dex]')
        ax3.set_xlim(-1.2,0.2)

        ax4.scatter(bins[4][:-1],ratio[4])
        ax4.set_xlabel(r'log$_{10}$(g)')
        ax4.set_xlim(1,4)

        ax5.axis('off')
        plt.tight_layout()
        plt.show()

    return hist, b

''' TRILEGAL selection cuts '''
TRILEGAL_C3 = prop.det_prob(TRILEGAL_C3,'numax',3090.0,135.1)
TRILEGAL_C6 = prop.det_prob(TRILEGAL_C6,'numax',3090.0,135.1)

TRILEGAL_C3 = TRILEGAL_C3[ (TRILEGAL_C3['numax'] > 10) & (TRILEGAL_C3['numax'] < 280) & \
(TRILEGAL_C3['Hmag'] > 7) & (TRILEGAL_C3['Hmag'] < 12) & (TRILEGAL_C3['JK'] > 0.5) & \
(TRILEGAL_C3['prob_s'] > 0.95)]
TRILEGAL_C6 = TRILEGAL_C6[ (TRILEGAL_C6['numax'] > 10) & (TRILEGAL_C6['numax'] < 280) & \
(TRILEGAL_C6['Vcut'] > 9) & (TRILEGAL_C6['Vcut'] < 15) & (TRILEGAL_C6['JK'] > 0.5) & \
(TRILEGAL_C6['prob_s'] > 0.95)]

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

TRI3 = TRILEGAL_C3[TRILEGAL_C3['logAge'] > np.log10(9e9)]
TRI6 = TRILEGAL_C6[TRILEGAL_C6['logAge'] > np.log10(9e9)]
TRI3.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/TRILEGAL_C3_old_stars',index=False)
TRI6.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/TRILEGAL_C6_old_stars',index=False)
print('Trilegal saved out')

YC3,SC3,BC3,EC3,YEC3,SEC3,besa3,YC6,SC6,BC6,EC6,besa6,YEC6,SEC6 = prop.selection_function(sel_list,sel_numax)

Y6_dnu = pd.DataFrame()
Y6_dnu = YC6[YC6['dnu'] != np.nan]
Y6_dnu.reset_index(drop=True)
print(len(Y6_dnu), len(YC6))

plt.figure()
hist, bins, patches = plt.hist(YC6['Radius'],histtype='step',bins=50,label=r'All Data')#,normed=True)
plt.hist(Y6_dnu['Radius'],bins=bins,histtype='step',label=r'With $\Delta\nu$')#,normed=True)
plt.xlabel(r'Radius [R$_{\odot}$]')
plt.legend()
plt.tight_layout(rect=[0, 0.03, 1, 0.95])


''' Add detection flags to data/save out values for comparisons '''
# YC3,BC3,SC3 = prop.individ(YC3,BC3,SC3,'K2P2_C3')
# YC6,BC6,SC6 = prop.individ(YC6,BC6,SC6,'K2P2_C6')
# YEC3,EC3,SEC3 = prop.individ(YEC3,EC3,SEC3,'EVEREST_C3')
# YEC6,EC6,SEC6 = prop.individ(YEC6,EC6,SEC6,'EVEREST_C6')

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
GES.to_csv(ext_GA+'GA/K2Poles/Gaia_ESO/GES_full.csv',index=False,na_rep='Inf')
print( "Gaia-ESO saved out")

''' Merging of APOGEE data with single asteroseismic dets '''
YA6,SA6,BA6,EA6 = dat.APO_merge(seismo6_list,APO6,seismo6_name)
APO = pd.concat([BA6,SA6,YA6],ignore_index=True)
APO = APO.drop_duplicates(subset=['EPIC'])
APO = APO.fillna(value='NaN',method=None)
APO = APO.reset_index(drop=True)
APO.to_csv(ext_GA+'GA/K2Poles/APO_LAMOST/APOGEE_full.csv',index=False,na_rep='Inf')
print( "APOGEE saved out")

''' Merging of LAMOST data with single asteroseismic dets '''
YL3 = pd.merge(YC3,LAMOST3,how='inner',on=['EPIC'])
BL3 = pd.merge(BC3,LAMOST3,how='inner',on=['EPIC'])
SL3 = pd.merge(SC3,LAMOST3,how='inner',on=['EPIC'])
LST = pd.concat([BL3,SL3,YL3],ignore_index=True)
LST = LST.drop_duplicates(subset=['EPIC'])
LST = LST.fillna(value='NaN',method=None)
LST = LST.reset_index(drop=True)
LST.to_csv(ext_GA+'GA/K2Poles/APO_LAMOST/LAMOST_full_C3.csv',index=False,na_rep='Inf')
print( "LAMOST C3 saved out ", len(LST))

YL6 = pd.merge(YC6,LAMOST6,how='inner',on=['EPIC'])
BL6 = pd.merge(BC6,LAMOST6,how='inner',on=['EPIC'])
SL6 = pd.merge(SC6,LAMOST6,how='inner',on=['EPIC'])
LAMOST = pd.concat([BL6,SL6,YL6],ignore_index=True)
LAMOST = LAMOST.drop_duplicates(subset=['EPIC'],keep='first')
LAMOST = LAMOST.fillna(value='NaN',method=None)
LAMOST = LAMOST.reset_index(drop=True)
LAMOST.to_csv(ext_GA+'GA/K2Poles/APO_LAMOST/LAMOST_full_C6.csv',index=False,na_rep='Inf')
print( "LAMOST C6 saved out ", len(LAMOST))


''' Complete asteroseismic lists '''
camp3_0 = pd.concat([YC3,SC3,BC3],ignore_index=True)
camp3_0 = camp3_0.drop_duplicates(subset=['EPIC'])
camp3_0 = camp3_0.reset_index(drop=True)
camp3_0 = camp3_0.fillna(value='NaN',method=None)

camp6_0 = pd.concat([YC6,SC6,BC6],ignore_index=True)
camp6_0 = camp6_0.drop_duplicates(subset=['EPIC'])
camp6_0 = camp6_0.reset_index(drop=True)
camp6_0 = camp6_0.fillna(value='NaN',method=None)

print(len(camp3_0),len(camp6_0))


''' Merging of multiple fields for comparison '''
cols_to_use = YC3.columns.difference(SC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YS_C3 = pd.merge(SC3,YC3[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = YC6.columns.difference(SC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YS_C6 = pd.merge(SC6,YC6[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = YC3.columns.difference(BC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YB_C3 = pd.merge(BC3,YC3[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = YC6.columns.difference(BC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
YB_C6 = pd.merge(BC6,YC6[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = SC3.columns.difference(BC3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
BS_C3 = pd.merge(BC3,SC3[cols_to_use],how='inner',on=['EPIC'])
cols_to_use = SC6.columns.difference(BC6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
BS_C6 = pd.merge(BC6,SC6[cols_to_use],how='inner',on=['EPIC'])
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

for i in range(0,len(merged),1):
    merged[i] = prop.sigma_clip(merged[i],nmx1[i],nmx2[i],dnu1[i],dnu2[i], \
                err_nmx1[i],err_nmx2[i],err_dnu1[i],err_dnu2[i],3)
''' Merged sigma clipped parameters. Turn on for 1-to-1 plots, saving out and
other general use. Toggle off for sigma comparative plots unless using a smaller
sigma subset to calculate the std. dev. '''
YS_C3 = merged[0][0]
Z0 = merged[0][1]
YS_C6 = merged[1][0]
Z1 = merged[1][1]
YB_C3 = merged[2][0]
Z2 = merged[2][1]
YB_C6 = merged[3][0]
Z3 = merged[3][1]
BS_C3 = merged[4][0]
Z4 = merged[4][1]
BS_C6 = merged[5][0]
Z5 = merged[5][1]
EB_C3 = merged[6][0]
Z6 = merged[6][1]
EB_C6 = merged[7][0]
Z7 = merged[7][1]
EE_C3 = merged[8][0]
Z8 = merged[8][1]
EE_C6 = merged[9][0]
Z9 = merged[9][1]
YY_C3 = merged[10][0]
Z10 = merged[10][1]
YY_C6 = merged[11][0]
Z11 = merged[11][1]
SS_C3 = merged[12][0]
Z12 = merged[12][1]
SS_C6 = merged[13][0]
Z13 = merged[13][1]
YSE_C3 = merged[14][0]
Z14 = merged[14][1]
YSE_C6 = merged[15][0]
Z15 = merged[15][1]
SE_C3 = merged[16][0]
Z16 = merged[16][1]
SE_C6 = merged[17][0]
Z17 = merged[17][1]

print( "Sigma clip complete")

''' Merging RAVE data with asteroseismic '''
RAVE_list = [RAVE3,RAVE6,RAVE3,RAVE6,RAVE3,RAVE6,RAVE6,RAVE6,RAVE3,RAVE3]
YSR3,YSR6,BSR3,BSR6,YBR3,YBR6 = dat.RAVE_merge([YS_C3,YS_C6,BS_C3,BS_C6,YB_C3,YB_C6], \
                                RAVE_list,['YSR3','YSR6','BSR3','BSR6','YBR3','YBR6'])

RAVE3 = pd.concat([YBR3,YSR3,BSR3],ignore_index=True)
RAVE3 = RAVE3.drop_duplicates(subset=['EPIC'])
RAVE3 = RAVE3.fillna(value='NaN',method=None)
RAVE3 = RAVE3.reset_index(drop=True)
RAVE3.to_csv(ext_GA+'GA/K2Poles/RAVE_C3.csv',index=False,na_rep='Inf')
RAVE6 = pd.concat([YBR6,YSR6,BSR6],ignore_index=True)
RAVE6 = RAVE6.drop_duplicates(subset=['EPIC'])
RAVE6 = RAVE6.fillna(value='NaN',method=None)
RAVE6 = RAVE6.reset_index(drop=True)
RAVE6.to_csv(ext_GA+'GA/K2Poles/RAVE_C6.csv',index=False,na_rep='Inf')
# GR = pd.merge(RAVE3,GES,how='inner',on=['EPIC'])
# RAVEGES = pd.concat([RAVE3,GES],ignore_index=True)
# RAVEGES = RAVEGES.drop_duplicates(subset=['EPIC'])
# print( len(RAVEGES))

APORAVE = pd.merge(APO,RAVE6,how='inner',on=['EPIC'])
# print( len(APORAVE))

''' Concatenated lists of stars in each campaign, K2P2 and Everest.
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
camp3.to_csv(ext_GA+'GA/K2Poles/matlab_in/C3_070218.csv',index=False)

camp6 = pd.concat([YB_C6,YS_C6,BS_C6],ignore_index=True)
camp6 = camp6.drop_duplicates(subset=['EPIC'])
camp6 = camp6.reset_index(drop=True)
camp6 = camp6.fillna(value='NaN',method=None)

camp6 = prop.single_seismo(camp6,['Bnumax','nmx','Snumax'],'NUMAX')
camp6 = prop.single_seismo(camp6,['e_Bnumax','nmx_err','e_Snumax'],'NUMAX_err')
camp6 = prop.single_seismo(camp6,['BDnu','dnu','SDnu'],'DNU')
camp6 = prop.single_seismo(camp6,['e_BDnu','dnu_err','e_SDnu'],'DNU_err')
# camp6 = prop.met_filter(camp6)
camp6.to_csv(ext_GA+'GA/K2Poles/matlab_in/C6_070218.csv',index=False)

print(len(camp3),len(camp6))

''' Data Flag for EPIC parametric values - which [Fe/H] to use? '''
spectro_EPICS_3 = camp3[camp3['stpropflag'] != 'rpm']
spectro_EPICS_6 = camp6[camp6['stpropflag'] != 'rpm']
spectro_EPICS_3 = spectro_EPICS_3[spectro_EPICS_3['stpropflag'] != 'rav']
spectro_EPICS_6 = spectro_EPICS_6[spectro_EPICS_6['stpropflag'] != 'rav']
print(len(spectro_EPICS_6['stpropflag']))

GAP_camp3 = pd.merge(GAP3_v2,camp3,how='inner',on=['EPIC'])
GAP_camp6 = pd.merge(GAP6_v2,camp6,how='inner',on=['EPIC'])
print(len(GAP_camp3),len(GAP_camp6))

fig,ax = plt.subplots(2)
ax[0].hist(GAP3_v2['Hmag'],bins=np.linspace(7,12,50),histtype='step',label=r'GAP')
ax[0].hist(camp3['Hmag'],bins=np.linspace(7,12,50),histtype='step',label=r'Seismic')
ax[0].hist(GAP_camp3['Hmag_x'],bins=np.linspace(7,12,50),histtype='step',label=r'Combined')
ax[0].legend()
ax[0].set_xlabel(r'H')
ax[1].hist(GAP6_v2['Vcut'],bins=np.linspace(9,15,50),histtype='step',label=r'GAP')
ax[1].hist(camp6['Vcut'],bins=np.linspace(9,15,50),histtype='step',label=r'Seismic')
ax[1].hist(GAP_camp6['Vcut_x'],bins=np.linspace(9,15,50),histtype='step',label=r'Combined')
ax[1].set_xlabel(r'V')
plt.tight_layout()

cols_to_use = camp3.columns.difference(RAVE3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
RC3 = pd.merge(RAVE3,camp3[cols_to_use],how='inner',on=['EPIC'])
RC3.to_csv(ext_GA+'GA/K2Poles/matlab_in/RC3_070218.csv',index=False,na_rep='Inf')
cols_to_use = camp6.columns.difference(RAVE6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
RC6 = pd.merge(RAVE6,camp6[cols_to_use],how='inner',on=['EPIC'])
RC6.to_csv(ext_GA+'GA/K2Poles/matlab_in/RC6_070218.csv',index=False,na_rep='Inf')
print("RAVE saved out", len(RC3), len(RC6))

''' Merging of GES data with multiple asteroseismic dets '''
cols_to_use = camp3.columns.difference(GES3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
GES = pd.merge(GES3,camp3[cols_to_use],how='inner',on=['EPIC'])
GES.to_csv(ext_GA+'GA/K2Poles/matlab_in/GES_070218.csv',index=False,na_rep='Inf')
print( "Gaia-ESO saved out", len(GES))

''' Merging of LAMOST data with multiple asteroseismic dets '''
cols_to_use = camp3.columns.difference(LAMOST3.columns)
cols_to_use = cols_to_use.union(['EPIC'])
L3 = pd.merge(LAMOST3,camp3[cols_to_use],how='inner',on=['EPIC'])
L3.to_csv(ext_GA+'GA/K2Poles/matlab_in/LAMOST3_070218.csv',index=False,na_rep='Inf')
print( "LAMOST3 saved out ", len(LAMOST3))
cols_to_use = camp6.columns.difference(LAMOST6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
L6 = pd.merge(LAMOST6,camp6[cols_to_use],how='inner',on=['EPIC'])
L6.to_csv(ext_GA+'GA/K2Poles/matlab_in/LAMOST6_070218.csv',index=False,na_rep='Inf')
print( "LAMOST6 saved out ", len(LAMOST6))

''' Merging of APOGEE data with multiple asteroseismic dets '''
cols_to_use = camp6.columns.difference(APO6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
APO = pd.merge(APO6,camp6[cols_to_use],how='inner',on=['EPIC'])
APO.to_csv(ext_GA+'GA/K2Poles/APOGEE_070218.csv',index=False,na_rep='Inf')
print( "APOGEE saved out", len(APO))

print(len(pd.merge(L3,GES,how='inner',on=['EPIC'])))
print(len(pd.merge(L3,RC3,how='inner',on=['EPIC'])))
print(len(pd.merge(L6,RC6,how='inner',on=['EPIC'])))
print(len(pd.merge(L6,APO,how='inner',on=['EPIC'])))
print(len(L6),len(RC6))
cols_to_use = RC6.columns.difference(L6.columns)
cols_to_use = cols_to_use.union(['EPIC'])
LR6 = pd.merge(L6,RC6[cols_to_use],how='inner',on=['EPIC'])
LR6.to_csv(ext_GA+'GA/K2Poles/LAMOST_RAVE_C6.csv',index=False,na_rep='Inf')
print('LAMOST/RAVE C6 saved out')

E3 = pd.concat([EE_C3,SE_C3,YSE_C3],ignore_index=True)
E3 = E3.drop_duplicates(subset=['EPIC'])
E3 = E3.reset_index(drop=True)
E_camp3 = pd.merge(camp3,EC3,how='inner',on=['EPIC'])

E6 = pd.concat([EE_C6,SE_C6,YSE_C6],ignore_index=True)
E6 = E6.drop_duplicates(subset=['EPIC'])
E6 = E6.reset_index(drop=True)
E_camp6 = pd.merge(camp6,EC6,how='inner',on=['EPIC'])

E_camp6 = dat.n_epics(E_camp6,oc)

K2PEV = n1 = n2 = n3 = x = pd.DataFrame()
x = camp6[camp6['occ'] == 2.0]
K2PEV = camp6[~camp6['EPIC'].isin(E_camp6['EPIC'])].dropna()
n1 = K2PEV[K2PEV['occ'] == 1.0]
n2 = K2PEV[K2PEV['occ'] == 2.0]
n3 = K2PEV[K2PEV['occ'] == 3.0]
print( len(K2PEV), len(n1), len(n2), len(n3))
print( 'Percentage Single = ', 100*len(n1)/len(K2PEV))
print( 'Percentage Double = ', 100*len(n2)/len(K2PEV))
print( 'Percentage Triple = ', 100*len(n3)/len(K2PEV))

print( "Processing Complete")


''' Check length of data sets '''
# print( len(BC3), len(YC3), len(SC3), len(BC6), len(YC6), len(SC6))
# print( len(BR3), len(YR3), len(SR3), len(BR6), len(YR6), len(SR6))
# print( len(YEC3), len(YC6), len(EE_C6))
# print( C3.columns.values)

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
plt.figure()
plt.subplot(2,2,1)
hist1, xb1, yb1, im1 = plt.hist2d(BC3['JK'],BC3['Hmag'],bins=49,cmap=colormaps.parula)#,normed=True)
cbar = plt.colorbar()
cbar.set_label(r'Number', rotation=270, fontsize=20, labelpad=25)
cbar.ax.tick_params(labelsize=20)
plt.ylabel(r'H',fontsize=20, labelpad=20)
plt.ylim(max(BC3['Hmag']),min(BC3['Hmag']))
plt.xlim(min(BC3['JK']),max(BC3['JK']))
plt.xlabel(r'$\nu_{\rm{max}}$',fontsize=20, labelpad=10)
plt.title(r'Benoit C3',fontsize=20)
plt.tick_params(labelsize=15)
# plt.tight_layout()

plt.subplot(2,2,2)
hist2, xb2, yb2, im2 = plt.hist2d(camp3['JK'],camp3['Hmag'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
cbar = plt.colorbar()
cbar.set_label(r'Number', rotation=270, fontsize=20, labelpad=25)
cbar.ax.tick_params(labelsize=20)
plt.ylabel(r'H',fontsize=20, labelpad=20)
plt.ylim(max(yb1),min(yb1))
plt.xlim(min(xb1),max(xb1))
# plt.ylim(max(TRILEGAL_C6['Vmag']),min(TRILEGAL_C6['Vmag']))
# plt.xlim(min(TRILEGAL_C6['JK']),max(TRILEGAL_C6['JK']))
plt.xlabel(r'J-K',fontsize=20, labelpad=10)
plt.title(r'C3 Multi det',fontsize=20)
plt.tick_params(labelsize=15)
# plt.tight_layout()

hist = hist2-hist1
# hist = np.zeros((49,49))
# print(hist.item(0))
# h = 49**2
# hist = [[(1 - hist2[i][j]/hist1[i][j])*100 for j in np.arange(49)] for i in np.arange(49)]
# print(np.divide(hist1,hist2))
# np.savetxt('/home/bmr165/GA/K2Poles/hist.txt',hist)
# plt.figure()
plt.subplot(2,2,3)
plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
cbar = plt.colorbar()
cbar.set_label(r'$N_{\rm{obs}} - N_{\rm{sim}}$', rotation=270, fontsize=20, labelpad=25)
cbar.ax.tick_params(labelsize=20)
plt.ylabel(r'H',fontsize=20, labelpad=20)
plt.xlabel(r'J-K',fontsize=20, labelpad=10)
plt.title(r'Multi det - Benoit (C3)',fontsize=20)
plt.tick_params(labelsize=15)

plt.tight_layout()

# a = np.sum(hist2)
# b = np.sum(hist1)
# c = np.float(b/a)*100
# print( a, b, c)

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
# plt.figure()
# plt.hist(BC3['Grad']/1000,bins=75)
# plt.xlabel(r'Galactocentric Distance [kpc]')


# plt.show()
