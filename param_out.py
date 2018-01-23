''' Reading in of PARAM output files for K2 Poles project '''

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.mlab as mlab
import pandas as pd
import sys
from pandas import DataFrame, read_csv
from astropy.io import fits
from astropy.table import Table
import numpy as np
from scipy.stats import norm
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde
from pyqt_fit import kde
import scipy.odr.odrpack as odrpack
import colormaps
import K2_plot as k2p
import K2_properties as prop
import K2_constants as const
import K2_data as dat
import PARAM as p
from matplotlib.ticker import NullFormatter
# import seaborn as sn

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

''' APOKASC in '''
APO = pd.read_csv('/home/bmr135/GA/K2Poles/APOKASC4BEN.txt')

''' Isotrack from which to extract isochrone in '''
# iso = pd.read_csv('/home/bmr135/GA/K2Poles/ptcri_ovh0.20.ovhe0.50.z0.01756.y0.26618.LOW',skiprows=7,names=['age','logg','logT','dnu','numax','idx1','idx2'],delimiter=r'\s+')
iso = pd.read_csv('/home/bmr135/GA/K2Poles/ptcri_ovh0.20.ovhe0.50.z0.03123.y0.27994.LOW',skiprows=7,names=['age','logL','logT','dnu','DP','idx1','idx2'],delimiter=r'\s+')

df = ['m06','m08','m080','m10','m11','m12','m120','m13','m14','m15','m155','m16','m165','m17','m175','m18','m20']
index = [99,94,94,94,92,91,91,92,90,99,95,93,95,95,88,92,87]
cumul = [99,193,287,381,473,564,655,747,837,936,1031,1124,1219,1314,1402,1494,1581]

for i in range(0,len(df),1):
    df[i] = pd.DataFrame(np.zeros((index[i],7)),columns=['age','logL','logT','dnu','DP','idx1','idx2'])

# print np.shape(df[0])
# print iso.iloc[0]
df[0] = iso.iloc[0:index[0]]
df[1] = iso.iloc[cumul[0]:cumul[1]]
df[2] = iso.iloc[cumul[1]:cumul[2]]
df[3] = iso.iloc[cumul[2]:cumul[3]]
df[4] = iso.iloc[cumul[3]:cumul[4]]
df[5] = iso.iloc[cumul[4]:cumul[5]]
df[6] = iso.iloc[cumul[5]:cumul[6]]
df[7] = iso.iloc[cumul[6]:cumul[7]]
df[8] = iso.iloc[cumul[7]:cumul[8]]
df[9] = iso.iloc[cumul[8]:cumul[9]]
df[10] = iso.iloc[cumul[9]:cumul[10]]
df[11] = iso.iloc[cumul[10]:cumul[11]]
df[12] = iso.iloc[cumul[11]:cumul[12]]
df[13] = iso.iloc[cumul[12]:cumul[13]]
df[14] = iso.iloc[cumul[13]:cumul[14]]
df[15] = iso.iloc[cumul[14]:cumul[15]]
df[16] = iso.iloc[cumul[15]:cumul[16]]

#
# for i in range(0,len(iso),1):
#     if i < cumul[0]:
#         df[0].iloc[i] = iso.iloc[i]
#         # print iso.iloc[i]
#         # print df[0].iloc[i]
#     if (i >= cumul[0]) & (i < cumul[1]):
#         df[1].iloc[i-cumul[0]] = iso.iloc[i]
#     if (i >= cumul[1]) & (i < cumul[2]):
#         df[2].iloc[i-cumul[1]] = iso.iloc[i]
#     if (i >= cumul[2]) & (i < cumul[3]):
#         df[3].iloc[i-cumul[2]] = iso.iloc[i]
#     if (i >= cumul[3]) & (i < cumul[4]):
#         df[4].iloc[i-cumul[3]] = iso.iloc[i]
#     if (i >= cumul[4]) & (i < cumul[5]):
#         df[5].iloc[i-cumul[4]] = iso.iloc[i]
#     if (i >= cumul[5]) & (i < cumul[6]):
#         df[6].iloc[i-cumul[5]] = iso.iloc[i]
#     if (i >= cumul[6]) & (i < cumul[7]):
#         df[7].iloc[i-cumul[6]] = iso.iloc[i]
#     if (i >= cumul[7]) & (i < cumul[8]):
#         df[8].iloc[i-cumul[7]] = iso.iloc[i]
#     if (i >= cumul[8]) & (i < cumul[9]):
#         df[9].iloc[i-cumul[8]] = iso.iloc[i]
#     if (i >= cumul[9]) & (i < cumul[10]):
#         df[10].iloc[i-cumul[9]] = iso.iloc[i]
#     if (i >= cumul[10]) & (i < cumul[11]):
#         df[11].iloc[i-cumul[10]] = iso.iloc[i]
#     if (i >= cumul[11]) & (i < cumul[12]):
#         df[12].iloc[i-cumul[11]] = iso.iloc[i]
#     if (i >= cumul[12]) & (i < cumul[13]):
#         df[13].iloc[i-cumul[12]] = iso.iloc[i]
#     if (i >= cumul[13]) & (i < cumul[14]):
#         df[14].iloc[i-cumul[13]] = iso.iloc[i]
#     if (i >= cumul[14]) & (i < cumul[15]):
#         df[15].iloc[i-cumul[14]] = iso.iloc[i]
#     if (i >= cumul[15]) & (i < cumul[16]):
#         df[16].iloc[i-cumul[15]] = iso.iloc[i]

# print df[0]

''' Reading in of PARAM output files for the K2 GAP (no spectro) '''
BC3G = BC6G = EC6G = SC3G = SC6G = YC3G = YC6G = C3G = C6G = pd.DataFrame()
files_in = [BC3G,BC6G,EC6G,SC3G,SC6G,YC3G,YC6G]
file_names = ['BC3_GAP','BC6_GAP','EC6_GAP','SC3_GAP','SC6_GAP','YC3_GAP','YC6_GAP']
BC3G,BC6G,EC6G,SC3G,SC6G,YC3G,YC6G = p.param_in(files_in,file_names)
C3G,C6G = p.param_in([C3G,C6G],['C3_comp','C6_comp'])

# print len(C3G), len(C6G)

''' Read in pre-PARAM run files. Edit names depending upon use of data with or
without spectroscopy '''
BC3 = BC6 = EC6 = SC3 = SC6 = YC3 = YC6 = C3 = C6 = pd.DataFrame()
files_in2 = [BC3,BC6,EC6,SC3,SC6,YC3,YC6]
file_names2 = ['BC3_noSpec_ready','BC6_noSpec_ready','EC6_noSpec_ready','SC3_noSpec_ready', \
                'SC6_noSpec_ready','YC3_noSpec_ready','YC6_noSpec_ready']
BC3,BC6,EC6,SC3,SC6,YC3,YC6 = p.param_inputs(files_in2,file_names2)
C3,C6 = p.param_inputs_csv([C3,C6],['C3','C6'])



''' Spectroscopic run read in '''
BCs3 = BCs6 = ECs6 = SCs3 = SCs6 = YCs3 = YCs6 = pd.DataFrame()
files_in3 = [BCs3,BCs6,SCs3,SCs6,YCs3,YCs6]
file_names3 = ['Benoit_C3','Benoit_C6','Savita_C3','Savita_C6', \
               'Yvonne_C3','Yvonne_C6']
BCs3,BCs6,SCs3,SCs6,YCs3,YCs6 = p.param_in(files_in3,file_names3)
RAVE3 = pd.read_csv('/home/bmr135/GA/K2Poles/RAVE_C3.csv')
RAVE6 = pd.read_csv('/home/bmr135/GA/K2Poles/RAVE_C6.csv')
GES = pd.read_csv('/home/bmr135/Dropbox/GES-K2/param_tables/params_abund.txt',delimiter=r'\s+')
# print len(RAVE3), len(RAVE6)
# print BCs3.columns.values

besa3, besa6 = dat.BESANCON() # Read in Besacon simulations
besa3,besa6 = p.selection_function([besa3,besa6],['numax','numax'])
# besa3 = besa3[besa3['zgal'] < 6]
besa6 = besa6[besa6['zgal'] < 6]

''' Uncertainty calculation '''
uncerts = ['e_age1','e_age2','e_mass1','e_mass2','e_logg1','e_logg2','e_rad1','e_rad2', \
           'e_logrho1','e_logrho2','e_Kepler1','e_Kepler2','e_dist1','e_dist2']
p_in = ['age','mass','logg','logrho','rad','Kepler','dist']
p_ins = ['age','massp','logg','logrho','rad','Kepler','dist']
BCG3,BC6G,EC6G,SC3G,SC6G,YC3G,YC6G = p.uncerts(p_in,uncerts,files_in) # caclulate uncertainties on p_ins parameters
BCs3,BCs6,SCs3,SCs6,YCs3,YCs6 = p.uncerts(p_ins,uncerts,files_in3)
C3G,C6G = p.uncerts(p_ins,uncerts,[C3G,C6G])

CS3 = pd.merge(RAVE3,C3G,how='inner',on=['EPIC'])
CS3 = CS3.reset_index(drop=True)

CS6 = pd.merge(RAVE6,C6G,how='inner',on=['EPIC'])
CS6 = CS6.reset_index(drop=True)

# print len(CS3), len(CS6)



# print besa3.columns.values
# print C3.columns.values
# print C3G.columns.values

# prop.temp_JK(C3,C3['[Fe/H]'],C3['JK'])
# prop.temp_JK(C6,C6['[Fe/H]'],C6['JK'])

# print C3G['EPIC']
# print C3['EPIC']
Matched3 = pd.merge(C3,C3G,how='inner',on='EPIC')
Matched3s = pd.merge(C3,CS3,how='inner',on='EPIC')

Matched6 = pd.merge(C6,C6G,how='inner',on='EPIC')
Matched6s = pd.merge(C6,CS6,how='inner',on='EPIC')


# print Matched3.columns.values
# print GES.columns.values
# print Matched['mass_y'], SC3G['mass']
# Matched['spectro_mass'] = (Matched['Bnumax']/const.solar_Numax)**3 * (Matched['BDnu']/const.solar_Dnu)**-4 * \
#                     (Matched['Teff_N_K']/const.solar_Teff)**1.5

Matched3 = Matched3[Matched3['massp'] > 0]
# Matched3 = Matched3[Matched3['rad'] < 8]
# Matched3s = Matched3s[Matched3s['mass_y'] > 0]
# Matched3s = Matched3s[Matched3s['rad'] < 8]
CS3 = CS3[CS3['massp'] > 0]
# CS3 = CS3[CS3['rad'] < 8]
# CS3 = CS3[CS3['Teff_y'] > 1]
# CS3 = CS3[CS3['Teff_N_K'] > 1]
APO = APO[APO['mass'] > 0]
# APO = APO[APO['radius'] < 8]
Matched6 = Matched6[Matched6['massp'] > 0]
# Matched6 = Matched6[Matched6['rad'] < 8]
# Matched6s = Matched6s[Matched6s['mass_y'] > 0]
# Matched6s = Matched6s[Matched6s['rad'] < 8]
CS6 = CS6[CS6['massp'] > 0]
# CS6 = CS6[CS6['rad'] < 8]
# CS6 = CS6[CS6['Teff'] > 1]
# CS6 = CS6[CS6['Teff_N_K'] > 1]

# print len(CS3), len(CS6)



''' Least Squares fitting '''
#for least squares with x and y errors
# mpar, cpar, empar, ecpar = p.least_squares_2err(Matched)
#for least squares with errors on 1 axis
# mpar, cpar, empar, ecpar = p.least_squares_1err(Matched)
# x = np.arange(0,max(Matched['mass_x'])+0.1,0.1)
#
# plt.figure()
# plt.errorbar(Matched['mass_x'], Matched['mass_y'], yerr=[Matched['e_mass1'],Matched['e_mass2']], fmt='o', alpha=0.5)# xerr=[Matched['e_mass1_x'],Matched['e_mass2_x']],
# plt.xlabel(r'Original Mass')
# plt.ylabel(r'PARAM Mass')
# plt.plot([0,max(x)],[0,max(x)],c='r')
# plt.plot([0,6],[0.0+cpar[0],(6*mpar[0])+cpar[0]],color='r',linestyle='--')
# plt.plot(x,(x*mpar[0])+cpar[0],color='r',linestyle='--')
# plt.plot(x,(x*(mpar[0]+empar[0]))+(cpar[0]+ecpar[0]),color='m',linestyle='--')
# plt.plot(x,(x*(mpar[0]-empar[0]))+(cpar[0]-ecpar[0]),color='m',linestyle='--')
# plt.fill_between(x,(x*(mpar[0]-empar[0]))+(cpar[0]-ecpar[0]),(x*(mpar[0]+empar[0]))+(cpar[0]+ecpar[0]), \
    #            facecolor='m',edgecolor='red', alpha=0.5)
#
# print CS3.columns.values
#
# sel = pd.DataFrame()
# sel = Matched[Matched['L'] > 50]
# sel = Matched[Matched['Teff'] > 5000]
# a=np.arange(0.3,4.5,0.1)
# plt.figure()
# plt.scatter(Matched['mass_y'],Matched['age']*10**9)
# plt.scatter(sel['mass_y'],sel['age']*10**9,color='r')
# plt.plot(a,10**9*a**-3.2,linestyle='--')
# plt.xlabel(r'Mass')
# plt.ylabel(r'Age')

BESA = pd.concat([besa3,besa6],ignore_index=True)
BESA = BESA.reset_index(drop=True)


Matched3['mass_seismo'] = (Matched3['NUMAX']/3090)**3 * (Matched3['DNU']/135.5)**-4 * (Matched3['Teff_JK']/5777)**1.5
Matched6['mass_seismo'] = (Matched6['NUMAX']/3090)**3 * (Matched6['DNU']/135.5)**-4 * (Matched6['Teff_JK']/5777)**1.5

Matched3['mass_per'] = np.absolute((Matched3['massp'] - Matched3['mass_seismo'])/Matched3['mass_seismo'])
Matched6['mass_per'] = np.absolute((Matched6['massp'] - Matched6['mass_seismo'])/Matched6['mass_seismo'])

Matched3 = Matched3[Matched3['mass_per'] <= 0.3]
Matched6 = Matched6[Matched6['mass_per'] <= 0.3]

''' Mass vs Vertical Height '''
x = Matched3['dist']*np.sin(Matched3['Glat']*np.pi/180)*1e-3
CS3['xs'] = CS3['dist']*np.sin(CS3['Glat_x']*np.pi/180)*1e-3
CS3 = CS3[CS3['xs'] < .0]
CS3 = CS3.reset_index(drop=True)
plt.figure()
plt.scatter(BESA['IniMass'],BESA['zgal'],alpha=0.2,label=r'Sim',color='grey')
plt.scatter(APO['mass'],APO['Z']*1e-3,alpha=0.2,label=r'APOKASC',color='g')
# plt.scatter(besa3['IniMass'],besa3['zgal'],alpha=0.2,label=r'Sim C3',color='grey')
plt.scatter(Matched3['mass_seismo'],x,alpha=0.7,label=r'C3')#,c=Matched3['age']
# plt.scatter(CS3['mass'],CS3['xs'],alpha=1.0,label=r'C3') # ,c=np.log10(CS3['age'])
# plt.colorbar()
# plt.ylim(0,max(x)+1)
# plt.title(r'C3')
plt.xlabel(r'Mass, [$M_{\odot}$]')
# plt.ylabel(r'Z [kpc]')
# plt.legend()

y = Matched6['dist']*np.sin(Matched6['Glat']*np.pi/180)*1e-3
CS6['ys'] = CS6['dist']*np.sin(CS6['Glat_x']*np.pi/180)*1e-3
CS6 = CS6[CS6['ys'] > .0]
CS6 = CS6.reset_index(drop=True)
# plt.figure()
# plt.scatter(besa6['IniMass'],besa6['zgal'],alpha=0.2,label=r'Sim C6',color='grey')
plt.scatter(Matched6['mass_seismo'],y,alpha=0.7,label=r'C6',color='r')#,c=Matched6['age']
# plt.scatter(CS6['mass'],CS6['ys'],alpha=1.0,color='r',label=r'C6') # ,c=np.log10(CS6['age'])
# plt.colorbar()
# plt.ylim(0,max(y)+1)
# plt.title(r'C6')
# plt.xlabel(r'Mass, [$M_{\odot}$]')
plt.ylabel(r'Z [kpc]')
plt.legend(prop={'size':15})
plt.tight_layout()

''' Mass vs Temp '''
# plt.figure()
# plt.scatter(CS3['mass_y'],CS3['Teff_N_K'],label=r'RAVE C3')
# # plt.scatter(CS3['mass_y'],CS3['Teff_y'],color='r',label=r'C3 MAST')
# plt.scatter(CS6['mass_y'],CS6['Teff_N_K'],marker='s',color='g',label=r'RAVE C6')
# # plt.scatter(CS6['mass_y'],CS6['Teff'],color='m',marker='s',label=r'C6 MAST')
# plt.xlabel(r'Mass, [$M_{\odot}$]')
# plt.ylabel(r'Teff [K]')
# plt.legend()
# plt.tight_layout()

''' PARAM vs Scaling Mass '''
# CS3['mass_seismo'] = (CS3['Bnumax']/3090)**3 * (CS3['BDnu']/135.4)**-4 * (CS3['Teff_N_K']/5777)**1.5
# CS6['mass_seismo'] = (CS6['Bnumax']/3090)**3 * (CS6['BDnu']/135.4)**-4 * (CS6['Teff_N_K']/5777)**1.5
#
# plt.figure()
# plt.scatter(CS3['massp'],CS3['mass_seismo'],label=r'C3',alpha=0.5)
# plt.scatter(CS6['massp'],CS6['mass_seismo'],color='r',label=r'C6',alpha=0.5)
# plt.plot([0,3],[0,3])
# plt.xlabel(r'Mass - PARAM, [$M_{\odot}$]')
# plt.ylabel(r'Mass - Scaling, [$M_{\odot}$]')
# plt.legend()
# plt.tight_layout()


plt.figure()
plt.scatter(Matched3['massp'],Matched3['mass_seismo'],label=r'C3',alpha=0.5)
plt.scatter(Matched6['massp'],Matched6['mass_seismo'],color='r',label=r'C6',alpha=0.5)
plt.plot([0,3],[0,3])
plt.xlabel(r'Mass - PARAM, [$M_{\odot}$]')
plt.ylabel(r'Mass - Scaling, [$M_{\odot}$]')
plt.legend()
plt.tight_layout()

''' MAST vs RAVE Teff '''
# plt.figure()
# plt.scatter(CS3['Teff_y'],CS3['Teff_N_K'],label=r'C3',alpha=0.5)
# plt.scatter(CS6['Teff_y'],CS6['Teff_N_K'],color='r',label=r'C6',alpha=0.5)
# plt.plot([0,10000],[0,10000])
# plt.xlabel(r'Teff - MAST, [K]')
# plt.ylabel(r'Teff - RAVE, [K]')
# plt.legend()
# plt.tight_layout()

''' HRD '''
# plt.figure()
# plt.scatter(BESA['Teff'],BESA['L'],alpha=0.05,color='grey')
# for i in range(0,len(CS3),1):
#     if CS3['massp'][i] > 0.8:
#         plt.scatter(CS3['Teff_N_K'][i],CS3['L_y'][i],alpha=0.4)
#     elif CS3['massp'][i] <= 0.8:
#         plt.scatter(CS3['Teff_N_K'][i],CS3['L_y'][i],color='r')
# for i in range(0,len(CS3),1):
#     if CS6['massp'][i] > 0.8:
#         plt.scatter(CS6['Teff_N_K'][i],CS6['L_y'][i],alpha=0.4)
#     elif CS6['massp'][i] <= 0.8:
#         plt.scatter(CS6['Teff_N_K'][i],CS6['L_y'][i],color='g')
# plt.xlim(6000,3000)
# plt.xlabel(r'Teff [K]')
# plt.ylabel(r'L [L$_{\odot}$]')

Matched3['L'] = Matched3['rad']**2 * (Matched3['Teff_JK']/5777)**4
Matched6['L'] = Matched6['rad']**2 * (Matched6['Teff_JK']/5777)**4
APO['L'] = APO['radius']**2 * (APO['teff']/5777)**4
plt.figure()
plt.scatter(APO['teff'],APO['L'],color='grey',alpha=0.2,label=r'APOKASC')
plt.scatter(Matched3['Teff_JK'],Matched3['L'],color='r',label=r'C3')
plt.scatter(Matched6['Teff_JK'],Matched6['L'],color='g',label=r'C6')
plt.xlim(6000,3000)
plt.xlabel(r'Teff [K]')
plt.ylabel(r'L [L$_{\odot}$]')
plt.legend(prop={'size':15})
plt.tight_layout()


BESA['numax'] = BESA['IniMass'] * BESA['Radius']**-2 * (BESA['Teff']/5777)**-0.5 * 3090

''' Kiel Diagram '''
# plt.figure()
# plt.scatter(BESA['Teff'],BESA['numax'],alpha=0.05,color='grey')
# for i in range(0,len(CS3),1):
#     if CS3['mass'][i] > 0.8:
#         plt.scatter(CS3['Teff_N_K'][i],CS3['numax'][i],alpha=0.4)
#     elif CS3['mass'][i] <= 0.8:
#         plt.scatter(CS3['Teff_N_K'][i],CS3['numax'][i],color='r')
# for i in range(0,len(CS3),1):
#     if CS6['mass'][i] > 0.8:
#         plt.scatter(CS6['Teff_N_K'][i],CS6['Bnumax'][i],alpha=0.4)
#     elif CS6['mass'][i] <= 0.8:
#         plt.scatter(CS6['Teff_N_K'][i],CS6['Bnumax'][i],color='g')
# plt.xlim(6000,3000)
# plt.ylim(280,0)
# plt.xlabel(r'Teff [K]')
# plt.ylabel(r'$\nu_{\rm{max}}$ [$\mu$Hz]')

''' alpha vs Fe '''

RAVE3 = RAVE3[RAVE3['Alpha_c_1'] > -1.5]
RAVE6 = RAVE6[RAVE6['Alpha_c_1'] > -1.5]
GES = GES[GES['ALPHA'] > -1.5]

nullfmt = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.55
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_FEH = [left, bottom_h-0.1, width, 0.2]
rect_ALP = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(figsize=(10, 5))

axScatter = plt.axes(rect_scatter)
axFEH = plt.axes(rect_FEH)
axALP = plt.axes(rect_ALP)

# no labels
axFEH.xaxis.set_major_formatter(nullfmt)
axFEH.set_yticks([])
axALP.yaxis.set_major_formatter(nullfmt)
axALP.set_xticks([])

# the scatter plot:
axScatter.scatter(RAVE3['Met_N_K'],RAVE3['Alpha_c_1'],alpha=0.3,label=r'RAVE C3')
axScatter.scatter(RAVE6['Met_N_K'],RAVE6['Alpha_c_1'],alpha=0.3,label=r'RAVE C6',color='r')
axScatter.scatter(GES['FEH'],GES['ALPHA'],alpha=0.3,label=r'Gaia-ESO C3',color='g')

axScatter.set_xlabel(r'[Fe/H]')
axScatter.set_ylabel(r'[$\alpha$/Fe]')

# now determine nice limits by hand:
# binwidth = 0.25
# xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
# lim = (int(xymax/binwidth) + 1) * binwidth

# axScatter.set_xlim((-lim, lim))
# axScatter.set_ylim((-lim, lim))

bin_list = np.arange(-3, 1.5, 0.1)
axFEH.hist(RAVE3['Met_N_K'][np.isfinite(RAVE3['Met_N_K'])], bins=bin_list, normed=1, facecolor='b', label='RAVE3', histtype='step')
axFEH.hist(RAVE6['Met_N_K'][np.isfinite(RAVE6['Met_N_K'])], bins=bin_list, normed=1, facecolor='r', label='RAVE6', histtype='step')
axFEH.hist(GES['FEH'], bins=bin_list, normed=1, facecolor='g', label='Gaia-ESO', histtype='step')
axALP.hist(RAVE3['Alpha_c_1'][np.isfinite(RAVE3['Alpha_c_1'])], bins=bin_list, normed=1, facecolor='b', label='RAVE3',orientation='horizontal', histtype='step')
axALP.hist(RAVE6['Alpha_c_1'][np.isfinite(RAVE6['Alpha_c_1'])], bins=bin_list, normed=1, facecolor='r', label='RAVE6',orientation='horizontal', histtype='step')
axALP.hist(GES['ALPHA'], bins=bin_list, normed=1, facecolor='g', label='Gaia-ESO',orientation='horizontal', histtype='step')


# axHisty.hist(y, bins=bins, orientation='horizontal')

axFEH.set_xlim(axScatter.get_xlim())
axALP.set_ylim(axScatter.get_ylim())


# fig = plt.figure(figsize=(8,4.5))
# gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
# ax1 = fig.add_subplot(gs[0])
#
axScatter.legend(loc=3,prop={'size':15})
# # fig.ylim(-1,1)
#
# ax2 = fig.add_subplot(gs[1], sharex=ax1)
# bin_list = np.linspace(-3,1.5,51)
# n,bins,patches = ax2.hist(RAVE3['Met_N_K'][np.isfinite(RAVE3['Met_N_K'])], bins=bin_list, normed=1, facecolor='b',alpha=0.25, label='RAVE3')
# n1,bins1,patches1 = ax2.hist(RAVE6['Met_N_K'][np.isfinite(RAVE6['Met_N_K'])], bins=bin_list, normed=1, facecolor='r',alpha=0.25, label='RAVE6')
# n2,bins2,patches2 = ax2.hist(GES['FEH'], bins=bin_list, normed=1, facecolor='g',alpha=0.25, label='Gaia-ESO')

# fig.xlabel(r'[Fe/H]')
# fig.ylabel(r'[$\alpha$/Fe]')
# fig.subplots_adjust(wspace=0.0)
# plt.tight_layout()

''' numax vs dnu w/ mass colour bar and isochrones '''
m = [0.6,0.8,0.8,1.0,1.1,1.2,1.2,1.3,1.4,1.5,1.55,1.6,1.65,1.7,1.75,1.8,2.0]
for i in range(0,len(df),1):
    X = df[i]
    # print X
    X['R'] = np.sqrt((10**X['logL'])/((10**X['logT'])/5777)**4)
    X['numax'] = 3090 * m[i] * X['R']**-2 * ((10**X['logT'])/5777)**-0.5
    df[i] = X

# print df[0]['R']

K2 = pd.concat([Matched3,Matched6],ignore_index=True)
K2 = K2.reset_index(drop=True)

plt.figure()
plt.scatter(APO['numax'],APO['Dnu'],color='grey',alpha=0.2,label=r'APOKASC')
# plt.scatter(K2['NUMAX'],K2['DNU'],color='r',label=r'K2',alpha=0.3,facecolors='none')
plt.scatter(Matched3['NUMAX'],Matched3['DNU'],color='b',label=r'C3',alpha=0.3)#,facecolors='none')
plt.scatter(Matched6['NUMAX'],Matched6['DNU'],color='r',label=r'C6',alpha=0.3)#,facecolors='none')
plt.plot(df[0]['numax'],df[0]['dnu'],label=r'M$=0.6$',linewidth=2,color='k',linestyle='-.')
plt.plot(df[3]['numax'],df[3]['dnu'],label=r'M$=1.0$',linewidth=2,linestyle='-',color='k')
plt.plot(df[9]['numax'],df[9]['dnu'],label=r'M$=1.5$',linewidth=2,linestyle='--',color='k')
# plt.plot(df[14]['numax'],df[14]['dnu'],label=r'M$=1.75$',linewidth=2)
plt.plot(df[16]['numax'],df[16]['dnu'],label=r'M$=2.0$',linewidth=2,linestyle=':',color='k')

plt.xlim(0,280)
plt.ylim(0,20)
plt.xlabel(r'$\nu_{\rm{max}}$ [$\mu$Hz]',fontsize='20')
plt.ylabel(r'$\Delta\nu$ [$\mu$Hz]',fontsize='20')
plt.legend(loc=4,prop={'size':15})
plt.tight_layout()

plt.show()
