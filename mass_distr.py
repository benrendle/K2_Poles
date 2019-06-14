''' Mass Distribution Analysis for PARAM outputs '''
'''
Required to enter a number depending upon the set-up desired:
    - 0: Non-sampled simulation with clump
    - 1: Non-sampled simulation without clump
    - 2: Sampled simulation with clump
    - 3: Sampled simulation without clump
'''

import pandas as pd
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import scipy
from scipy import stats
import scipy.odr.odrpack as odrpack
from sklearn import linear_model, datasets
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import seaborn as sns
import colormaps
# from path import path
from pyqt_fit import kde
import sim_pert as sp
import copy
import mass_distr_functs as mdf
import K2_properties as prop
import K2_data as dat
import K2_constants as const
import subprocess
import TAR as tt
import time
import matplotlib.backends.backend_pdf
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
from matplotlib import lines
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from astropy.io import fits, ascii
from astropy.table import Table
from matplotlib.collections import LineCollection

mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams["font.family"] = "serif"


def ransac_fit(df,df1,param,label,f):#,uncert):
    ''' Using RANSAC fitting algorithm to fit the trends observed between different
        spectroscopic pipelines.

    Data in:
        - DataFrames containing the parameters to be compared (should already be
          the same length)
        - Parameters to be compared (string)
        - Uncertainty if applicable (string)
        - Axes labels (string)
        - f, factor to extend plot limits to in the x-axis (float)
    '''
    df['Diff'] = df[param[0]] - df1[param[1]]
    print(np.std(df['Diff']))
    # df['sig'] = np.sqrt(df[uncert[0]]**2 + df1[uncert[1]]**2)
    X = (df[param[0]])
    X = X.values.reshape(-1,1)
    y = (df['Diff'])

    ''' Fit line using all data '''
    lr = linear_model.LinearRegression()
    lr.fit(X, y)

    ''' Robustly fit linear model with RANSAC algorithm '''
    a = np.zeros((1000,2))
    b = np.zeros((1000,np.shape(df)[0])) # Size of inlier mask
    for i in range(1000):
        ransac = linear_model.RANSACRegressor(min_samples=5)
        ransac.fit(X, y)
        inlier_mask = ransac.inlier_mask_
        outlier_mask = np.logical_not(inlier_mask)

        ''' Predict data of estimated models '''
        line_X = np.arange(X.min()-f, X.max()+f, 0.01)#[:, np.newaxis]
        line_X = line_X.reshape(-1,1)
        line_y = lr.predict(line_X)
        line_y_ransac = ransac.predict(line_X)
        ''' Compare estimated coefficients '''
        a[i,0], a[i,1] = ransac.estimator_.coef_,ransac.estimator_.intercept_
        b[i,:] = inlier_mask

    ''' Make inliers boolean and find optimal fit using the median of the data '''
    b = b.astype('bool')
    c = np.median(a[:,0])
    medIdx = (np.abs(a[:,0] - c)).argmin()
    print(a[medIdx,0],a[medIdx,1])
    # print(np.median(a[:,0]),a[medIdx,0],b[medIdx])

    ''' Plot the inliers, outliers and optimal fit to the data '''
    plt.figure()
    lw = 2
    plt.plot(line_X, a[medIdx,0]*line_X + a[medIdx,1], color='k', linewidth=3, label=r'%.5s$*x$ + %.5s'%(a[medIdx,0],a[medIdx,1]), alpha=0.5)
    # plt.errorbar(X, y, yerr=df['sig'], color='gold', marker='.', label='Outliers',fmt='o')
    # plt.errorbar(X[b[medIdx]], y[b[medIdx]], yerr=df['sig'][b[medIdx]], color='yellowgreen', marker='.', label='Inliers',fmt='o')
    plt.scatter(X, y, color='gold', label='Outliers')
    plt.scatter(X[b[medIdx]], y[b[medIdx]], color='yellowgreen', label='Inliers')
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.xlim(X.min()-f, X.max()+f)
    plt.legend()
    df['outlier_flag'] = 1.
    df['outlier_flag'][b[medIdx]] = 0.
    # print(df['outlier_flag'])
    plt.show()
    sys.exit()


if __name__ == '__main__':

    '''
    PHASE 1: READ IN PRELIMINARY PARAM RUNS AND CROSS MATCH WITH INPUT DATA
    '''

    ext = '/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/'
    #
    ''' Read any results from tar files '''
    # a = [\
    #     ['.tgz', 'folder/', 'file.in.mo'], \
    #     ]
    #
    # z=[pd.DataFrame()]*len(a)
    # for i in range(len(a)):
    #     tar = tt.TAR(ext,a[i][0],a[i][1],a[i][2],r'\s+')
    #     z[i] = tar()
    #
    # df = z
    #
    ''' Read in files for updates to the results '''
    # APK = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/APOKASC4BEN.txt')
    # APK2 = pd.read_csv(ext+'APOKASC/param_in_June222017_150813.in.mo',delimiter=r'\s+')
    # #
    # TRI3 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C3_self')
    # # TRI3 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C3_old_stars.csv')
    # TRI6 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C6_self')
    #
    # ''' Read in files used with a Y enriched grid '''
    # C3en = pd.read_csv(ext+'K2.R3/C3_Andrea_27072018.in.mo',delimiter=r'\s+')
    # C6en = pd.read_csv(ext+'K2.R7/C6_Andrea_27072018.in.mo',delimiter=r'\s+')
    # RC3en = pd.read_csv(ext+'K2.R2/RC3_Andrea_27072018.in.mo',delimiter=r'\s+')
    # RC6en = pd.read_csv(ext+'K2.R5/RC6_Andrea_27072018.in.mo',delimiter=r'\s+')
    # GESen = pd.read_csv(ext+'K2.R6/GES_Andrea_27072018.in.mo',delimiter=r'\s+')
    # AP3en = pd.read_csv(ext+'K2.R9/AP3_Andrea_27072018.in.mo',delimiter=r'\s+')
    # AP6en = pd.read_csv(ext+'K2.R1/AP6_Andrea_27072018.in.mo',delimiter=r'\s+')
    # L3en = pd.read_csv(ext+'K2.R4/L3_Andrea_27072018.in.mo',delimiter=r'\s+')
    # L6en = pd.read_csv(ext+'K2.R8/L6_Andrea_27072018.in.mo',delimiter=r'\s+')
    #
    # ''' [Fe/H] corrected files '''
    # C3_feh = pd.read_csv(ext+'C3_FeH_200818/C3_FeH_20082018.in.mo',delimiter=r'\s+')
    # C6_feh = pd.read_csv(ext+'C6_FeH_200818/C6_FeH_20082018.in.mo',delimiter=r'\s+')
    #
    # ''' Skymapper (Luca Casagrande) [Fe/H] and Teff input results '''
    # C3_Luca = pd.read_csv(ext+'C3_Luca_impSig/C3_Luca.in.mo',delimiter=r'\s+')
    # Luca3 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_inputs/Poles/C3_Luca.in',delimiter=r'\s+')
    # # /media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/C6_Luca_impSig/C6_Luca.in.mo
    # C6_Luca = pd.read_csv(ext+'C6_Luca_impSig/C6_Luca.in.mo',delimiter=r'\s+')
    # Luca6 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_inputs/Poles/C6_Luca.in',delimiter=r'\s+')
    # Luca3en = pd.read_csv(ext+'/K2.R10/C3_Luca.in.mo',delimiter=r'\s+')
    # # Luca6en = pd.read_csv(ext+'K2.R11/C6_Luca.in.mo',delimiter=r'\s+')
    #
    #
    # ''' Kepler Simulation if required '''
    # Kep_Sim = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Standard_kepler_field/k1.6_K15.all.out.txt',delimiter=r'\s+')
    #
    # ''' Application of K2 Selection Function to Kepler Simulation '''
    # Kep_Sim['Teff'] = 10**(Kep_Sim['logTe'])
    # Kep_Sim['g'] = 10**(Kep_Sim['logg'])
    # Kep_Sim['L'] = 10**(Kep_Sim['logL'])
    # Kep_Sim['radius'] = np.sqrt(Kep_Sim['Mass'] / (Kep_Sim['g']/4.441))
    # Kep_Sim['JK'] = Kep_Sim['Jmag'] - Kep_Sim['Kmag']
    # Kep_Sim['Vcut'] = Kep_Sim['Kmag'] + 2*(Kep_Sim['JK']+0.14) + 0.382*np.exp(2*(Kep_Sim['JK']-0.2))
    # # Kep_Sim = prop.det_prob_Kepler(Kep_Sim,'numax',3090.0,135.1)
    # # Kep_Sim = Kep_Sim[ (Kep_Sim['numax'] > 10) & (Kep_Sim['numax'] < 280) & \
    # # (Kep_Sim['Vcut'] > 9) & (Kep_Sim['Vcut'] < 15) & (Kep_Sim['JK'] > 0.5)]
    #
    # ''' [Fe/H] correct PARAM runs '''
    # C3_New = pd.read_csv(ext+'C3_Andrea_27072018.in.csv.mo',delimiter=r'\s+')
    # C6_New = pd.read_csv(ext+'C6_Andrea_27072018.in.csv.mo',delimiter=r'\s+')
    # GES_T = pd.read_csv(ext+'GES_Andrea_27072018.csv.mo',delimiter=r'\s+')
    # RC3_T = pd.read_csv(ext+'RC3_Andrea_27072018.csv.mo',delimiter=r'\s+')
    # RC6_T = pd.read_csv(ext+'RC6_Andrea_27072018.csv.mo',delimiter=r'\s+')
    # L3_T = pd.read_csv(ext+'L3_Andrea_27072018.csv.mo',delimiter=r'\s+')
    # L6_T = pd.read_csv(ext+'L6_Andrea_27072018.csv.mo',delimiter=r'\s+')
    # AP3_T = pd.read_csv(ext+'AP3_Andrea_27072018.csv.mo',delimiter=r'\s+')
    # AP6_T = pd.read_csv(ext+'AP6_Andrea_27072018.csv.mo',delimiter=r'\s+')
    # print(len(GES_T))
    # ''' All spectroscopic parameter values run at once with PARAM '''
    # C_three = pd.read_csv(ext+'C3AS_040718.in.mo',delimiter=r'\s+')
    # C_six = pd.read_csv(ext+'C6AS_040718.in.mo',delimiter=r'\s+')
    # # C3_New = pd.read_csv(ext+'C3_070718.in.mo',delimiter=r'\s+')
    # # C6_New = pd.read_csv(ext+'C6_070718.in.mo',delimiter=r'\s+')
    #
    #
    # ''' Additional PARAM inputs that aren't output values '''
    # C3_New = mdf.p_in('C3_Andrea_27072018',C3_New,'C3')
    # C6_New = mdf.p_in('C6_Andrea_27072018',C6_New,'C6')
    # GES = mdf.p_in('GES_Andrea_27072018',GES_T,'C3')
    # print(len(GES))
    # RC3 = mdf.p_in('RC3_Andrea_27072018',RC3_T,'C3')
    # RC3 = RC3.drop_duplicates(subset=['#Id'])
    # RC3 = RC3.dropna(subset=['#Id'])
    # L3 = mdf.p_in('L3_Andrea_27072018',L3_T,'C3')
    # AP3 = mdf.p_in('AP3_Andrea_27072018',AP3_T,'C3')
    # AP6 = mdf.p_in('AP6_Andrea_27072018',AP6_T,'C6')
    # RC6 = mdf.p_in('RC6_Andrea_27072018',RC6_T,'C6')
    # L6 = mdf.p_in('L6_Andrea_27072018',L6_T,'C6')
    #
    # C3en = mdf.p_in('C3_Andrea_27072018',C3en,'C3')
    # C6en = mdf.p_in('C6_Andrea_27072018',C6en,'C6')
    # GESen = mdf.p_in('GES_Andrea_27072018',GESen,'C3')
    # RC3en = mdf.p_in('RC3_Andrea_27072018',RC3en,'C3')
    # RC3en = RC3en.drop_duplicates(subset=['#Id'])
    # RC3en = RC3en.dropna(subset=['#Id'])
    # L3en = mdf.p_in('L3_Andrea_27072018',L3en,'C3')
    # AP3en = mdf.p_in('AP3_Andrea_27072018',AP3en,'C3')
    # AP6en = mdf.p_in('AP6_Andrea_27072018',AP6en,'C6')
    # RC6en = mdf.p_in('RC6_Andrea_27072018',RC6en,'C6')
    # L6en = mdf.p_in('L6_Andrea_27072018',L6en,'C6')
    #
    # C3_feh = mdf.p_in('C3_Andrea_27072018',C3_feh,'C3')
    # C6_feh = mdf.p_in('C6_Andrea_27072018',C6_feh,'C6')
    #
    # C_three = mdf.p_in('All_Spec_280618',C_three,'C3')
    # C_three = C_three.drop_duplicates(subset=['#Id'])
    # C_six = mdf.p_in('All_Spec_280618',C_six,'C6')
    # C_six = C_six.drop_duplicates(subset=['#Id'])
    #
    #
    #
    # ''' Addition of RAVE [alpha/Fe] '''
    # R3alp = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/matlab_in/RC3_25042018_154829.csv')
    # R6alp = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/matlab_in/RC6_25042018_154829.csv')
    # # print(R3alp.columns.values)
    # RC3 = pd.merge(RC3,R3alp[['#Id','ALPHA']],how='inner',on=['#Id'])
    # RC6 = pd.merge(RC6,R6alp[['#Id','ALPHA']],how='inner',on=['#Id'])
    #
    # # C6 = pd.concat([C6_1,C6_2],ignore_index=True)
    #
    # TRI3['age'] = (10**TRI3['logAge'])/1e9
    # TRI3['Vcut'] = TRI3['Kmag'] + 2*(TRI3['JK']+0.14) + 0.382*np.exp(2*(TRI3['JK']-0.2))
    # TRI6['age'] = (10**TRI6['logAge'])/1e9
    #
    # cols_to_use = APK.columns.difference(APK2.columns)
    # cols_to_use = cols_to_use.union(['#Id'])
    # APK2 = pd.merge(APK2,APK[cols_to_use],how='inner',on=['#Id'])
    # APK2['dist'] = APK['dist']
    # APK2['Z'] = APK['Z']
    # APK2 = APK2[APK2['dist'] > -99.9]
    #
    # ''' Adding coordinates to K2 simulations '''
    # mdf.sim_coords(TRI3,3)
    # mdf.sim_coords(TRI6,6)
    # mdf.sim_coords(Kep_Sim,'K')
    # mdf.sim_dist(TRI3)
    # mdf.sim_dist(TRI6)
    # mdf.sim_dist(Kep_Sim)
    # Kep_Sim = Kep_Sim[Kep_Sim['dist'] < 2e4]
    # mdf.vert_dist(TRI3)
    # mdf.vert_dist(TRI6)
    # mdf.vert_dist(Kep_Sim)
    #
    # mdf.vert_dist(C3_New)
    # mdf.vert_dist(C6_New)
    # mdf.vert_dist(APK2)
    # mdf.vert_dist(GES)
    # mdf.vert_dist(AP3)
    # mdf.vert_dist(AP6)
    # mdf.vert_dist(RC3)
    # mdf.vert_dist(RC6)
    # mdf.vert_dist(L3)
    # mdf.vert_dist(L6)
    # mdf.vert_dist(C_three)
    # mdf.vert_dist(C_six)
    #
    # mdf.vert_dist(C3en)
    # mdf.vert_dist(C6en)
    # mdf.vert_dist(GESen)
    # mdf.vert_dist(AP3en)
    # mdf.vert_dist(AP6en)
    # mdf.vert_dist(RC3en)
    # mdf.vert_dist(RC6en)
    # mdf.vert_dist(L3en)
    # mdf.vert_dist(L6en)
    #
    # mdf.vert_dist(C3_feh)
    # mdf.vert_dist(C6_feh)
    # print(len(GES))
    #
    # ''' Spectro tests '''
    # ges,c3 = mdf.spectro(GES,C3_New)
    # ap3,c3a = mdf.spectro(AP3,C3_New)
    # ap6,c6a = mdf.spectro(AP6,C6_New)
    # rc3,c3_R = mdf.spectro(RC3,C6_New)
    # rc6,c6_R = mdf.spectro(RC6,C6_New)
    # l3, c3_l = mdf.spectro(L3,C3_New)
    # l6, c6_l = mdf.spectro(L6,C6_New)
    # GR3,RG3 = mdf.spectro(GES,RC3)
    # AR6,RA6 = mdf.spectro(AP6,RC6)
    # RL6,LR6 = mdf.spectro(RC6,L6)
    # AL6,LA6 = mdf.spectro(AP6,L6)
    # AR3,RA3 = mdf.spectro(AP3,RC3)
    # AG3,GA3 = mdf.spectro(AP3,GES)
    # print(len(ges))
    #
    # ''' Gaia-ESO Uncertainty Comps '''
    # param = ['mass','rad','age']
    # x = [['mass',0.0],['rad',0.0],['age',0.0]]
    # err1 = mdf.uncerts(C3_New,param,x)
    # err2 = mdf.uncerts(GES,param,x)
    # err3 = mdf.uncerts(ges,param,x)
    # err4 = mdf.uncerts(c3,param,x)
    #
    # ''' Combine all spectroscopic data for each campaign '''
    # c_three = pd.concat([ap3,rc3,ges,l3],ignore_index=True)
    # c_three = c_three.drop_duplicates(subset=['#Id'],keep='first')
    # c_three.reset_index(drop=True)
    # c_six = pd.concat([ap6,rc6,l6],ignore_index=True)
    # c_six = c_six.drop_duplicates(subset=['#Id'],keep='first')
    # c_six.reset_index(drop=True)
    # AS = pd.concat([c_three,c_six],ignore_index=True)
    # AS.reset_index(drop=True)
    # # c_three.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/C3_spec')
    #
    #
    # # K2 = pd.concat([C3,C6],ignore_index=True)
    # # K2.reset_index(drop=True)
    # # K2.reset_index(drop=True)
    # K2_New = pd.concat([C3_New,C6_New],ignore_index=True)
    # K2_New.reset_index(drop=True)
    # TRI = pd.concat([TRI3,TRI6],ignore_index=True)
    # TRI.reset_index(drop=True)
    #
    # # K2['lgs'] = np.log10(27400 * (K2['nmx']/3090) * np.sqrt(K2['Teff']/5777))
    # K2_New['lgs'] = np.log10(27400 * (K2_New['nmx']/3090) * np.sqrt(K2_New['Teff']/5777))
    # AS['lgs'] = np.log10(27400 * (AS['nmx']/3090) * np.sqrt(AS['Teff']/5777))
    # # mdf.vert_dist(AS)
    # # AS['feh'] = AS['feh'] + 3
    #
    # ''' Photometric data for spectroscopic K2 stars '''
    # cols_to_use = c_three.columns.difference(C3_New.columns)
    # cols_to_use = cols_to_use.union(['#Id'])
    # C3_AS = pd.merge(C3_New,c_three[cols_to_use],how='inner',on=['#Id'])
    # C3_AS = C3_AS.drop_duplicates(subset=['#Id'],keep='first')
    # C3_AS.reset_index(drop=True)
    #
    # cols_to_use = c_six.columns.difference(C6_New.columns)
    # cols_to_use = cols_to_use.union(['#Id'])
    # C6_AS = pd.merge(C6_New,c_six[cols_to_use],how='inner',on=['#Id'])
    # C6_AS = C6_AS.drop_duplicates(subset=['#Id'],keep='first')
    # C6_AS.reset_index(drop=True)
    #
    # K2_AS = pd.concat([C3_AS,C6_AS],ignore_index=True)
    # K2_AS.reset_index(drop=True)
    #
    # # K2_dnu = pd.DataFrame()
    # # K2_dnu = K2[K2['dnu'] > 0]
    # # K2_dnu.reset_index(drop=True)
    #
    # AS_red = AS[AS['logAge'] < 9.8]
    # AS_red = AS_red[AS_red['logAge'] > 9.6]
    # AS_red.reset_index(drop=True)
    #
    # ''' File lengths '''
    # # print(len(C3_New),len(C6_New))
    # # print(len(c_three),len(c_six))
    # # print(len(RC3),len(RC6))
    # print(len(GES))
    # # print(len(L3),len(L6))
    # # print(len(AP3),len(AP6))
    # # print(len(K2_New),len(AS))
    # # print(len(C3_Luca),len(C6_Luca))
    # # print(len(APK2))
    # print("-----")
    #
    # ''' Quality Flagging '''
    # ext_load = '/home/bmr135/K2_Poles/Mass_Distr_In/'
    # # ext_load = '/home/ben/K2_Poles/Mass_Distr_In/'
    # C3_skew = pd.read_csv(ext_load+'Additional_Tests/C3_EPIC_skew')
    # C6_skew = pd.read_csv(ext_load+'Additional_Tests/C6_EPIC_skew')
    #
    # C3_New['Mass_Flag'] = 0
    # C3_New['Radius_Flag'] = 0
    # C3_New['Age_Flag'] = 0
    # C3_New['mass_comp'] = abs((C3_New['mass'] - C3_New['m_seismo'])/C3_New['mass'])
    # for i in range(len(C3_New)):
    #     if C3_New['mass_comp'].iloc[i] > 0.5:
    #         C3_New['Mass_Flag'].iloc[i] == 1
    #     if abs((C3_New['rad'].iloc[i] - C3_New['r_seismo'].iloc[i])/C3_New['rad'].iloc[i]) > 0.5:
    #         C3_New['Radius_Flag'].iloc[i] == 1
    # for i in range(len(C3_New)):
    #     for j in range(len(C3_skew)):
    #         if C3_New['#Id'].iloc[i] == C3_skew['#ID'].iloc[j]:
    #             C3_New['Age_Flag'].iloc[i] = 1
    #
    # C6_New['Mass_Flag'] = 0
    # C6_New['Radius_Flag'] = 0
    # C6_New['Age_Flag'] = 0
    # for i in range(len(C6_New)):
    #     if abs((C6_New['mass'].iloc[i] - C6_New['m_seismo'].iloc[i])/C6_New['mass'].iloc[i]) > 0.5:
    #         C6_New['Mass_Flag'].iloc[i] == 1
    #     if abs((C6_New['rad'].iloc[i] - C6_New['r_seismo'].iloc[i])/C6_New['rad'].iloc[i]) > 0.5:
    #         C6_New['Radius_Flag'].iloc[i] == 1
    # for i in range(len(C6_New)):
    #     for j in range(len(C6_skew)):
    #         if C6_New['#Id'].iloc[i] == C6_skew['#ID'].iloc[j]:
    #             C6_New['Age_Flag'].iloc[i] = 1
    #
    # K2_New['Mass_Flag'] = 0
    # K2_New['Radius_Flag'] = 0
    # K2_New['Age_Flag'] = 0
    # for i in range(len(K2_New)):
    #     if abs((K2_New['mass'].iloc[i] - K2_New['m_seismo'].iloc[i])/K2_New['mass'].iloc[i]) > 0.5:
    #         K2_New['Mass_Flag'].iloc[i] == 1
    #     if abs((K2_New['rad'].iloc[i] - K2_New['r_seismo'].iloc[i])/K2_New['rad'].iloc[i]) > 0.5:
    #         K2_New['Radius_Flag'].iloc[i] == 1
    # for i in range(len(K2_New)):
    #     for j in range(len(C3_skew)):
    #         if K2_New['#Id'].iloc[i] == C3_skew['#ID'].iloc[j]:
    #             K2_New['Age_Flag'].iloc[i] = 1
    #     for j in range(len(C3_skew)):
    #         if K2_New['#Id'].iloc[i] == C6_skew['#ID'].iloc[j]:
    #             K2_New['Age_Flag'].iloc[i] = 1
    #
    # ''' Luca Photometry Additional Values '''
    # C3_Luca = pd.merge(C3_Luca,Luca3[['#Id','GLON','GLAT','feh','efeh','teff','eteff','numax','Dnu']],how='inner',on=['#Id'])
    # C3_Luca = C3_Luca.reset_index(drop=True)
    # Luca3en = pd.merge(Luca3en,Luca3[['#Id','GLON','GLAT','feh','efeh','teff','eteff','numax','Dnu']],how='inner',on=['#Id'])
    # Luca3en = Luca3en.reset_index(drop=True)
    # C6_Luca = pd.merge(C6_Luca,Luca6[['#Id','GLON','GLAT','feh','efeh','teff','eteff','numax','Dnu']],how='inner',on=['#Id'])
    # C6_Luca = C6_Luca.reset_index(drop=True)
    # # Luca6en = pd.merge(Luca6en,Luca6[['#Id','GLON','GLAT','feh','efeh','teff','eteff','numax','Dnu']],how='inner',on=['#Id'])
    # # Luca6en = Luca6en.reset_index(drop=True)
    # C3_Luca['Z'] = C3_Luca['dist']*np.sin(C3_Luca['GLAT']*np.pi/180)*1e-3
    # Luca3en['Z'] = Luca3en['dist']*np.sin(Luca3en['GLAT']*np.pi/180)*1e-3
    # C6_Luca['Z'] = C6_Luca['dist']*np.sin(C6_Luca['GLAT']*np.pi/180)*1e-3
    # # Luca6en['Z'] = Luca6en['dist']*np.sin(Luca6en['GLAT']*np.pi/180)*1e-3
    # C3_Luca['m_seismo'] = (C3_Luca['numax']/3090)**3 * (C3_Luca['Dnu']/135.1)**-4 * (C3_Luca['teff']/5777)**1.5
    # C3_Luca['r_seismo'] = (C3_Luca['numax']/3090) * (C3_Luca['Dnu']/135.1)**2 * (C3_Luca['teff']/5777)**0.5
    # Luca3en['m_seismo'] = (Luca3en['numax']/3090)**3 * (Luca3en['Dnu']/135.1)**-4 * (Luca3en['teff']/5777)**1.5
    # Luca3en['r_seismo'] = (Luca3en['numax']/3090) * (Luca3en['Dnu']/135.1)**2 * (Luca3en['teff']/5777)**0.5
    # C6_Luca['m_seismo'] = (C6_Luca['numax']/3090)**3 * (C6_Luca['Dnu']/135.1)**-4 * (C6_Luca['teff']/5777)**1.5
    # C6_Luca['r_seismo'] = (C6_Luca['numax']/3090) * (C6_Luca['Dnu']/135.1)**2 * (C6_Luca['teff']/5777)**0.5
    # # Luca6en['m_seismo'] = (Luca6en['numax']/3090)**3 * (Luca6en['Dnu']/135.1)**-4 * (Luca6en['teff']/5777)**1.5
    # # Luca6en['r_seismo'] = (Luca6en['numax']/3090) * (Luca6en['Dnu']/135.1)**2 * (Luca6en['teff']/5777)**0.5
    #
    # C3_Luca['percentage'] = 100*((C3_Luca['m_seismo']/C3_Luca['mass'])-1)
    # C6_Luca['percentage'] = 100*((C6_Luca['m_seismo']/C6_Luca['mass'])-1)
    # Luca3en['percentage'] = 100*((Luca3en['m_seismo']/Luca3en['mass'])-1)
    # # Luca6en['percentage'] = 100*((Luca6en['m_seismo']/Luca6en['mass'])-1)
    # c_three['percentage'] = 100*((c_three['m_seismo']/c_three['mass'])-1)
    # c_six['percentage'] = 100*((c_six['m_seismo']/c_six['mass'])-1)
    #
    #
    # ''' Merging RAVE and RAVE with additional Y grid '''
    # ind1, ind2 = pd.DataFrame(),pd.DataFrame()
    # ind1['#Id'] = RC3['#Id']
    # RC3en = pd.merge(RC3en,ind1,how='inner',on=['#Id'])
    # RC3en.reset_index(drop=True)
    # ind2['#Id'] = RC6['#Id']
    # RC6en = pd.merge(RC6en,ind2,how='inner',on=['#Id'])
    # RC6en.reset_index(drop=True)
    #
    # ASen = pd.concat([RC3en,RC6en,AP3en,AP6en,GESen,L3en,L6en],ignore_index=True)
    # ASen = ASen.drop_duplicates(subset=['#Id'],keep='first')
    # ASen = ASen.reset_index(drop=True)
    #
    # ''' Photometric C3/C6 values with Spectroscopic Coverage '''
    # # C3_pa = pd.merge(C3_New,c_three[['#Id']],time("%d%m%Y"),index=False)
    # C3_pa = pd.read_csv(ext_load+'Normal/Jun_2018/C3_AS_21062018')
    # C6_pa = pd.read_csv(ext_load+'Normal/Jun_2018/C6_AS_21062018')
    #
    # ''' Addition of RAVE [alpha/Fe] '''
    # # R3alp = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/matlab_in/RC3_25042018_154829.csv')
    # # R6alp = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/matlab_in/RC6_25042018_154829.csv')
    # # RC3 = pd.merge(RC3,R3alp[['#Id','ALPHA']],how='inner',on=['#Id'])
    # # RC6 = pd.merge(RC6,R6alp[['#Id','ALPHA']],how='inner',on=['#Id'])
    # # RC3.to_csv(ext_load+'RC3_19042018',index=False)
    # # RC6.to_csv(ext_load+'RC6_19042018',index=False)
    # # print(C6_New['Z'])
    # # print(K2_New.columns.values)
    # # sys.exit()
    #
    # ''' Form a single [alpha/Fe] column '''
    # for i in range(len(AS)):
    #         if AS['alpha'].iloc[i] == -99.9:
    #             AS['alpha'].iloc[i] = AS['ALPHA'].iloc[i]
    #
    # for i in range(len(c_three)):
    #     for j in range(len(RC3)):
    #         if c_three['#Id'].iloc[i] == RC3['#Id'].iloc[j]:
    #             c_three['alpha'].iloc[i] = RC3['ALPHA'].iloc[j]
    #
    # for i in range(len(c_six)):
    #     for j in range(len(RC6)):
    #         if c_six['#Id'].iloc[i] == RC6['#Id'].iloc[j]:
    #             c_six['alpha'].iloc[i] = RC6['ALPHA'].iloc[j]
    #
    # TRI = pd.concat([TRI3,TRI6],ignore_index=True)
    # TRI.reset_index(drop=True)
    #
    # GAP3, GAP6 = dat.K2_GAP()
    # GAP3 = GAP3[['EPIC','[Fe/H]']]
    # GAP6 = GAP6[['EPIC','[Fe/H]']]
    #
    # K2_GAP = pd.concat([GAP3,GAP6],ignore_index=True)
    # K2_GAP.reset_index(drop=True)
    # K2_GAP.rename(columns={'EPIC':'#Id'},inplace=True)
    # K2_New = pd.merge(K2_New,K2_GAP,how='inner',on=['#Id'])
    # K2_New = K2_New.reset_index(drop=True)
    # C3_New = pd.merge(C3_New,K2_GAP,how='inner',on=['#Id'])
    # C3_New = C3_New.reset_index(drop=True)
    # # sys.exit()
    #
    # ''' PARAM based mass cuts if desired '''
    # # C3_New = C3_New[abs((C3_New['mass'] - C3_New['m_seismo'])/C3_New['mass']) < 0.1]
    # # C3_New = C3_New.reset_index(drop=True)
    # # C6_New = C6_New[abs((C6_New['mass'] - C6_New['m_seismo'])/C6_New['mass']) < 0.1]
    # # C6_New = C6_New.reset_index(drop=True)
    # # print(len(C3_New),len(C6_New))
    #
    # ''' Error propagation - Galactic Radius and Z '''
    # APK2['sig_Z'] = 1e-3 * abs(APK2['dist_68U'] - APK2['dist_68L'])/2
    # C3_New['sig_Z'] = 1e-3 * abs(C3_New['dist_68U'] - C3_New['dist_68L'])/2
    # C6_New['sig_Z'] = 1e-3 * abs(C6_New['dist_68U'] - C6_New['dist_68L'])/2
    # K2_New['sig_Z'] = 1e-3 * abs(K2_New['dist_68U'] - K2_New['dist_68L'])/2
    # AS['sig_Z'] = 1e-3 * abs(AS['dist_68U'] - AS['dist_68L'])/2
    #
    # APK2['sig_Gal_Rad'] = np.sqrt((2*(APK2['sig_Z']/APK2['X'])**2) + (2*(APK2['sig_Z']/APK2['Y'])**2))
    # C3_New['sig_Gal_Rad'] = np.sqrt((2*(C3_New['sig_Z']/C3_New['X'])**2) + (2*(C3_New['sig_Z']/C3_New['Y'])**2))
    # C6_New['sig_Gal_Rad'] = np.sqrt((2*(C6_New['sig_Z']/C6_New['X'])**2) + (2*(C6_New['sig_Z']/C6_New['Y'])**2))
    # K2_New['sig_Gal_Rad'] = np.sqrt((2*(K2_New['sig_Z']/K2_New['X'])**2) + (2*(K2_New['sig_Z']/K2_New['Y'])**2))
    # AS['sig_Gal_Rad'] = np.sqrt((2*(AS['sig_Z']/AS['X'])**2) + (2*(AS['sig_Z']/AS['Y'])**2))
    #
    # ''' Data arrays '''
    # data = [C3_New, C6_New, K2_New, K2_AS, AS, c_three, c_six, RC3, RC6, L3, L6, GES, AP3, AP6]
    # data_en = [C3en, C6en, RC3en, RC6en, GESen, AP3en, AP6en, L3en, L6en]
    # data_feh = [C3_feh, C6_feh]
    # data_Luca = [C3_Luca, Luca3en, C6_Luca]
    #
    # ''' Parameter Uncertainties '''
    # data = mdf.uncert(['age','mass'],data)
    # # print(len(AP3),len(AP6))
    # # print(len(K2_New),len(AS))
    # # print(len(C3_Luca),len(C6_
    # data_en = mdf.uncert(['age','mass'],data_en)
    # data_feh = mdf.uncert(['age','mass'],data_feh)
    # data_Luca = mdf.uncert(['age','mass'],data_Luca)
    # APK2 = mdf.uncert(['age','mass'],[APK2])[0]
    #
    # ''' File lengths '''
    # # print(len(C3_New),len(C6_New))
    # # print(len(c_three),len(c_six))
    # # print(len(RC3),len(RC6))
    # print(len(GES))
    # # print(len(L3),len(L6))
    # # print(len(AP3),len(AP6))
    # # print(len(K2_New),len(AS))
    # # print(len(C3_Luca),len(C6_Luca))
    # # print(len(APK2))
    # print("-----")
    #
    # ''' Gaia cuts from detection probability tests '''
    # K2_Gaia = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/K2_det_prob_gaia')
    # K2_Gaia.rename(columns={'EPIC':'#Id'}, inplace=True)
    # data = mdf.gaiaCut(data,K2_Gaia)
    # data_en = mdf.gaiaCut(data_en,K2_Gaia)
    # data_feh = mdf.gaiaCut(data_feh,K2_Gaia)
    # data_Luca = mdf.gaiaCut(data_Luca,K2_Gaia)
    #
    # C3_New, C6_New, K2_New, K2_AS, AS, c_three, c_six, RC3, RC6, L3, L6, GES, AP3, AP6 = data
    # C3en, C6en, RC3en, RC6en, GESen, AP3en, AP6en, L3en, L6en = data_en
    # C3_feh, C6_feh = data_feh
    # C3_Luca, Luca3en, C6_Luca = data_Luca
    #
    # ''' File lengths '''
    # # print(len(C3_New),len(C6_New))
    # # print(len(c_three),len(c_six))
    # # print(len(RC3),len(RC6))
    # print(len(GES))
    # # print(len(L3),len(L6))
    # # print(len(AP3),len(AP6))
    # # print(len(K2_New),len(AS))
    # # print(len(C3_Luca),len(C6_Luca))
    # # sys.exit()
    # # ''' Save files with uncertainties included '''
    # ext_save = '/home/bmr135/K2_Poles/Mass_Distr_In/'
    # # ext_save = '/home/ben/K2_Poles/Mass_Distr_In/'
    # C3_New.to_csv(ext_save+'Normal/Feb_2019/C3_'+time.strftime("%d%m%Y"),index=False)
    # C6_New.to_csv(ext_save+'Normal/Feb_2019/C6_'+time.strftime("%d%m%y"),index=False)
    # K2_New.to_csv(ext_save+'Normal/Feb_2019/K2_'+time.strftime("%d%m%y"),index=False)
    # K2_AS.to_csv(ext_save+'Normal/Feb_2019/K2_AS_'+time.strftime("%d%m%y"),index=False)
    # AS.to_csv(ext_save+'Normal/Feb_2019/AS_'+time.strftime("%d%m%y"),index=False)
    # APK2.to_csv(ext_save+'Normal/Feb_2019/APOKASC_'+time.strftime("%d%m%y"),index=False)
    # c_three.to_csv(ext_save+'Normal/Feb_2019/Spec_C3_'+time.strftime("%d%m%y"),index=False)
    # c_six.to_csv(ext_save+'Normal/Feb_2019/Spec_C6_'+time.strftime("%d%m%y"),index=False)
    # RC3.to_csv(ext_save+'Normal/Feb_2019/RC3_'+time.strftime("%d%m%y"),index=False)
    # RC6.to_csv(ext_save+'Normal/Feb_2019/RC6_'+time.strftime("%d%m%y"),index=False)
    # L3.to_csv(ext_save+'Normal/Feb_2019/L3_'+time.strftime("%d%m%y"),index=False)
    # L6.to_csv(ext_save+'Normal/Feb_2019/L6_'+time.strftime("%d%m%y"),index=False)
    # GES.to_csv(ext_save+'Normal/Feb_2019/GES_'+time.strftime("%d%m%y"),index=False)
    # AP3.to_csv(ext_save+'Normal/Feb_2019/AP3_'+time.strftime("%d%m%y"),index=False)
    # AP6.to_csv(ext_save+'Normal/Feb_2019/AP6_'+time.strftime("%d%m%y"),index=False)
    # TRI3.to_csv(ext_save+'Normal/Feb_2019/TRI3_'+time.strftime("%d%m%y"),index=False)
    # TRI6.to_csv(ext_save+'Normal/Feb_2019/TRI6_'+time.strftime("%d%m%y"),index=False)
    #
    # # C3en.to_csv(ext_save+'Y_enhanced_grid/C3_'+time.strftime("%d%m%y"),index=False)
    # # C6en.to_csv(ext_save+'Y_enhanced_grid/C6_'+time.strftime("%d%m%y"),index=False)
    # # RC3en.to_csv(ext_save+'Y_enhanced_grid/RC3_'+time.strftime("%d%m%y"),index=False)
    # # RC6en.to_csv(ext_save+'Y_enhanced_grid/RC6_'+time.strftime("%d%m%y"),index=False)
    # # GESen.to_csv(ext_save+'Y_enhanced_grid/GES_'+time.strftime("%d%m%y"),index=False)
    # # AP3en.to_csv(ext_save+'Y_enhanced_grid/AP3_'+time.strftime("%d%m%y"),index=False)
    # # AP6en.to_csv(ext_save+'Y_enhanced_grid/AP6_'+time.strftime("%d%m%y"),index=False)
    # # L3en.to_csv(ext_save+'Y_enhanced_grid/L3_'+time.strftime("%d%m%y"),index=False)
    # # L6en.to_csv(ext_save+'Y_enhanced_grid/L6_'+time.strftime("%d%m%y"),index=False)
    #
    # C3_feh.to_csv(ext_save+'Additional_Tests/C3_FeH_'+time.strftime("%d%m%y"),index=False)
    # C6_feh.to_csv(ext_save+'Additional_Tests/C6_FeH_'+time.strftime("%d%m%y"),index=False)
    #
    # C3_Luca.to_csv(ext_save+'Luca_photom/C3_gaia'+time.strftime("%d%m%y"),index=False)
    # C6_Luca.to_csv(ext_save+'Luca_photom/C6_gaia'+time.strftime("%d%m%y"),index=False)
    # # Luca3en.to_csv(ext_save+'Luca_photom/C3en'+time.strftime("%d%m%y"),index=False)
    # # Luca6en.to_csv(ext_save+'Luca_photom/C6en'+time.strftime("%d%m%y"),index=False)
    #
    # sys.exit()


    '''
    PHASE 2 - READ IN RESULTS FOR ANY PARAMETER CUTS AND PLOTTING
    '''

    ext_load = '/home/bmr135/K2_Poles/Mass_Distr_In/'
    # ext_load = '/home/ben/K2_Poles/Mass_Distr_In/'
    ''' Read in processed results files '''
    C3_New = pd.read_csv(ext_load+'Normal/Feb_2019/C3_270219')
    C6_New = pd.read_csv(ext_load+'Normal/Feb_2019/C6_270219')
    K2_New = pd.read_csv(ext_load+'Normal/Feb_2019/K2_270219')
    K2_AS = pd.read_csv(ext_load+'Normal/Feb_2019/K2_AS_270219')
    AS = pd.read_csv(ext_load+'Normal/Feb_2019/AS_270219')
    APK2 = pd.read_csv(ext_load+'Normal/Feb_2019/APOKASC_270219')
    c_three = pd.read_csv(ext_load+'Normal/Feb_2019/Spec_C3_270219')
    c_six = pd.read_csv(ext_load+'Normal/Feb_2019/Spec_C6_270219')
    RC3 = pd.read_csv(ext_load+'Normal/Feb_2019/RC3_270219')
    RC6 = pd.read_csv(ext_load+'Normal/Feb_2019/RC6_270219')
    L3 = pd.read_csv(ext_load+'Normal/Feb_2019/L3_270219')
    L6 = pd.read_csv(ext_load+'Normal/Feb_2019/L6_270219')
    GES = pd.read_csv(ext_load+'Normal/Feb_2019/GES_270219')
    AP3 = pd.read_csv(ext_load+'Normal/Feb_2019/AP3_270219')
    AP6 = pd.read_csv(ext_load+'Normal/Feb_2019/AP6_270219')
    TRI3 = pd.read_csv(ext_load+'Normal/Feb_2019/TRI3_270219')
    TRI6 = pd.read_csv(ext_load+'Normal/Feb_2019/TRI6_270219')
    TRI = pd.concat([TRI3,TRI6],ignore_index=True)

    C3_Sky = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Gaia/matched_GAP_SkyMapper.csv')
    # C3_Sky = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/Gaia/matched_GAP_SkyMapper.csv')
    C3_Sky = pd.merge(C3_Sky,C3_New[['#Id']],how='inner',on=['#Id'])
    # Gaia_IDs = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Gaia/GAP_Gaia_Luca.csv')
    # Gaia_IDs = pd.merge(Gaia_IDs, K2_New[['#Id']], how='inner', on=['#Id'])
    # Gaia_IDs.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Gaia/GAP_Gaia_Luca.csv',index=False)

    ''' Read in files used with a Y enriched grid '''
    C3en = pd.read_csv(ext_load+'Y_enhanced_grid/C3_281118')
    C6en = pd.read_csv(ext_load+'Y_enhanced_grid/C6_281118')
    RC3en = pd.read_csv(ext_load+'Y_enhanced_grid/RC3_281118')
    RC6en = pd.read_csv(ext_load+'Y_enhanced_grid/RC6_281118')
    GESen = pd.read_csv(ext_load+'Y_enhanced_grid/GES_281118')
    AP3en = pd.read_csv(ext_load+'Y_enhanced_grid/AP3_281118')
    AP6en = pd.read_csv(ext_load+'Y_enhanced_grid/AP6_281118')
    L3en = pd.read_csv(ext_load+'Y_enhanced_grid/L3_281118')
    L6en = pd.read_csv(ext_load+'Y_enhanced_grid/L6_281118')

    ASen = pd.concat([RC3en,RC6en,GESen,AP3en,AP6en,L3en,L6en],ignore_index=True)
    ASen = ASen.drop_duplicates(subset=['#Id']).reset_index(drop=True)

    ''' [Fe/H] corrected files '''
    C3_feh = pd.read_csv(ext_load+'Additional_Tests/C3_FeH_281118')
    C6_feh = pd.read_csv(ext_load+'Additional_Tests/C6_FeH_281118')

    ''' Skymapper (Luca Casagrande) [Fe/H] and Teff input results

        - IS -> Improved Simga: PARAM run determined with improved uncertainties
    '''
    C3_Luca = pd.read_csv(ext_load+'Luca_photom/C3_gaia270219')
    C3_LucaIS = pd.read_csv(ext_load+'Luca_photom/C3_impSig280918')
    C6_Luca = pd.read_csv(ext_load+'Luca_photom/C6_gaia270219')
    C6_LucaIS = pd.read_csv(ext_load+'Luca_photom/C6_impSig280918')
    Luca3en = pd.read_csv(ext_load+'Luca_photom/C3en050918')
    # Luca6en = pd.read_csv(ext_load+'Luca_photom/C6en050918')
    Luca6en = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/K2.R7/C6_Andrea_27072018.in.mo',delimiter=r'\s+')
    # Luca6en = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/param_outputs/Poles/K2.R7/C6_Andrea_27072018.in.mo',delimiter=r'\s+')
    Luca6en = mdf.p_in('C6_Andrea_27072018',Luca6en,'C6')
    mdf.vert_dist(Luca6en)
    # Luca3 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_inputs/Poles/C3_Luca.in')
    # Luca6 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_inputs/Poles/C6_Luca.in')
    Lucaen = pd.concat([Luca3en,Luca6en],ignore_index=True)
    Lucaen = Lucaen.drop_duplicates(subset=['#Id']).reset_index(drop=True)
    # sys.exit()

    ''' Brief investigation of Teff relation between the EPIC and SkyMapper. EPIC Teff
        values are reasonably good, with a comparable scatter to that shown in the
        Huber et al. 2016 paper describing how the parameters were derived. '''
    # print(C3_New.columns.values)
    # print(C3_LucaIS.columns.values)
    # sys.exit()
    C3 = pd.merge(C3_New[['#Id','feh']],C3_LucaIS[['#Id','feh']],on=['#Id'])
    # ransac_fit(C3,C3,['feh_x','feh_y'],[r'$T_{\rm{eff}}$ - EPIC, C3',r'$T_{\rm{eff}}$ - SkyMapper'],0.1)


    C3_Luca['Glon'] = C3_Luca['GLON']
    C3_Luca['Glat'] = C3_Luca['GLAT']
    mdf.vert_dist(C3_Luca)
    C3_Luca['sig_Z'] = 1e-3 * abs(C3_Luca['dist_68U'] - C3_Luca['dist_68L'])/2
    C3_Luca['sig_Gal_Rad'] = np.sqrt((2*(C3_Luca['sig_Z']/C3_Luca['X'])**2) + (2*(C3_Luca['sig_Z']/C3_Luca['Y'])**2))

    C6_Luca['Glon'] = C6_Luca['GLON']
    C6_Luca['Glat'] = C6_Luca['GLAT']
    mdf.vert_dist(C6_Luca)
    C6_Luca['sig_Z'] = 1e-3 * abs(C6_Luca['dist_68U'] - C6_Luca['dist_68L'])/2
    C6_Luca['sig_Gal_Rad'] = np.sqrt((2*(C6_Luca['sig_Z']/C6_Luca['X'])**2) + (2*(C6_Luca['sig_Z']/C6_Luca['Y'])**2))

    C3_LucaIS['Glon'] = C3_LucaIS['GLON']
    C3_LucaIS['Glat'] = C3_LucaIS['GLAT']
    mdf.vert_dist(C3_LucaIS)
    C3_LucaIS['sig_Z'] = 1e-3 * abs(C3_LucaIS['dist_68U'] - C3_LucaIS['dist_68L'])/2
    C3_LucaIS['sig_Gal_Rad'] = np.sqrt((2*(C3_LucaIS['sig_Z']/C3_LucaIS['X'])**2) + (2*(C3_LucaIS['sig_Z']/C3_LucaIS['Y'])**2))

    C6_LucaIS['Glon'] = C6_LucaIS['GLON']
    C6_LucaIS['Glat'] = C6_LucaIS['GLAT']
    mdf.vert_dist(C6_LucaIS)
    C6_LucaIS['sig_Z'] = 1e-3 * abs(C6_LucaIS['dist_68U'] - C6_LucaIS['dist_68L'])/2
    C6_LucaIS['sig_Gal_Rad'] = np.sqrt((2*(C6_LucaIS['sig_Z']/C6_LucaIS['X'])**2) + (2*(C6_LucaIS['sig_Z']/C6_LucaIS['Y'])**2))

    ''' Data arrays '''
    data = [C3_New, C6_New, K2_New, K2_AS, AS, APK2, c_three, c_six, RC3, RC6, L3, L6, GES, AP3, AP6]
    data_en = [C3en, C6en, RC3en, RC6en, GESen, AP3en, AP6en, L3en, L6en]
    data_feh = [C3_feh, C6_feh]
    data_Luca = [C3_Luca, Luca3en, C6_Luca, C3_LucaIS, C6_LucaIS]

    # print(AS['alpha'].iloc[360:])
    # sys.exit()
    ''' Parameter cuts '''
    data = mdf.sig_cut(data, 'age', 0.5, 19.9)
    data_en = mdf.sig_cut(data_en, 'age', 0.5, 19.9)
    data_feh = mdf.sig_cut(data_feh, 'age', 0.5, 19.9)
    data_Luca = mdf.sig_cut(data_Luca, 'age', 0.5, 19.9)

    C3_New, C6_New, K2_New, K2_AS, AS, APK2, c_three, c_six, RC3, RC6, L3, L6, GES, AP3, AP6 = data
    C3en, C6en, RC3en, RC6en, GESen, AP3en, AP6en, L3en, L6en = data_en
    C3_feh, C6_feh = data_feh
    C3_Luca, Luca3en, C6_Luca, C3_LucaIS, C6_LucaIS = data_Luca
    Luca1 = pd.concat([C3_Luca, C6_Luca],ignore_index=True)

    ''' Gaia Radii Calculations '''
    # Luca = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_x_gaiadr2.csv')
    # cols = ['EPIC', 'logg', 'feh', 'teff', 'EBV']
    # Luca['EBV'] = 0
    # g2 = Luca[cols]
    # g2.to_csv('input.sample.all', header=None, sep=r' ', index=False)
    # p = subprocess.Popen(['./bcall'])
    # p.wait()
    # p = subprocess.Popen(["mv", "output.file.all", "/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_Gaia_BCs.csv"])
    # p.wait()
    # BC = pd.read_table('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_Gaia_BCs.csv', sep=r'\s+', header=0)
    # BC = BC[['ID', 'BC_1', 'BC_2']]
    # BC = BC.rename(columns={'ID':'EPIC','BC_1':'BC_K','BC_2':'BC_G'})
    # data_BC = pd.merge(Luca, BC, on=['EPIC'], how='inner')
    # data_BC.to_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_Gaia_BC_full.csv', index=False)
    # Luca = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_Gaia_BC_full.csv')
    # Gap3 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/GAP3')
    # Gap6 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/GAP6')
    # Gap = pd.concat([Gap3,Gap6],ignore_index=True)
    # Luca = pd.merge(Luca,Gap[['EPIC','Kmag']],how='inner',on=['EPIC'])
    # # print(Luca.columns.values)
    # sys.exit()
    # Luca['Kabs'] = Luca['Kmag'] - 5*np.log10(Luca['dist_ABJ']*100)
    # Luca['Kabs2'] = Luca['Kmag'] - 5*np.log10(Luca['dist_1P']*100)
    # Luca['A_K'] = Luca['Av']*0.31
    # T = (5777/Luca['teff'])**4 # Temperature term
    # Mm = (4.75-Luca['BC_K']-Luca['Kabs']+Luca['A_K'])/2.5 # Magnitude term including bolometric correction and extinction
    # Mm2 = (4.75-Luca['BC_K']-Luca['Kabs2']+Luca['A_K'])/2.5 # Magnitude term including bolometric correction and extinction
    # Luca['Rgaia'] = np.sqrt(T * 10**(Mm))
    # Luca['Rgaia2'] = np.sqrt(T * 10**(Mm2))
    # Luca['Lumo'] = Luca['rad']**2 * (Luca['teff']/const.solar_Teff)**4
    # Luca['gLumo'] = Luca['Rgaia']**2 * (Luca['teff']/const.solar_Teff)**4
    # Luca = Luca.dropna(subset=['Rgaia'])
    # Luca = Luca.reset_index(drop=True)
    # # print(Luca[['Rgaia','Kabs','Ks']],Mm,Luca.columns.values)
    # Luca.to_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_Gaia_BC_full.csv', index=False)
    # sys.exit()

    Luca = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_Gaia_BC_full.csv')
    # Luca = pd.read_csv('/home/ben/K2_Poles/Mass_Distr_In/Gaia/SM_Gaia_BC_full.csv')
    Luca = Luca.rename(columns={'EPIC':'#Id'})
    Luca = Luca[(Luca['age'] > 0.) & (Luca['age'] < 19.) & (Luca['sig_age']/Luca['age'] < 0.5)]
    Luca = Luca.rename(columns={'GLON':'Glon','GLAT':'Glat'})
    Luca = mdf.vert_dist(Luca)
    Luca['sig_Z'] = 1e-3 * abs(Luca['dist_68U'] - Luca['dist_68L'])/2
    Luca['sig_Gal_Rad'] = np.sqrt((2*(Luca['sig_Z']/Luca['X'])**2) + (2*(Luca['sig_Z']/Luca['Y'])**2))
    Luca = Luca.reset_index(drop=True)
    # print(np.mean(abs(Luca['dist_68U']-Luca['dist_68L'])/Luca['dist']))
    # sys.exit()
    Luca['sig_A_K'] = 0.15*(Luca['Av_68U']-Luca['Av_68L'])/2 # Conversion from Schaifers and Voigt 1982
    Luca['sig_BC'] = 0.03295*Luca['BC_K'] # From varying inputs to Luca's code. Teff had greatest effect at +/-100K limits
    Luca['sig_Kbol'] = (0.05*Luca['sig_dist_ABJ'])/(Luca['dist_ABJ']*np.log(10))
    Luca['sig_Rgaia'] = np.sqrt( ((2*Luca['eteff'])/Luca['teff'])**2 + ((-0.2*np.log(10))**2)*(Luca['sig_A_K']**2 + Luca['sig_BC']**2 + Luca['sig_Kbol']**2) )*Luca['Rgaia']
    Luca['sig_rad'] = abs(Luca['rad_68U']-Luca['rad_68L'])/2
    # print(Luca['sig_A_K'])
    # print(Luca['sig_BC'])
    # print(Luca['sig_Kbol'])
    # print((2*Luca['eteff'])/Luca['teff'])
    # print(np.median(Luca['sig_Rgaia']/Luca['Rgaia']))
    # print(np.median(Luca['sig_rad']/Luca['rad']))
    # print(len(Luca),len(C3_Luca),len(C6_Luca),len(Luca[Luca['Z']>0.]))
    # sys.exit()

    AS = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/AS_Gaia_BC_full.csv')
    # AS = pd.read_csv('/home/ben/K2_Poles/Mass_Distr_In/Gaia/AS_Gaia_BC_full.csv')
    AS_EPIC = pd.concat([RC3[['#Id']],RC6[['#Id']],GES[['#Id']],AP3[['#Id']],AP6[['#Id']]],ignore_index=True).drop_duplicates(subset=['#Id']).reset_index(drop=True)
    AS = pd.merge(AS,AS_EPIC,how='inner',on=['#Id'])
    AS = AS[(AS['age'] > 0.) & (AS['age'] < 19.) & (AS['sig_age']/AS['age'] < 0.5)]
    AS = mdf.vert_dist(AS)
    AS['sig_Z'] = 1e-3 * abs(AS['dist_68U'] - AS['dist_68L'])/2
    AS['sig_Gal_Rad'] = np.sqrt((2*(AS['sig_Z']/AS['X'])**2) + (2*(AS['sig_Z']/AS['Y'])**2))
    AS['APO']=0
    # for i in range(len(AS)):
    #     for j in range(len(AP3)):
    #         if AS['#Id'].iloc[i] == AP3['#Id'].iloc[j]:
    #             AS['APO'].iloc[i] = 1
    # for i in range(len(AS)):
    #     for j in range(len(AP6)):
    #         if AS['#Id'].iloc[i] == AP6['#Id'].iloc[j]:
    #             AS['APO'].iloc[i] = 1
    # AS = AS[AS['APO'] == 1]
    # AS.to_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/AS_Gaia_BC_APO.csv',index=False)
    # print(len(AS))
    # a = AS[AS['#Id'] < 208000000]
    # b = AS[AS['#Id'] > 208000000]
    # print(len(a),len(b))
    # sys.exit()

    # f = plt.figure()
    # plt.hist(AS['Rgaia'],bins=np.linspace(0,20,50),histtype='step',label=r'R$_{Gaia}$',normed=True,linewidth=2)
    # plt.hist(AS['rad'],bins=np.linspace(0,20,50),histtype='step',label=r'R$_{PARAM}$',normed=True,linewidth=2)
    # plt.xlim(3.5,18)
    # plt.xlabel(r'Radius [R$_{\odot}$]',fontsize=15)
    # plt.legend(prop={'size':10})
    # f = plt.figure()
    # plt.hist(AS['gLumo'],bins=50,histtype='step',label=r'R$_{Gaia}$',normed=True,linewidth=2)
    # plt.hist(AS['Lumo'],bins=50,histtype='step',label=r'R$_{PARAM}$',normed=True,linewidth=2)
    # # plt.xlim(3.5,18)
    # plt.xlabel(r'Luminosity [L$_{\odot}$]',fontsize=15)
    # plt.legend(prop={'size':10})
    # plt.show()
    # sys.exit()

    LucaIS = pd.concat([C3_LucaIS, C6_LucaIS],ignore_index=True)
    LucaIS = LucaIS[LucaIS['age'] > 0.]
    LucaIS = LucaIS[LucaIS['sig_age']/LucaIS['age'] < 0.4]
    LucaIS = LucaIS.reset_index(drop=True)

    ''' GOLDEN STANDARD K2 SAMPLE
    - Combined C3 and C6 sample
    - Age < 19 - Removal of boundary build up
    - Sig_Fractional_Age < 0.35 or 0.3
    - Fractional Delta(Mass) < (2 * median(sigma_mass)). Median calc'd prior to cuts
    - Gaia raduis/parallax flag to be introduced by Saniya
    '''
    C3_LucaIS['field'] = 3
    C6_LucaIS['field'] = 6
    gold = pd.concat([C3_LucaIS, C6_LucaIS],ignore_index=True)
    # print(len(gold[gold['field']==3]),len(gold[gold['field']==6]))
    med = np.median(gold['sig_mass'])
    gold = mdf.sig_cut([gold],'age',0.35,19.)
    gold = gold[0]
    # print(len(gold[gold['field']==3]),len(gold[gold['field']==6]))
    gold['mass_rat'] = abs(gold['m_seismo']-gold['mass'])/gold['mass']
    gold = gold[gold['mass_rat'] < 2.*med]
    gold = gold.reset_index(drop=True)
    ASg = pd.merge(AS,gold[['#Id']],how='inner',on=['#Id'])
    ASg = ASg.reset_index(drop=True)
    # print(ASg)
    # print(len(gold[gold['field']==3]),len(gold[gold['field']==6]))
    # gold.to_csv(ext_load+'GOLD2',index=False)
    # print(gold.columns.values)
    # sys.exit()

    # sns.set(style="whitegrid")
    # f, ax = plt.subplots()
    # sns.lineplot(data=gold,x="numax",y=0.22*gold['numax']**0.8)
    # sns.scatterplot(x="numax", y="Dnu", data=gold, ax=ax)
    # plt.show()
    # sys.exit()

    ''' Save out figures to a single pdf '''
    # pdf = matplotlib.backends.backend_pdf.PdfPages("K2_Gold1_KielMass.pdf")

    ''' File lengths '''
    # print(len(C3_New),len(C6_New))
    # print(len(c_three),len(c_six))
    # print(len(RC3),len(RC6))
    # print(len(GES))
    # print(len(L3),len(L6))
    # print(len(AP3),len(AP6))
    # print(len(K2_New),len(AS))
    # sys.exit()

    C3_New_c = C3_New[(C3_New['rad'] > 10.0) & (C3_New['rad'] < 12.0)]
    C6_New_c = C6_New[(C6_New['rad'] > 10.0) & (C6_New['rad'] < 12.0)]

    ''' Scaling Mass Cut '''
    mCut = 110.0
    C3_Luca = C3_Luca[abs(C3_Luca['percentage']) < mCut]
    C3_Luca = C3_Luca.reset_index(drop=True)
    C6_Luca = C6_Luca[abs(C6_Luca['percentage']) < mCut]
    C6_Luca = C6_Luca.reset_index(drop=True)
    c_three = c_three[abs(c_three['percentage']) < mCut]
    c_three = c_three.reset_index(drop=True)
    c_six = c_six[abs(c_six['percentage']) < mCut]
    c_six = c_six.reset_index(drop=True)

    # LucaAS = pd.merge(Luca[['#Id','teff','age']],AS[['#Id','Teff','age']],how='inner',on=['#Id'])
    # LucaAS['dT'] = LucaAS.Teff - LucaAS.teff
    # LucaAS['da'] = LucaAS.age_y - LucaAS.age_x
    # print(np.median(LucaAS['dT']), max(LucaAS['dT']),np.mean(LucaAS['dT']))
    # print(np.median(LucaAS['da']), max(LucaAS['da']),np.mean(LucaAS['da']))
    # sys.exit()

    Luca['m_seismo'] = (Luca['numax']/3090)**3 * (Luca['Dnu']/135.1)**-4 * (Luca['teff']/5777)**1.5
    Luca['percentage'] = 100*((Luca['m_seismo']/Luca['mass'])-1)
    AS['m_seismo'] = (AS['nmx']/3090)**3 * (AS['dnu']/135.1)**-4 * (AS['Teff']/5777)**1.5
    AS['percentage'] = 100*((AS['m_seismo']/AS['mass'])-1)
    APK2['m_seismo'] = (APK2['numax']/3090)**3 * (APK2['Dnu']/135.1)**-4 * (APK2['Teff']/5777)**1.5
    APK2['percentage'] = 100*((APK2['m_seismo']/APK2['mass'])-1)
    Luca = Luca[abs(Luca['percentage']) < mCut]
    Luca = Luca.reset_index(drop=True)
    AS = AS[abs(AS['percentage']) < mCut]
    AS = AS.reset_index(drop=True)
    APK2 = APK2[abs(APK2['percentage']) < mCut]
    APK2 = APK2.reset_index(drop=True)

    ''' Clump Cut '''
    # C3_Luca = C3_Luca[np.logical_or((C3_Luca['rad'] < 10.0) , (C3_Luca['rad'] > 11.5))] # Crude clump
    # C6_Luca = C6_Luca[np.logical_or((C6_Luca['rad'] < 10.0) , (C6_Luca['rad'] > 11.5))] # Crude clump
    # c_three = c_three[np.logical_or((c_three['rad'] < 10.0) , (c_three['rad'] > 11.5))] # Crude clump
    # c_six = c_six[np.logical_or((c_six['rad'] < 10.0) , (c_six['rad'] > 11.5))] # Crude clump
    # print(len(C3_Luca),len(C6_Luca),len(c_three),len(c_six))


    ''' Binning Campaign fields by Z '''
    z31 = C3_Luca[abs(C3_Luca['Z']) < 0.5]
    z32 = C3_Luca[(abs(C3_Luca['Z']) >= 0.5) & (abs(C3_Luca['Z'] < 1.0))]
    z33 = C3_Luca[(abs(C3_Luca['Z']) >= 1.0) & (abs(C3_Luca['Z'] < 1.5))]
    z34 = C3_Luca[(abs(C3_Luca['Z']) >= 1.5) & (abs(C3_Luca['Z'] < 2.0))]
    z35 = C3_Luca[(abs(C3_Luca['Z']) >= 2.0) & (abs(C3_Luca['Z'] < 2.5))]
    z36 = C3_Luca[abs(C3_Luca['Z']) > 2.5]

    z61 = C6_Luca[C6_Luca['Z'] < 0.5]
    z62 = C6_Luca[(C6_Luca['Z'] >= 0.5) & (C6_Luca['Z'] < 1.0)]
    z63 = C6_Luca[(C6_Luca['Z'] >= 1.0) & (C6_Luca['Z'] < 1.5)]
    z64 = C6_Luca[(C6_Luca['Z'] >= 1.5) & (C6_Luca['Z'] < 2.0)]
    z65 = C6_Luca[(C6_Luca['Z'] >= 2.0) & (C6_Luca['Z'] < 2.5)]
    z66 = C6_Luca[C6_Luca['Z'] > 2.5]

    # z = pd.read_csv('/home/bmr135/K2_Poles/Alt5_Zgt1.csv')
    # z = pd.merge(z,C3_Luca[['#Id','Z']],how='inner',on=['#Id'])
    # z['GES'] = 0
    # z['RAVE'] = 0
    # z['APOGEE'] = 0
    # for i in range(len(z)):
    #     for j in range(len(GES)):
    #         if z['#Id'].iloc[i] == GES['#Id'].iloc[j]:
    #             z['GES'].iloc[i] = 1
    # for i in range(len(z)):
    #     for j in range(len(RC3)):
    #         if z['#Id'].iloc[i] == RC3['#Id'].iloc[j]:
    #             z['RAVE'].iloc[i] = 1
    # for i in range(len(z)):
    #     for j in range(len(AP3)):
    #         if z['#Id'].iloc[i] == AP3['#Id'].iloc[j]:
    #             z['APOGEE'].iloc[i] = 1
    # # gaia = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Gaia/GAP_Gaia_DR2.csv')
    # # gz = pd.merge(z[['#Id']],gaia[['#Id','RA','Dec','radial_velocity','radial_velocity_error']],how='left',on=['#Id'])
    # print(z)
    # z.to_csv('Alt5_Zgt1.csv',index=False)
    # sys.exit()
    ''' Quick fire comps with Y enhanced grid results '''
    # C3en = C3en[C3en['mass'] > -90]
    # C3en.reset_index(drop=True)
    # ins = [C3_New, C6_New, C3en, C6en]
    # label = ['C3-Phot', 'C6-Phot', 'C3-Diff', 'C6-diff']
    #
    # for i in range(2):
    #     # p = 'age'
    #     # p1 = p+'_x'
    #     # p2 = p+'_y'
    #     # p3 = p+'_xy'
    #     # xtit1 = r'Age [Gyr], '
    #     # xtit2 = r'Phot Age [Gyr]'
    #     # xtit3 = r'Phot Age [Gyr]'
    #     # ytit2 = r'Phot-Diff Age [Gyr]'
    #     # ytit3 = r'$\Delta$A [(A$_{\rm{O}}$-A$_{\rm{Diff}}$)/A$_{\rm{O}}$]'
    #     # upper = 20
    #     # save = 'Age_'
    #
    #     # p = 'mass'
    #     # p1 = p+'_x'
    #     # p2 = p+'_y'
    #     # p3 = p+'_xy'
    #     # xtit1 = r'Mass [M$_{\odot}$], '
    #     # xtit2 = r'Phot Mass [M$_{\odot}$]'
    #     # xtit3 = r'Phot Mass [M$_{\odot}$]'
    #     # ytit2 = r'Phot-Diff Mass [M$_{\odot}$]'
    #     # ytit3 = r'$\Delta$M [(M$_{\rm{O}}$-M$_{\rm{Diff}}$)/M$_{\rm{O}}$]'
    #     # upper = 2.75
    #     # save = 'Mass_'
    #
    #     # p = 'rad'
    #     # p1 = p+'_x'
    #     # p2 = p+'_y'
    #     # p3 = p+'_xy'
    #     # xtit1 = r'Radius [R$_{\odot}$], '
    #     # xtit2 = r'Luca Radius [R$_{\odot}$]'
    #     # xtit3 = r'Luca Radius [R$_{\odot}$]'
    #     # ytit2 = r'Luca-Diff Radius [R$_{\odot}$]'
    #     # ytit3 = r'$\Delta$R [(R$_{\rm{O}}$-R$_{\rm{Diff}}$)/R$_{\rm{O}}$]'
    #     # upper = 20
    #     # save = 'Radius_'
    #
    #     a = int(len(ins)/2.)
    #
    #     fig, (ax1,ax2,ax3) = plt.subplots(3,sharex=True)
    #     d1 = kde.KDE1D(ins[i][p])
    #     d2 = kde.KDE1D(ins[i+a][p])
    #     x1 = np.r_[min(ins[i][p]):max(ins[i][p]):1024j]
    #     x2 = np.r_[min(ins[i+a][p]):max(ins[i+a][p]):1024j]
    #     ax1.plot(x1,d1(x1),linewidth=2,label=r'Phot.')
    #     ax1.plot(x2,d2(x2),linewidth=2,label=r'Phot Diff.')
    #     ax1.set_yticks([])
    #     ax1.set_xlim(0,upper)
    #     # ax1.set_xlabel(xtit1+label[i])
    #     ax1.legend()
    #     # plt.savefig('/home/bmr135/Dropbox/K2Poles/pop_trends/Luca_Tests/Full/KDE/'+save+label[i]+'.png')
    #
    #
    #     param = pd.DataFrame()
    #     param = pd.merge(ins[i][['#Id',p]],ins[i+a][['#Id',p]],how='inner',on=['#Id'])
    #     # plt.figure()
    #     ax2.plot([0,upper],[0,upper],color='k',linestyle='--')
    #     ax2.scatter(param[p1],param[p2],label=label[i])
    #     # ax2.set_xlabel(xtit2)
    #     ax2.set_ylabel(ytit2)
    #     ax2.legend()
    #     ax2.set_xlim(0,upper)
    #     ax2.set_ylim(0,upper)
    #     # plt.savefig('/home/bmr135/Dropbox/K2Poles/pop_trends/Luca_Tests/Full/1to1/'+save+label[i]+'.png')
    #
    #     # plt.figure()
    #     param[p3] = (param[p1] - param[p2])/param[p1]
    #     ax3.plot([0,upper],[0,0],color='k',linestyle='--')
    #     ax3.scatter(param[p1],param[p3],label=label[i])
    #     ax3.set_xlabel(xtit3)
    #     ax3.set_ylabel(ytit3)
    #     ax3.legend()
    #     ax3.set_xlim(0,upper)
    #
    #     plt.show()
    #     # plt.savefig('/home/bmr135/Dropbox/K2Poles/pop_trends/Luca_Tests/Full/delta_param/'+save+label[i]+'.png')
    #
    # sys.exit()

    # plt.figure()
    # param = pd.merge(C3_Luca[['#Id','feh']],c_three[['#Id','feh']],how='inner',on=['#Id'])
    # plt.scatter(param['feh_y'],param['feh_y']-param['feh_x'])
    # plt.show()

    ''' Age Distribution '''

    # fig, ((ax,ax1),(ax2,ax3),(ax4,ax5)) = plt.subplots(3,2,sharex='col',sharey=True)
    # ax.hist(z31['age'],bins=np.linspace(0,20,40),histtype='step',label=r'Z $<$ 0.5',linewidth=2)
    # ax.set_yticks([])
    # ax.legend()
    # ax1.hist(z32['age'],bins=np.linspace(0,20,40),histtype='step',label=r'0.5 $<$ Z $<$ 1.0',linewidth=2)
    # ax1.legend()
    # ax2.hist(z33['age'],bins=np.linspace(0,20,40),histtype='step',label=r'1.0 $<$ Z $<$ 1.5',linewidth=2)
    # ax2.set_yticks([])
    # ax2.legend()
    # ax3.hist(z34['age'],bins=np.linspace(0,20,40),histtype='step',label=r'1.5 $<$ Z $<$ 2.0',linewidth=2)
    # ax3.legend()
    # ax4.hist(z35['age'],bins=np.linspace(0,20,40),histtype='step',label=r'2.0 $<$ Z $<$ 2.5',linewidth=2)
    # ax4.set_yticks([])
    # ax4.legend()
    # ax5.hist(z36['age'],bins=np.linspace(0,20,40),histtype='step',label=r'Z $>$ 2.5',linewidth=2)
    # ax5.legend()
    # ax4.set_xlabel(r'Age [Gyr], C3$_{Luca}$ - Mass $< 2.25$')
    # ax5.set_xlabel(r'Age [Gyr], C3$_{Luca}$ - Mass $< 2.25$')
    # plt.tight_layout()
    # # plt.show()
    #
    # fig, ((ax,ax1),(ax2,ax3),(ax4,ax5)) = plt.subplots(3,2,sharex='col',sharey=True)
    # ax.hist(z61['age'],bins=np.linspace(0,20,40),histtype='step',label=r'Z $<$ 0.5',linewidth=2)
    # ax.set_yticks([])
    # ax.legend()
    # ax1.hist(z62['age'],bins=np.linspace(0,20,40),histtype='step',label=r'0.5 $<$ Z $<$ 1.0',linewidth=2)
    # ax1.legend()
    # ax2.hist(z63['age'],bins=np.linspace(0,20,40),histtype='step',label=r'1.0 $<$ Z $<$ 1.5',linewidth=2)
    # ax2.set_yticks([])
    # ax2.legend()
    # ax3.hist(z64['age'],bins=np.linspace(0,20,40),histtype='step',label=r'1.5 $<$ Z $<$ 2.0',linewidth=2)
    # ax3.legend()
    # ax4.hist(z65['age'],bins=np.linspace(0,20,40),histtype='step',label=r'2.0 $<$ Z $<$ 2.5',linewidth=2)
    # ax4.set_yticks([])
    # ax4.legend()
    # ax5.hist(z66['age'],bins=np.linspace(0,20,40),histtype='step',label=r'Z $>$ 2.5',linewidth=2)
    # ax5.legend()
    # ax4.set_xlabel(r'Age [Gyr], C6$_{Luca}$ - Mass $< 2.25$')
    # ax5.set_xlabel(r'Age [Gyr], C6$_{Luca}$ - Mass $< 2.25$')
    # plt.tight_layout()
    # plt.show()
    #
    # sys.exit()

    ''' Age Uncertainty Comparison - 1 to 1 plot '''
    # ages = pd.merge(AS[['#Id','age']],K2_AS[['#Id','age']],how='inner',on=['#Id'])
    # plt.figure()
    # plt.scatter(ages['age_x'],ages['age_y'])
    # plt.xlabel(r'Spectroscopic Age [Gyr]')
    # plt.ylabel(r'Photometric Age [Gyr]')
    # plt.show()

    # sys.exit()
    # # ax1.hist([10**RC6['logAge']/1e9,10**AP6['logAge']/1e9,10**L6['logAge']/1e9], \
    # #         bins=np.linspace(0,20,40),stacked=True,label=[r'RAVE',r'APOGEE',r'LAMOST'])
    # # ax1.hist(APK2['age'],bins=np.linspace(0,20,40),histtype='step',label=r'APOKASC',normed=True,linewidth=2)
    #

    ''' Age KDEs '''
#     d3 = kde.KDE1D(c_three['age'])
#     # d6 = kde.KDE1D(c_six['age'])
#     dA6 = kde.KDE1D(AP3['age'])
#     dR6 = kde.KDE1D(RC3['age'])
#     dL6 = kde.KDE1D(GES['age'])
#     x03 = np.r_[min(c_three['age']):max_age:1024j]
#     # x06 = np.r_[min(c_six['age']):max_age:1024j]
#     x1 = np.r_[min(AP3['age']):max(AP3['age']):1024j]
#     x2 = np.r_[min(RC3['age']):max(RC3['age']):1024j]
#     x3 = np.r_[min(GES['age']):max(GES['age']):1024j]
#     fig, ax1 = plt.subplots(1)
#     ax1.plot(x03,d3(x03),linewidth=2,label=r'C3 Spec.')
#     # ax1.plot(x06,d6(x06),linewidth=2,label=r'C6 Spec')
#     ax1.plot(x1,dA6(x1),linewidth=2,label=r'APOGEE')
#     ax1.plot(x2,dR6(x2),linewidth=2,label=r'RAVE')
#     ax1.plot(x3,dL6(x3),linewidth=2,label=r'GES')
#     ax1.legend()
# #
#     # fig, ax1 = plt.subplots(1)
#     # ax1.scatter(AP3['age'],(AP3['sig_age']/AP3['age']),label=r'APOGEE')
#     # ax1.scatter(RC3['age'],(RC3['sig_age']/RC3['age']),label=r'RAVE')
#     # ax1.scatter(GES['age'],(GES['sig_age']/GES['age']),label=r'Gaia-ESo')
#     ax1.set_xlabel(r'Age [Gyr], Spec.')
#     # ax1.set_ylabel(r'Fract. $\sigma_{\rm{Age}}$')
#     # ax1.legend()
#     plt.show()
#     sys.exit()
#     # K2_gd_age = AS[(AS['sig_age']/AS['age']) < 0.35]
#     # fig, ax1 = plt.subplots(1)
#     # ax1.hist(K2_gd_age['age'],bins=np.linspace(0,20,40),histtype='step',label=r'AS',normed=True,linewidth=2)
#     # ax1.set_xlabel(r'Age [Gyr]')
#     # ax1.legend()
#     # plt.show()

    ''' [alpha/Fe] vs Age '''
    '''
    - APOGEE: alpha
    - GES: alpha
    - LAMOST: alpha
    - RAVE: ALPHA
    '''
    # AP3 = AP3[AP3['alpha'] > -99.9]
    # AP6 = AP6[AP6['alpha'] > -99.9]
    # GES = GES[GES['alpha'] > -99.9]
    # RC3 = RC3[RC3['ALPHA'] > -99.9]
    # RC6 = RC6[RC6['ALPHA'] > -99.9]
    # # L6 = L6[L6['alpha'] > -99.9]
    # fig, ax = plt.subplots(1)
    # ax.scatter(AP3['age'],AP3['alpha'],label=r'APOGEE3')
    # ax.scatter(AP6['age'],AP6['alpha'],label=r'APOGEE6')
    # ax.scatter(RC3['age'],RC3['ALPHA'],label=r'RAVE3')
    # ax.scatter(RC6['age'],RC6['ALPHA'],label=r'RAVE6')
    # ax.scatter(GES['age'],GES['alpha'],label=r'GES3')
    # # ax.scatter(L6['age'],L6['alpha'],label=r'LAMOST6')
    # ax.set_xlabel(r'Age [Gyr]')
    # ax.set_ylabel(r'[$\alpha$/Fe]')
    # ax.legend()
    # plt.show()

    ''' Alpha vs Age Combined '''
    # fig, ax = plt.subplots(1)
    # hist2, xb2, yb2, im2 = plt.hist2d(AS['age'],AS['alpha'],bins=[np.linspace(0,20,41),np.linspace(-0.3,0.7,11)],cmap=colormaps.parula)#,normed=True)
    # cbar = plt.colorbar()
    # cbar.set_label(r'Number', rotation=270, fontsize=15, labelpad=25)
    # ax.set_ylabel(r'[$\alpha$/Fe]',fontsize=15, labelpad=20)
    # ax.set_xlabel(r'Age [Gyr]',fontsize=15, labelpad=10)
    # plt.tight_layout()
    #
    # fig, ax = plt.subplots(1)
    # ax.scatter(AS['age'],AS['alpha'])
    # ax.set_ylabel(r'[$\alpha$/Fe]',fontsize=15, labelpad=20)
    # ax.set_xlabel(r'Age [Gyr]',fontsize=15, labelpad=10)
    # plt.tight_layout()
    #
    # plt.show()
    # sys.exit()

    ''' Parameter distribution plots '''
    # K2_RGB = pd.DataFrame()
    APK2 = APK2[APK2['radius'] > 0.0]
    # APK2 = APK2[(APK2['rad'] < 9.5) | (APK2['rad'] > 12.0)] # Crude clump
    # K2_New = K2_New[(K2_New['rad'] < 9.5) | (K2_New['rad'] < 12.0)] # Crude clump
    # C3_Luca = C3_Luca[(C3_Luca['rad'] < 9.5) | (C3_Luca['rad'] > 12.0)]
    # C6_Luca = C6_Luca[(C6_Luca['rad'] < 9.5) | (C6_Luca['rad'] > 12.0)]
    # Luca  = Luca[(Luca['rad'] < 9.5) | (Luca['rad'] > 12.0)]
    # gold = gold[(gold['rad'] < 9.5) | (gold['rad'] > 12.0)]

    APK2 = APK2.reset_index(drop=True)
    K2_New = K2_New.reset_index(drop=True)
    C3_Luca = C3_Luca.reset_index(drop=True)
    C6_Luca = C6_Luca.reset_index(drop=True)
    Luca = Luca.reset_index(drop=True)
    gold = gold.reset_index(drop=True)
    # # sys.exit()
    # # AS = AS[AS['logAge'] > 9.9]
    # fig, axes = plt.subplots(3,2)
    # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
    # # plt.suptitle(r'log$_{10}$(Age) $> 10.0$', fontsize=15)
    # # plt.suptitle(r'$\forall$ R', fontsize=15)
    # ax0.hist(K2_AS['mass'],bins=np.linspace(0.5,2.5,20),histtype='step',label=r'Photom.',normed=True)
    # ax0.hist(AS['mass'],bins=np.linspace(0.5,2.5,20),histtype='step',label=r'Spec.',normed=True)
    # ax0.legend(prop={'size':10})
    # ax0.set_xlabel(r'Mass [M$_{\odot}$]')
    # ax0.set_yticks([])
    # ax1.hist(K2_AS['age'],bins=np.linspace(0,20,40),histtype='step',label=None,normed=True)
    # ax1.hist(AS['age'],bins=np.linspace(0,20,40),histtype='step',label=None,normed=True)
    # ax1.set_xlabel(r'Age [Gyr]')
    # ax1.set_yticks([])
    # ax2.hist(K2_AS['rad'],bins=np.linspace(3,20,50),histtype='step',label=None,normed=True)
    # ax2.hist(AS['rad'],bins=np.linspace(3,20,50),histtype='step',label=None,normed=True)
    # ax2.set_xlabel(r'Radius [R$_{\odot}$]')
    # ax2.set_yticks([])
    # ax3.hist(K2_AS['Z'],bins=np.linspace(-8,8,130),histtype='step',label=None,normed=True)
    # ax3.hist(AS['Z'],bins=np.linspace(-8,8,130),histtype='step',label=None,normed=True)
    # ax3.set_xlabel(r'Z [kpc]')
    # ax3.set_yticks([])
    # ax4.hist(K2_AS['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'K2',normed=1)
    # ax4.hist(AS['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=None,normed=1)
    # ax4.set_xlabel(r'[Fe/H]')
    # ax4.set_yticks([])
    # ax5.axis('off')
    # plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.show()
    # # if save_out == 1:
    # #     plt.savefig(ext_fig+folder_loc+'Kep_K2_age_distr.png')

    multi = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/multimodal.txt')
    # multi = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/param_outputs/Poles/multimodal.txt')
    Lmul = pd.merge(Luca,multi,how='inner',on=['#Id'])
    Lmul = Lmul.reset_index(drop=True)
    # sys.exit()

    AS['sig_rad'] = (AS['rad_68U']-AS['rad_68L'])/2
    Luca['sig_rad'] = (Luca['rad_68U']-Luca['rad_68L'])/2
    Lmul['sig_rad'] = (Lmul['rad_68U']-Lmul['rad_68L'])/2

    err = pd.merge(Luca[['#Id','rad','sig_rad','sig_age','age']],AS[['#Id','rad','sig_rad']],how='inner',on=['#Id'])

    # fig = plt.subplots()
    Luca['frac_age'] = (Luca['sig_age']/Luca['age'])
    indices = Luca['frac_age'] >= 0.35
    Lmul['frac_age'] = (Lmul['sig_age']/Lmul['age'])
    # indices = Lmul['frac_age'] >= 0.35
    # plt.errorbar(Luca['rad'],Luca['rad'],xerr=Luca['sig_rad'],yerr=Luca['sig_rad'],fmt='.',alpha=0.2,ecolor='grey')
    # plt.scatter(Luca['rad'],Luca['sig_rad'],c=Luca['frac_age'],cmap=colormaps.parula)#,alpha=0.2)
    # # plt.scatter(Lmul['rad'][~indices],Lmul['sig_rad'][~indices],c=Lmul['frac_age'][~indices],cmap=colormaps.parula)#,alpha=0.2)
    # # plt.scatter(err['rad_x'],err['rad_x'],c=err['frac_age'],cmap=colormaps.parula)
    # cbar = plt.colorbar()
    # # plt.scatter(Luca['rad'][indices],Luca['rad'][indices],color='r')#,alpha=0.2)
    # plt.scatter(Lmul['rad'],Lmul['frac_age'],color='r')
    # plt.xlim(6.5,9.5)
    # # plt.ylim(7,9)
    # plt.tight_layout()
    # # plt.show()

    # trunc = Luca[(Luca['rad'] > 6.5) & (Luca['rad'] < 9.5)]
    # trunc['frac_rad'] = (trunc['sig_rad']/trunc['rad'])
    # trunc[['#Id','age','age_68L','age_68U','frac_age','rad','rad_68L','rad_68U','frac_rad']].to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/bump.txt',index=False)
    # sys.exit()

    # print(Luca.columns.values)
    # sys.exit()

    # print(APK2['evstate'])
    # sys.exit()
    # APK2 = APK2[APK2['evstate'] == 2]


    ''' TRILEGAL Ks mag distributions '''
    Kep_Sim = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Standard_kepler_field/k1.6_K15.all.out.txt',delimiter=r'\s+')
    # Kep_Sim = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/Standard_kepler_field/k1.6_K15.all.out.txt',delimiter=r'\s+')
    mdf.sim_dist(Kep_Sim)
    ''' Application of K2 Selection Function to Kepler Simulation '''
    # Kep_Sim['Teff'] = 10**(Kep_Sim['logTe'])
    # Kep_Sim['g'] = 10**(Kep_Sim['logg'])
    # Kep_Sim['L'] = 10**(Kep_Sim['logL'])
    # Kep_Sim['radius'] = np.sqrt(Kep_Sim['Mass'] / (Kep_Sim['g']/4.441))
    # Kep_Sim['JK'] = Kep_Sim['Jmag'] - Kep_Sim['Kmag']
    # Kep_Sim['Vcut'] = Kep_Sim['Kmag'] + 2*(Kep_Sim['JK']+0.14) + 0.382*np.exp(2*(Kep_Sim['JK']-0.2))
    # Kep_Sim = prop.det_prob_Kepler(Kep_Sim,'numax',3090.0,135.1)
    # Kep_Sim = Kep_Sim[ (Kep_Sim['numax'] > 10) & (Kep_Sim['numax'] < 280) & \
    # (Kep_Sim['Vcut'] > 9) & (Kep_Sim['Vcut'] < 15) & (Kep_Sim['JK'] > 0.5)]
    # Kep_Sim['radius'] = 10**Kep_Sim['logR']
    # print(Kep_Sim.columns.values)
    # sys.exit()
    # trig, ax = plt.subplots()
    TRI3['Kabs'] = TRI3['Kmag'] - 5*np.log10(TRI3['dist']/10)
    TRI6['Kabs'] = TRI6['Kmag'] - 5*np.log10(TRI6['dist']/10)
    Kep_Sim['Kabs'] = Kep_Sim['Kmag'] - 5*np.log10(Kep_Sim['dist']/10)
    TRI3a = TRI3[TRI3['#Gc'] == 2]
    TRI3a = TRI3a.reset_index(drop=True)
    TRI3b = TRI3[TRI3['#Gc'] != 2]
    TRI3b = TRI3b.reset_index(drop=True)

    # TRI3['L'] = 10**TRI3['logL']
    # TRI6['L'] = 10**TRI6['logL']

    mesa = pd.read_csv('MESA_track_example',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])
    mesa2 = pd.read_csv('MESA_track_example2',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])
    mesa3 = pd.read_csv('MESA_track_example3',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])
    mesa4 = pd.read_csv('MESA_track_example4',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])
    mesa5 = pd.read_csv('MESA_track_example5',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])

    for i in range(len(mesa2)):
        if mesa2.age.iloc[i] < 6e8:
            mesa2.age.iloc[i] += 9.722177696000e9
    for i in range(len(mesa4)):
        if mesa4.age.iloc[i] < 1e9:
            mesa4.age.iloc[i] += 2.200204761600e+10

    mesa4['age'] = mesa4['age']*10**-9
    mesa2['age'] = mesa2['age']*10**-9
    hrd1, ax = plt.subplots(1,figsize=(10,10))
    ax.plot(10**mesa4['logTe'],mesa4['logL'],color='k',alpha=0.1,label='__nolegend__')
    ax.plot(10**mesa2['logTe'],mesa2['logL'],color='k',alpha=0.1,label='__nolegend__')
    ax.scatter(10**mesa4['logTe'],mesa4['logL'],c=mesa4['age'],cmap=colormaps.parula,s=50,edgecolor='grey',vmin=0., vmax=20,label=r'0.8 M$_{\odot}$; -0.25 dex')
    d = ax.scatter(10**mesa2['logTe'],mesa2['logL'],c=mesa2['age'],cmap=colormaps.parula,s=50,edgecolor='grey',vmin=0.,vmax=20,label=r'1.0 M$_{\odot}$; -0.25 dex',marker='D')
    cbar = hrd1.colorbar(d)
    cbar.set_label(r'Age [Gyr]', rotation=270, fontsize=15, labelpad=15)
    ax.set_xlabel(r'$T_{\rm{eff}}$ [K]',fontsize=15)
    ax.set_ylabel(r'log$_{10}$($L$/L$_{\odot}$)',fontsize=15)
    ax.invert_xaxis()
    ax.legend(loc=2)
    # hrd1.savefig('/home/bmr135/Dropbox/Happy_Place/Chapters/Intro/mass_time.pdf',bbox_inches='tight')
    # plt.show()
    # sys.exit()


    TRI3a = TRI3a.sort_values(['Teff'])
    TRI3b = TRI3b.sort_values(['Teff'])
    ''' TRILEGAL HRD with lines of constant radius '''
    hrd, ax = plt.subplots(1,figsize=(10,10))
    ax.plot(10**mesa5['logTe'],mesa5['logL'],color='k',label=r'0.8 M$_{\odot}$; -0.5 dex',linestyle='--')
    ax.plot(10**mesa['logTe'],mesa['logL'],color='k',label=r'1.0 M$_{\odot}$; -0.5 dex')
    ax.plot(10**mesa4['logTe'],mesa4['logL'],color='m',label=r'0.8 M$_{\odot}$; -0.25 dex',linestyle='--')
    ax.plot(10**mesa2['logTe'],mesa2['logL'],color='m',label=r'1.0 M$_{\odot}$; -0.25 dex')
    d = ax.scatter(TRI3['Teff'],TRI3['logL'],c=TRI3['M_H'],cmap=colormaps.parula,label=None,vmin=-2.0, vmax=1.0,alpha=0.3,s=10)
    cbar = hrd.colorbar(d)
    cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=15)
    cbar.set_clim(-2.0,1.0)
    cbar.ax.set_yticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0])#,update_ticks=True)
    cbar.ax.set_yticklabels(['-2.0','-1.5','-1.0','-0.5','0.0','0.5','1.0'])

    Te = np.linspace(min(TRI3['Teff']-200),max(TRI3['Teff']+50),2324)
    ax.plot(Te,np.log10(25*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'5 R$_{\odot}$')
    ax.text(0.03, 0.345, '5 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,fontsize=10)
    ax.plot(Te,np.log10(100*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'10 R$_{\odot}$')
    ax.text(0.965, 0.43, '10 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,fontsize=10)
    ax.plot(Te,np.log10(121*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'11 R$_{\odot}$')
    ax.text(0.965, 0.484, '11 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,fontsize=10)
    ax.plot(Te,np.log10(144*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'12 R$_{\odot}$')
    ax.text(0.965, 0.54, '12 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,fontsize=10)
    ax.plot(Te,np.log10(169*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'13 R$_{\odot}$')
    ax.text(0.965, 0.59, '13 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,fontsize=10)
    ax.plot(Te,np.log10(196*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'14 R$_{\odot}$')
    ax.text(0.965, 0.63, '14 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,fontsize=10)
    ax.plot(Te,np.log10(225*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'15 R$_{\odot}$')
    ax.text(0.965, 0.67, '15 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,fontsize=10)
    ax.plot(Te,np.log10(400*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'20 R$_{\odot}$')
    ax.text(0.95, 0.85, '20 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes,fontsize=10)
    ax.invert_xaxis()
    ax.set_ylim(0.8,np.log10(175))
    ax.set_xlim(max(TRI3['Teff']+50),min(TRI3['Teff']-200))
    ax.set_xlabel(r'$T_{\rm{eff}}$ [K]',fontsize=15)
    ax.set_ylabel(r'log$_{10}$($L$/L$_{\odot}$)',fontsize=15)

    ax2 = plt.axes([0,0,0.1,0.1])
    ip = InsetPosition(ax, [0.51,0.07,0.48,0.12])
    ax2.set_axes_locator(ip)
    ax2.invert_xaxis()
    TRI3a = TRI3a[TRI3a['label'] == 4]
    TRI3b = TRI3b[TRI3b['label'] == 4]
    ax2.scatter(TRI3b['Teff'],TRI3b['logL'],c=TRI3b['M_H'],cmap=colormaps.parula,label=None,vmin=-2.0, vmax=1.0,alpha=0.3,s=10)
    ax2.scatter(TRI3a['Teff'],TRI3a['logL'],c=TRI3a['M_H'],cmap=colormaps.parula,label=None,vmin=-2.0, vmax=1.0,alpha=0.3,s=10,edgecolor='k')
    ax2.plot(Te,np.log10(100*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'10 R$_{\odot}$')
    ax2.text(0.16, 0.91, '10 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes,fontsize=10)
    ax2.plot(Te,np.log10(121*(Te/5777)**4),alpha=0.5,color='k',linestyle='--',label=r'11 R$_{\odot}$')
    ax2.text(0.8, 0.06, '11 R$_{\odot}$', horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes,fontsize=10)
    ax2.set_xlim(5300,4450)
    ax2.set_ylim(1.65,1.85)
    plt.close()
    # hrd.savefig('clump_rads.pdf',bbox_inches='tight')
    # #
    # plt.show()
    # sys.exit()
    APK2_alpha = APK2[APK2['alpha'] > 0.1]
    # APK2_alpha = APK2_alpha[APK2_alpha['age'] < 7.]
    K2_alpha = AS[AS['alpha'] > 0.1]
    # print(len(APK2_alpha)/len(APK2))
    # sys.exit()
    # tri = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C3_self')
    # print(max(tri.logg))
    # sys.exit()
    # tri = tri[tri['logg'] > 3.5]
    # fig1, ax = plt.subplots()
    # d = ax.scatter(tri['Mass'],tri['logAge'],c=tri['M_H'],cmap=colormaps.parula,label=None,vmin=-2.0, vmax=1.0,alpha=0.3,s=10)
    # cbar = hrd.colorbar(d)
    # cbar.set_label(r'[M/H]', rotation=270, fontsize=15, labelpad=15)
    # cbar.set_clim(-2.0,1.0)
    # cbar.ax.set_yticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0])#,update_ticks=True)
    # cbar.ax.set_yticklabels(['-2.0','-1.5','-1.0','-0.5','0.0','0.5','1.0'])
    # ax.set_xlabel(r'Mass [M$_{\odot}$]',fontsize=15)
    # ax.set_ylabel(r'log$_{10}$(Age)',fontsize=15)
    # plt.show()

    # ax.hist(Kep_Sim['Kabs'],bins=np.linspace(-4.2,0.75,60),label=r'APOKASC',color='grey',alpha=0.2,normed=True)
    # ax.hist(TRI3['Kabs'],bins=np.linspace(-4.2,0.75,60),histtype='step',label=r'TRI C3',normed=True,linewidth=2)
    # ax.hist(TRI6['Kabs'],bins=np.linspace(-4.2,0.75,60),histtype='step',label=r'TRI C6',normed=True,linewidth=2)
    # # ax.hist(APK2['Ks'],bins=np.linspace(-4.2,0.75,60),histtype='step',label=r'APOKASC, full',normed=True,linewidth=2)
    # # ax.hist(APK2_alpha['Ks'],bins=np.linspace(-4.2,0.75,60),histtype='step',label=r'APOKASC, $\alpha$-rich',normed=True,linewidth=2)
    # plt.yticks([])
    # # plt.xlim(3.5,18)
    # plt.xlabel(r'$M_{K}$',fontsize=15)
    # ax.legend(prop={'size':10})
    # trig.savefig('Kabs_APOKASC_sel_func.pdf',bbox_inches='tight')
    # plt.show()
    # sys.exit()

    ''' [Fe/H] - K2 vs Kelper '''
    # fig, ax1 = plt.subplots(1)
    # ax1.hist(APK2['feh'],bins=np.linspace(-2,0.75,30),label=r'APOKASC',normed=True,linewidth=2,color='grey',alpha=0.2)
    # ax1.hist(AS['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'K2 Spec.',normed=True,linewidth=2)
    # # ax1.hist(gold['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'K2, Gold',normed=True,linewidth=2)
    # ax1.hist(Luca['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'K2 SM',normed=True,linewidth=2)
    # # ax1.hist(TRI3['M_H'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'TRILEGAL C3',normed=True,linewidth=2)
    # ax1.set_xlabel(r'[Fe/H]',fontsize=20)
    # ax1.set_yticks([])
    # ax1.tick_params(labelsize=15)
    # ax1.legend()
    # plt.tight_layout()
    # fig.savefig('feh.pdf', bbox_inches='tight')
    # # pdf.savefig(fig)
    # #
    # # fig, ax1 = plt.subplots(1)
    # # ax1.hist(AS['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'K2 Spec.',normed=True,linewidth=2)
    # # ax1.hist(APK2['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'APOKASC',normed=True,linewidth=2)
    # # # ax1.hist(K2_New['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'K2',normed=True,linewidth=2)
    # # ax1.set_xlabel(r'[Fe/H]')
    # # ax1.legend()
    # plt.show()
    # sys.exit()

    ''' Percentage Difference Between PARAM and Scaling Relations '''
    # df = pd.DataFrame()
    # df['percentage'] = 100*((C3_Luca['m_seismo']/C3_Luca['mass'])-1)
    # plt.figure()
    # plt.scatter(Luca3en['mass'],Luca3en['percentage'])
    # plt.plot([0.5,max(Luca3en['mass'])+0.1],[0,0],linewidth=2,color='k')
    # plt.xlabel(r'Mass$_{\rm{PARAM}}$, C3 Luca Diff')
    # plt.xlim([0.5,max(Luca3en['mass'])+0.1])
    # plt.ylabel(r'$\%$ Difference between PARAM and Scaling-Relations')
    # plt.savefig(ext_fig+'Mass_percent_diff_spectro')

    # df = pd.DataFrame()
    # df['percentage'] = 100*((C6_Luca['m_seismo']/C6_Luca['mass'])-1)
    # plt.figure()
    # plt.scatter(C6_Luca['mass'],C6_Luca['percentage'])
    # plt.plot([0.5,max(C6_Luca['mass'])+0.1],[0,0],linewidth=2,color='k')
    # plt.xlabel(r'Mass$_{\rm{PARAM}}$, C6 Luca')
    # plt.xlim([0.5,max(C6_Luca['mass'])+0.1])
    # plt.ylabel(r'$\%$ Difference between PARAM and Scaling-Relations')
    # plt.savefig(ext_fig+'Mass_percent_diff_spectro')
    # plt.show()
    #
    # sys.exit()

    ''' Kepler Simulation if required '''
    # # Kep_Sim = pd.read_csv('/home/bmr135/GA/K2Poles/Standard_kepler_field/k1.6_K15.all.out.txt',delimiter=r'\s+')
    # Kep_Sim = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/Standard_kepler_field/k1.6_K15.all.out.txt',delimiter=r'\s+')
    # ''' Application of K2 Selection Function to Kepler Simulation '''
    # Kep_Sim['Teff'] = 10**(Kep_Sim['logTe'])
    # Kep_Sim['g'] = 10**(Kep_Sim['logg'])
    # Kep_Sim['L'] = 10**(Kep_Sim['logL'])
    # Kep_Sim['radius'] = np.sqrt(Kep_Sim['Mass'] / (Kep_Sim['g']/4.441))
    # Kep_Sim['JK'] = Kep_Sim['Jmag'] - Kep_Sim['Kmag']
    # Kep_Sim['Vcut'] = Kep_Sim['Kmag'] + 2*(Kep_Sim['JK']+0.14) + 0.382*np.exp(2*(Kep_Sim['JK']-0.2))
    # # Kep_Sim = prop.det_prob_Kepler(Kep_Sim,'numax',3090.0,135.1)
    # # Kep_Sim = Kep_Sim[ (Kep_Sim['numax'] > 10) & (Kep_Sim['numax'] < 280) & \
    # # (Kep_Sim['Vcut'] > 9) & (Kep_Sim['Vcut'] < 15) & (Kep_Sim['JK'] > 0.5)]


    ''' [Fe/H] vs [alpha/Fe] '''
    # mdf.met_comp([RC3,RC6,GES,APO],[r'RC3',r'RC6',r'GES',r'APOGEE'],'feh','alpha',['r','g','b','m'],[r'[Fe/H]',r'[$\alpha$/Fe]'])
    # plt.show()
    # mdf.met_comp([RC3],[r'RC3'],'feh','alpha',['r'],[r'[Fe/H]',r'[$\alpha$/Fe]'])
    # mdf.plt.show()
    # mdf.met_comp([RC6],[r'RC6'],'feh','alpha',['g'],[r'[Fe/H]',r'[$\alpha$/Fe]'])
    # plt.show()
    # mdf.met_comp([APO],[r'APOGEE'],'feh','alpha',['m'],[r'[Fe/H]',r'[$\alpha$/Fe]'])
    # plt.show()
    # mdf.met_comp([c_three],[r'C3 Spectro'],'feh','alpha',['y'],[r'[Fe/H]',r'[$\alpha$/Fe]'])
    # plt.show()
    # mdf.met_comp([GES],[r'GES'],'feh','alpha',['b'],[r'[Fe/H]',r'[$\alpha$/Fe]'])
    # plt.show()
    # mdf.met_comp([c_six],[r'C6 Spectro'],'feh','alpha',['orange'],[r'[Fe/H]',r'[$\alpha$/Fe]'])
    # plt.show()

    ''' Mulitple C3/C6 plots '''
    # param = ['feh','alpha','age']
    # label = [r'[Fe/H]',r'[$\alpha$/Fe]',r'Age [Gyr]']
    # for i in range (0,len(param),1):
    #     mdf.c3_6_param(AS,param[i],label[i])lt.subplots()
    # a = np.where(K2_camp_v2['prob_s'] >= 0.95)
    # b = np.where(K2_camp_v3['prob_s_gaia'] >= 0.95)
    # ax.scatter(K2_camp['JK'],K2_camp['Kabs'],alpha=0.5,label=r'GAP')
    # # ax.scatter(K2_camp_v2['JK'].iloc[a],K2_camp_v2['Kabs'].iloc[a],alpha=0.5,label=r'R$_{\rm{catalogue}}$')
    # ax.scatter(K2_camp_v3['JK'].iloc[b],K2_camp_v3['Kabs'].iloc[b],alpha=0.5,label=r'R$_{\rm{Gaia}}$')
    # ax.set_xlabel(r'J - K',fontsize=15)
    # ax.set_ylabel(r'K$_{\rm{abs}}$',fontsize=15)
    # ax.set_xlim(0.475,1.325)
    # ax.invert_yaxis()
    # ax.legend(loc=4)
    # plt.show()
    # fig.savefig('Det_prob_cut_HRDv2.pdf', bbox_inches='tight')
    # sys.exit()
    # plt.show()

    ''' K2 Age and Z Plots '''
    # plt.figure()
    # for j in range(itrs):
    #     if samp == 1:
    #         TRI = TRI.sample(n=len(K2))
    #         TRI = TRI.reset_index(drop=True)
    #     sim = copy.deepcopy(TRI)
    #     obs = copy.deepcopy(K2)
    #     alt_sim = mdf.sim_pert(obs,sim)
    #     # print(alt_sim.columns.values)
    #     thin_alt = alt_sim[alt_sim['#Gc'] == 1]
    #     thick_alt = alt_sim[alt_sim['#Gc'] == 2]
    #     param = ['logAge','Z']
    #     # param = ['age','Z']
    #     label = [r'$\rm{log}_{10}(\rm{Age})$',r'Z [kpc]']
    #     tag = [r'Tri Perturbed',r'K2 Observations']
    #     log=[0,0]
    #     ranges = [np.linspace(8.5,10.5,75),np.linspace(-8,8,130)]
    #     for i in range(2):
    #         plt.subplot(2,1,i+1)
    #         mdf.histo(K2,param[i],ranges[i],label[i],log[i],tag[1])
    #         mdf.histo(alt_sim,param[i],ranges[i],label[i],log[i],tag[0])
    #         plt.legend(prop={'size':15},loc=2)
    #         if clump == 1:
    #             plt.title(r'R $< 9$')
    # plt.show()
    #
    # # plt.title(r'$\forall$ R',fontsize=20)
    # # plt.tight_layout()
    # # plt.savefig(ext_fig+'K2_age_100_itr'+sys.argv[1]+'.png')
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'K2_age_Z.png')

    ''' Kepler Age and Z Plots '''
    # APK2['logAge'] = np.log10(APK2['age']*10**9)
    # plt.figure()
    # for j in range(itrs):
    #     if samp == 1:
    #         Kep_Sim = Kep_Sim.sample(n=len(APK2))
    #         Kep_Sim = Kep_Sim.reset_index(drop=True)
    #     sim_kep = copy.deepcopy(Kep_Sim)
    #     obs_kep = copy.deepcopy(APK2)
    #     alt_sim_kep = mdf.sim_pert(obs_kep,sim_kep)
    #     thin_alt_kep = alt_sim_kep[alt_sim_kep['#Gc'] == 1]
    #     thick_alt_kep = alt_sim_kep[alt_sim_kep['#Gc'] == 2]
    #     param_kep = ['logAge','Z']
    #     # param_kep = ['age','Z']
    #     label_kep = [r'$\rm{log}_{10}(\rm{Age})$',r'Z [kpc]']
    #     tag_kep = [r'Kepler Sim Perturbed',r'APOKASC Observations']
    #     log=[0,0]
    #     ranges_kep = [np.linspace(8.5,10.5,100),np.linspace(0,8,65)]
    #     for i in range(1):
    #         plt.subplot(1,1,i+1)
    #         mdf.histo(APK2,param_kep[i],ranges_kep[i],label_kep[i],log[i],tag_kep[1])
    #         mdf.histo(K2,param[i],ranges[i],label[i],log[i],tag[1])
    #         # mdf.histo(alt_sim_kep,param_kep[i],ranges_kep[i],label_kep[i],log[i],tag_kep[0])
    #         plt.legend(prop={'size':15},loc=2)
    #         if clump == 1:
    #             plt.title(r'R $< 9$')
    # # plt.savefig(ext_fig+'Kepler_APOKASC_age_100_itr'+sys.argv[1]+'.png')
    # plt.show()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'Kepler_APOKASC_age_Z.png')

    # thin = TRI[TRI['#Gc'] == 1]
    # thick = TRI[TRI['#Gc'] == 2]

    ''' Thick/Thin Age Comparison Histograms '''
    # plt.figure()
    # plt.hist([thick_alt['logAge'],thin_alt['logAge']],bins=ranges[0],stacked=True,label=[r'K2 Sim Thick',r'K2 Sim Thin'],normed=True)
    # plt.hist(K2['logAge'][np.isfinite(K2['logAge'])],bins=ranges[0],histtype='step',color='k',label=r'K2 Observations',normed=True)
    # plt.xlabel(r'log$_{10}$(Age)')
    # plt.title(r'$\forall$ R')
    # # plt.title(r'R $< 9$')
    # plt.legend(prop={'size':15},loc=2)
    # plt.show()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'K2_stacked_age.png')

    ''' PARAM vs Simulated vs Scaling Relation Mass '''
    # plt.figure()
    # bins = np.linspace(0.5,2.5,20)
    # mdf.histo(gold,'mass',bins,r'Mass [M$_{\odot}$]',0,r'K2 PARAM')
    # # mdf.histo(alt_sim,'Mass',bins,r'Mass [M$_{\odot}$]',0,r'TRILEGAL')
    # mdf.histo(gold,'m_seismo',bins,r'Mass [M$_{\odot}$]',0,r'K2 S-R')
    # plt.legend(prop={'size':15})
    # plt.show()
    # sys.exit()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'Mass_comp_param_sim_SR.png')

    ''' Spectro tests '''
    ges,c3 = mdf.spectro(GES,C3_New)
    ap3,c3a = mdf.spectro(AP3,C3_New)
    ap6,c6a = mdf.spectro(AP6,C6_New)
    rc3,c3_R = mdf.spectro(RC3,C3_New)
    rc6,c6_R = mdf.spectro(RC6,C6_New)
    l3, c3_l = mdf.spectro(L3,C3_New)
    l6, c6_l = mdf.spectro(L6,C6_New)
    GR3,RG3 = mdf.spectro(GES,RC3)
    AR6,RA6 = mdf.spectro(AP6,RC6)
    RL6,LR6 = mdf.spectro(RC6,L6)
    AL6,LA6 = mdf.spectro(AP6,L6)
    AR3,RA3 = mdf.spectro(AP3,RC3)
    AG3,GA3 = mdf.spectro(AP3,GES)
    # print(len(AR3),len(GR3), len(AG3))
    # print(len(AR6),len(RL6),len(AL6))

    ''' Teff/Met comparisons '''
    # # [Fe/H]
    # # ds = mdf.least_squares(c3,ges)
    # # ds1 = mdf.least_squares(c3_R,rc3)
    # # ds2 = mdf.least_squares(c6_R,rc6)
    # # ds3 = mdf.least_squares(c6a,ap6)
    # # ds3a = mdf.least_squares(c3a,ap3)
    # # ds4 = mdf.least_squares(c6_l,l6)
    #
    # df = mdf.least_squares(GR3,RG3)
    # df1 = mdf.least_squares(AR6,RA6)
    # df2 = mdf.least_squares(AR3,RA3)
    # df3 = mdf.least_squares(AG3,GA3)
    # df4 = mdf.least_squares(AL6,LA6)
    # df5 = mdf.least_squares(RL6,LR6)
    #
    # # Teff
    # # dsT = mdf.least_squares2(c3,ges)
    # # ds1T = mdf.least_squares2(c3_R,rc3)
    # # ds2T = mdf.least_squares2(c6_R,rc6)
    # # ds3T = mdf.least_squares2(c6a,ap6)
    # # ds3Ta = mdf.least_squares2(c3a,ap3)
    # # ds4T = mdf.least_squares2(c6_l,l6)
    #
    # dfT = mdf.least_squares2(GR3,RG3)
    # df1T = mdf.least_squares2(AR6,RA6)
    # df2T = mdf.least_squares2(AR3,RA3)
    # df3T = mdf.least_squares2(AG3,GA3)
    # df4T = mdf.least_squares2(AL6,LA6)
    # df5T = mdf.least_squares2(RL6,LR6)
    #
    # # x = np.linspace(4100,5500,100)
    # # plt.figure()
    # # plt.scatter(c3['Teff'],ges['Teff'],label=r'GES C3')
    # # plt.scatter(c3_R['Teff'],rc3['Teff'],label=r'RAVE C3',color='g')
    # # plt.scatter(c6_R['Teff'],rc6['Teff'],label=r'RAVE C6',color='r')
    # # plt.scatter(c6['Teff'],apo['Teff'],label=r'APOGEE C6',color='m')
    # # plt.xlabel(r'MAST T$_{\rm{eff}}$')
    # # plt.ylabel(r'Spectroscopic T$_{\rm{eff}}$')
    # # plt.plot([4100,5500],[4100,5500],c='k')
    # # plt.plot(x,(x*ds[0])+ds[1],linewidth=2,c='b')
    # # plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2,c='g')
    # # plt.plot(x,(x*ds2[0])+ds2[1],linewidth=2,c='r')
    # # plt.plot(x,(x*ds3[0])+ds3[1],linewidth=2,c='m')
    # # plt.xlim(4100,5500)
    # # plt.ylim(4100,5500)
    # # plt.legend(loc=4)
    # # plt.show()
    # #
    # # plt.figure()
    # # plt.subplot(2,2,1)
    # # plt.scatter(c3['feh'],ges['feh'],label=r'GES C3')
    # # plt.scatter(c3_R['feh'],rc3['feh'],label=r'RAVE C3',color='g')
    # # plt.scatter(c6_R['feh'],rc6['feh'],label=r'RAVE C6',color='r')
    # # plt.scatter(c6['feh'],apo['feh'],label=r'APOGEE C6',color='m')
    # # plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    # # plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    # # plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    # # plt.plot(x,(x*ds[0])+ds[1],linewidth=2,c='b')
    # # plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2,c='g')
    # # plt.plot(x,(x*ds2[0])+ds2[1],linewidth=2,c='r')
    # # plt.plot(x,(x*ds3[0])+ds3[1],linewidth=2,c='m')
    # # plt.xlim(-3.0,1.0)
    # # plt.ylim(-3.0,1.0)
    # # plt.tick_params(labelsize=15)
    # # plt.title(r'All Pipelines',fontsize=20)
    # #
    # # plt.subplot(2,2,2)
    # # plt.scatter(c3['feh'],ges['feh'],label=r'C3')
    # # plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    # # plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    # # plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    # # plt.plot(x,(x*ds[0])+ds[1],linewidth=2,c='b')
    # # plt.xlim(-3.0,1.0)
    # # plt.ylim(-3.0,1.0)
    # # plt.text(-.1, -2.75,r'Fit: %.6sx $+$ %.6s' %(ds[0][0],ds[1][0]), ha='center', va='center',fontsize=15)
    # # plt.tick_params(labelsize=15)
    # # plt.title(r'Gaia-ESO',fontsize=20)
    # #
    # # plt.subplot(2,2,3)
    # # plt.scatter(c3_R['feh'],rc3['feh'],label=r'C3',color='g')
    # # plt.scatter(c6_R['feh'],rc6['feh'],label=r'C6',color='r')
    # # plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    # # plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    # # plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    # # plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2,c='g')
    # # plt.plot(x,(x*ds2[0])+ds2[1],linewidth=2,c='r')
    # # plt.xlim(-3.0,1.0)
    # # plt.ylim(-3.0,1.0)
    # # plt.text(-.1, -2.5,r'Fit RC3: %.6sx $+$ %.6s' %(ds1[0][0],ds1[1][0]), ha='center', va='center',fontsize=15)
    # # plt.text(-.1, -2.75,r'Fit RC6: %.6sx $+$ %.6s' %(ds2[0][0],ds2[1][0]), ha='center', va='center',fontsize=15)
    # # plt.tick_params(labelsize=15)
    # # plt.title(r'RAVE',fontsize=20)
    # # plt.legend()
    # #
    # # plt.subplot(2,2,4)
    # # plt.scatter(c6['feh'],apo['feh'],label=r'C6',color='m')
    # # plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    # # plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    # # plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    # # plt.plot(x,(x*ds3[0])+ds3[1],linewidth=2,c='m')
    # # plt.xlim(-3.0,1.0)
    # # plt.ylim(-3.0,1.0)
    # # plt.tick_params(labelsize=15)
    # # plt.title(r'APOGEE',fontsize=20)
    # # plt.text(-.1, -2.75,r'Fit: %.6sx $+$ %.6s' %(ds3[0][0],ds3[1][0]), ha='center', va='center',fontsize=15)
    # #
    # # plt.tight_layout()
    # # plt.show()
    #
    # plt.figure()
    # plt.subplot(2,1,1)
    # plt.scatter(GR3['Teff'],RG3['Teff'])
    # plt.xlabel(r'Gaia-ESO T$_{\rm{eff}}$')
    # plt.ylabel(r'RAVE T$_{\rm{eff}}$')
    # plt.plot([4100,5500],[4100,5500],c='k')
    # # plt.plot(x,(x*df2[0])+df2[1],linewidth=2)
    # plt.xlim(4100,5500)
    # plt.title(r'C3')
    # plt.subplot(2,1,2)
    # plt.scatter(AR6['Teff'],RA6['Teff'])
    # plt.xlabel(r'APOGEE T$_{\rm{eff}}$')
    # plt.ylabel(r'RAVE T$_{\rm{eff}}$')
    # plt.plot([4100,5500],[4100,5500],c='k')
    # # plt.plot(x,(x*df3[0])+df3[1],linewidth=2)
    # plt.xlim(4100,5500)
    # plt.title(r'C6')
    #
    # plt.tight_layout()
    # plt.show()

    ''' Uncertainty asymmetry plots '''
    # plt.figure()
    # plt.errorbar((C3_Luca['age'])-1,3*C3_Luca.index.values,xerr=C3_Luca['sig_age'],fmt='o')
    # # plt.plot([-0.35,-0.35],[0,3*max(C3_Luca.index.values)],color='m')
    # plt.plot([13.8,13.8],[0,3*max(C3_Luca.index.values)],color='m',linewidth=2)
    # plt.show()
    # sys.exit()

    ''' Spectroscopic and Photometric Parameter Comparison Plots - FeH/Teff '''
    # x = np.linspace(4100,5500,1000)
    # plt.figure()
    # plt.subplot(3,1,1)
    # plt.scatter(AR3['Teff'],RA3['Teff']-AR3['Teff'])
    # plt.xlabel(r'T$_{\rm{eff}}$')
    # plt.ylabel(r'$\Delta$T$_{\rm{eff}}$ (RAVE - APOGEE)')
    # plt.plot([4100,5500],[0,0],c='k')
    # plt.plot(x,(x*df2T[0])+df2T[1],linewidth=2,color='r',linestyle='--')
    # plt.xlim(4100,5500)
    # plt.xticks([])
    # plt.title(r'C3')
    # plt.text(5300, -100,r'Offset: %.4sK' %(df2T[1][0]), ha='center', va='center')
    #
    # plt.subplot(3,1,2)
    # plt.scatter(GA3['Teff'],GA3['Teff']-AG3['Teff'])
    # # plt.xlabel(r'T$_{\rm{eff}}$')
    # plt.ylabel(r'$\Delta$T$_{\rm{eff}}$ (Gaia-ESO - APOGEE)')
    # plt.plot([4100,5500],[0,0],c='k')
    # plt.plot(x,(x*df3T[0])+df3T[1],linewidth=2,color='r',linestyle='--')
    # plt.xlim(4100,5500)
    # plt.xticks([])
    # plt.text(5300, -200,r'Offset: %.4sK' %(df3T[1][0]), ha='center', va='center')
    #
    # plt.subplot(3,1,3)
    # plt.scatter(RG3['Teff'],RG3['Teff']-GR3['Teff'])
    # plt.xlabel(r'T$_{\rm{eff}}$')
    # plt.ylabel(r'$\Delta$T$_{\rm{eff}}$ (RAVE - Gaia-ESO)')
    # plt.plot([4100,5500],[0,0],c='k')
    # plt.plot(x,(x*dfT[0])+dfT[1],linewidth=2,color='r',linestyle='--')
    # plt.xlim(4100,5500)
    # plt.text(5300, -200,r'Offset: %.4sK' %(dfT[1][0]), ha='center', va='center')
    #
    # # plt.tight_layout()
    # plt.show()
    #
    #
    # x = np.linspace(-1.5,.5,100)
    # plt.figure()
    # plt.subplot(3,1,1)
    # plt.scatter(RA3['feh'],RA3['feh']-AR3['feh'])
    # plt.xlabel(r'[Fe/H]')
    # plt.ylabel(r'$\Delta$[Fe/H] (RAVE - APOGEE)')
    # plt.plot([-1.5,0.5],[0,0],c='k')
    # plt.plot(x,(x*df2[0])+df2[1],linewidth=2,color='r',linestyle='--')
    # plt.xlim(-1.5,0.5)
    # plt.xticks([])
    # plt.title(r'C6')
    # plt.text(-1.25, -0.5,r'Offset: %.5sdex' %(df2[1][0]), ha='center', va='center')
    #
    # plt.subplot(3,1,2)
    # plt.scatter(GA3['feh'],GA3['feh']-AG3['feh'])
    # plt.xlabel(r'[Fe/H]')
    # plt.ylabel(r'$\Delta$[Fe/H] (Gaia-ESO - APOGEE)')
    # plt.plot([-1.5,0.5],[0,0],c='k')
    # plt.plot(x,(x*df3[0])+df3[1],linewidth=2,color='r',linestyle='--')
    # plt.xticks([])
    # plt.xlim(-1.5,0.5)
    # # plt.title(r'C3')
    # plt.text(-1.25, -0.5,r'Offset: %.5sdex' %(df3[1][0]), ha='center', va='center')
    #
    # plt.subplot(3,1,3)
    # plt.scatter(RG3['feh'],RG3['feh']-GR3['feh'])
    # plt.xlabel(r'[Fe/H]')
    # plt.ylabel(r'$\Delta$[Fe/H] (RAVE - Gaia-ESO)')
    # plt.plot([-1.5,0.5],[0,0],c='k')
    # plt.plot(x,(x*df[0])+df[1],linewidth=2,color='r',linestyle='--')
    # plt.xlim(-1.5,0.5)
    # # plt.title(r'C3')
    # plt.text(-1.25, -0.5,r'Offset: %.5sdex' %(df[1][0]), ha='center', va='center')
    #
    # # plt.tight_layout()
    # plt.show()

    ''' Metallicity offsets: EPIC to Spec. '''
    # EPIC = pd.DataFrame()
    # EPIC['#Id'] = K2_AS['#Id']
    # EPIC = pd.merge(EPIC,AS[['#Id']],how='inner',on=['#Id'])
    # AS1 = pd.merge(AS,EPIC,how='inner',on=['#Id'])
    # AS1 = AS1.reset_index(drop=True)
    # K2_AS = pd.merge(K2_AS,EPIC,how='inner',on=['#Id'])
    # K2_AS = K2_AS.reset_index(drop=True)
    #
    # K2_AS['delta_met'] = K2_AS['feh'] - AS1['feh']
    # K2_AS['delta_teff'] = K2_AS['Teff'] - AS1['Teff']
    #
    # EPIC1 = pd.DataFrame()
    # EPIC1['#Id'] = C3_Sky['#Id']
    # EPIC1 = pd.merge(EPIC1,AS[['#Id']],how='inner',on=['#Id'])
    # AS_Sky = pd.merge(AS,EPIC1,how='inner',on=['#Id'])
    # AS_Sky = AS_Sky.reset_index(drop=True)
    # C3_Sky = pd.merge(C3_Sky,AS_Sky[['#Id','feh']],how='inner',on=['#Id'])
    # C3_Sky = C3_Sky.reset_index(drop=True)
    #
    # C3_Sky['delta_met'] = C3_Sky['[Fe/H]'] - C3_Sky['feh']
    #
    # bins = np.linspace(-4.,4.5,18)
    # dfeh, edges, number = scipy.stats.binned_statistic(K2_AS['Z'],K2_AS['delta_met'],statistic='median',bins=bins)
    # stdFeh, stdEdges, stdNumber = scipy.stats.binned_statistic(K2_AS['Z'],K2_AS['delta_met'],statistic=lambda x: np.std(x),bins=bins)
    # # print(dfeh,edges)
    # # print(stdFeh,stdEdges)
    # # print(np.std(K2_AS['delta_met']),np.mean(K2_AS['delta_met']))
    # fig, ax = plt.subplots()
    # # ax.fill_between(edges[1:]-.25,dfeh-stdFeh,dfeh+stdFeh,color='orange',alpha=0.2)
    # a = ax.scatter(K2_AS['J']-K2_AS['Ks'],K2_AS['mbol'],c=K2_AS['delta_met'])
    # ax.plot(K2_AS['J']-K2_AS['Ks'],-13*(K2_AS['J']-K2_AS['Ks'])+9.,'g')
    # plt.gca().invert_yaxis()
    # plt.xlim(0.4,0.8)
    # plt.ylim(3,-1)
    # # ax.plot(edges[1:]-0.25,dfeh,color='orange',linewidth=2.)
    # cbar = fig.colorbar(a)
    # # cbar.set_label(r'Age [Gyr]', rotation=270, fontsize=15, labelpad=25)
    #
    # # from scipy.optimize import curve_fit
    # # K2_AS = K2_AS[K2_AS['Z'] < 0.]
    # # def func(x,a,c):
    # #     return a/x + c
    # # popt, pcov = curve_fit(func, K2_AS['Z'], K2_AS['delta_met'],sigma=K2_AS['feh_err'])
    # # print(popt)
    # # ax.plot(K2_AS['Z'],func(K2_AS['Z'],*popt),'g--')
    # import scipy.interpolate as interp
    # x = K2_AS['J']-K2_AS['Ks']
    # y = K2_AS['mbol']
    # xi, yi = np.linspace(x.min(),x.max(),500), np.linspace(y.min(),y.max(),500)
    # xi, yi = np.meshgrid(xi,yi)
    # rbf = interp.Rbf(x,y,K2_AS['delta_met'],function='linear')
    # zi = rbf(xi,yi)
    # # fig = plt.figure()
    # # cont = plt.contourf(xi,yi,zi,100,cmap=plt.cm.RdGy)
    # # plt.gca().invert_yaxis()
    # # cbar = plt.colorbar()
    # # cbar.set_label(r'$\Delta$ [Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # # plt.scatter(K2_AS['J']-K2_AS['Ks'],K2_AS['mbol'],alpha=0.07)
    # # plt.xlabel(r'J - K',fontsize=15)
    # # plt.ylabel(r'M$\_{\rm{bol}}$',fontsize=15)
    #
    # K2_New['feh_corr'] = 0
    # for j in range(len(K2_New)):
    #     JK = K2_New['J'].iloc[j] - K2_New['Ks'].iloc[j]
    #     Mbol = K2_New['mbol'].iloc[j]
    #     a = (np.abs(yi[:,0] - Mbol)).argmin()
    #     b = (np.abs(xi[0][:] - JK)).argmin()
    #     # print(b,(np.abs(xi[0][:] - JK)).argmin(),xi[0][b])
    #     # print(a,(np.abs(yi[:,0] - Mbol)).argmin(),yi[a][0])
    #     # print(zi[a][b])
    #     K2_New['feh_corr'].iloc[j] = zi[a][b]
    #
    # K2_New['corrected_feh'] = K2_New['feh'] - K2_New['feh_corr']
    #
    # # print(K2_New.columns.values)
    # plt.figure()
    # plt.scatter(K2_AS['feh'],K2_AS['delta_met'])
    # # cbar = plt.colorbar()
    # # plt.show()
    #
    # # plt.figure()
    # # plt.scatter(AS['Glon'],AS['Glat'],c=K2_AS['feh'],vmin=-3, vmax=1.1)
    # # cbar = plt.colorbar()
    #
    # plt.figure()
    # plt.scatter(C3_Sky['feh'],C3_Sky['delta_met'])#,c=C3_Sky['delta_met'],vmin=-1, vmax=1)
    # # cbar = plt.colorbar()
    # plt.show()
    # # print(np.median(abs(K2_AS['corrected_feh'] - AS['feh'])))
    # # K2_New[['#Id','corrected_feh']].to_csv('/home/bmr135/K2_Poles/Mass_Distr_In/K2_Photo_FeH_Corr',index=False)
    # # plt.show()
    # # sys.exit()

    ''' Specroscopic and Photometric Mass and Age Comparisons: C3 & C6 '''
    # data = [c_three,c3]
    # field = ['C3','C6']
    # p = ['mass','logAge']
    # ran = [np.linspace(0.5,2.5,20),np.linspace(8.5,10.5,100)]
    # xtag = [r'Mass [M$_{}\odot$]',r'$\rm{log}_{10}(\rm{Age})$']
    # for i in range(2):
    #     plt.figure()
    #     k = 0
    #     for j in p:
    #         plt.title(field[0])
    #         plt.subplot(2,1,k+1)
    #         mdf.histo(data[1],j,ran[k],xtag[k],0,r'K2 Photom')
    #         mdf.histo(data[0],j,ran[k],xtag[k],0,r'K2 Spectro')
    #         plt.legend(prop={'size':15},loc=k+1)
    #         k+=1
    #     if save_out == 1:
    #         plt.savefig(ext_fig+'Spectro_photom_distr_dir_comp_'+field[0]+'.png')
    # plt.show()

    ''' Mass/Age Distributions at Different Z Values '''
    # zm = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8]
    # za = [0.0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2]
    # mdf.Z_subplots(K2,TRI,'mass',np.linspace(0.5,2.5,20),r'Mass [M$_{\odot}$]',zm)
    # plt.show()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'Mass_as_funct_Z_K2.png')
    # mdf.Z_subplots(K2,TRI,'logAge',np.linspace(8.5,10.5,70),r'log$_{10}$(Age)',za)
    # plt.show()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'Age_as_funct_Z_K2.png')

    # print(TRI3['#Gc'])
    # sys.exit()

    ''' Radius Distributions '''
    APK2_alpha = APK2[APK2['alpha'] > 0.1]
    K2_alpha = AS[AS['alpha'] > 0.1]
    # TRI3a = TRI[TRI['#Gc'] == 1]
    # TRI3b = TRI[TRI['#Gc'] == 2]
    # TRI3c = TRI[TRI['#Gc'] == 3]
    K2_para = Luca[(Luca['parallax_error']/Luca['parallax'] < .1)]
    # print(len(Luca),len(K2_para))
    f, ax = plt.subplots()
    ax.hist(APK2['rad'],bins=np.linspace(0,20,50),label=r'APOKASC',color='grey',alpha=0.2,normed=True)
    # ax.hist(Kep_Sim['radius'],bins=np.linspace(0,20,100),label=r'APOKASC',color='grey',alpha=0.2,normed=True)
    # mdf.histo(APK2_alpha,'rad',np.linspace(0,20,100),r'Radius [R$_{\odot}$]',0,r'$\alpha$-rich Kepler')
    ax.hist(APK2_alpha['rad'],bins=np.linspace(0,20,100),histtype='step',label=r'APOKASC $\alpha$-rich',color='r',normed=True,linewidth=2,linestyle='--')
    # mdf.histo(AS,'rad',np.linspace(4,20,80),r'Radius [R$_{\odot}$]',0,r'K2 Spec.')

    # ax.hist(K2_para['rad'],bins=np.linspace(0,20,50),histtype='step',normed=True,label=r'K2 SM',linewidth=2)
    # plt.hist(K2_para['Rgaia'],bins=np.linspace(0,20,50),histtype='step',normed=True,label=r'K2$_{\rm{Gaia}}$',linewidth=2,color='m')
    # plt.hist(K2_para['Rgaia2'],bins=np.linspace(0,20,50),histtype='step',label=r'K2$_{Gaia}$: 1/$\varpi$',linewidth=2,color='k',linestyle='--')

    # ax.hist(AS['Rgaia'],bins=np.linspace(0,20,50),histtype='step',label=r'K2 Spec.',normed=True,linewidth=2)
    # ax.hist(AS['rad'],bins=np.linspace(0,20,50),histtype='step',label=r'K2 Spec.',normed=True,linewidth=2)
    ax.hist(Luca['rad'],bins=np.linspace(0,20,100),histtype='step',label=r'K2 SM',normed=True,linewidth=2)
    ax.hist(K2_alpha['rad'],bins=np.linspace(0,20,100),histtype='step',label=r'K2 $\alpha$-rich',normed=True,linewidth=2)
    # ax.hist(K2_alpha['Rgaia'],bins=np.linspace(0,20,50),histtype='step',label=r'$\alpha$-rich K2, Gaia',normed=True,linewidth=2)
    # ax.hist(TRI3a['radius'],bins=np.linspace(0,20,100),histtype='step',label=r'TRI: Thin',normed=True,linewidth=2)
    # ax.hist(TRI3b['radius'],bins=np.linspace(0,20,100),histtype='step',label=r'TRI: Thick',normed=True,linewidth=2)
    # mdf.histo(Luca,'rad',np.linspace(0,20,100),r'Radius [R$_{\odot}$]',0,r'K2')
    # mdf.histo(gold,'rad',np.linspace(0,20,100),r'Radius [R$_{\odot}$]',0,r'K2, Gold')
    ax.hist(APK2_alpha['rad'],bins=np.linspace(0,20,100),histtype='step',label='__nolegend__',color='r',normed=True,linewidth=2,linestyle='--')
    ax.text(0.1, 0.90, '(A)', horizontalalignment='center',verticalalignment='center',fontsize=15, transform=ax.transAxes,)
    ax.set_yticks([])
    ax.set_xlim(3.5,18)
    ax.set_xlabel(r'Radius [R$_{\odot}$]',fontsize=15)
    ax.legend(prop={'size':10})
    plt.close()
    # f.savefig('Radius_PARAM.pdf', bbox_inches='tight')
    # # pdf.savefig(f)
    # plt.show()
    # sys.exit()

#     # f, ax = plt.subplots()
#     # AS = AS[AS['alpha'] > -90]
#     # APK2 = APK2[APK2['alpha'] > -90]
#     # APK2 = APK2[APK2['mass'] > -90]
#     # APK2 = APK2[APK2['mass'] < 1.5]
#     # a = APK2[APK2['evstate'] == 1]
#     # b = APK2[APK2['evstate'] == 2]
#     # # APK2 = APK2[APK2['evstate'] == 2]
#     # # t = ax.scatter(APK2['rad'],APK2['Teff'],c=APK2['mass'],cmap=colormaps.parula,label=r'APOKASC')
#     # t = ax.scatter(b['rad'],b['Teff'],c=b['mass'],cmap=colormaps.parula,label=r'APOKASC')
#     # # t = ax.scatter(TRI3['radius'],TRI3['Teff'],c=TRI3['M_H'],cmap=colormaps.parula,label=r'TRILEGAL')
#     # # t = ax.scatter(AS['rad'],AS['feh'],c=AS['alpha'],cmap=colormaps.parula,label=r'K2 Spec.')
#     # ax.set_xlabel(r'Radius [R$_{\odot}$]',fontsize='15')
#     # ax.set_ylabel(r'$T_{\rm{eff}}$ [K]',fontsize='15')
#     # # ax.legend(prop={'size':15})
#     # ax.set_xlim(7,16)
#     # cbar = fig.colorbar(t,ax=ax)
#     # cbar.set_label(r'Mass', rotation=270, fontsize=15, labelpad=15)
#     # # f.savefig('TRI_Radius_Teff_FeH.pdf', bbox_inches='tight')
#     #
#     # # plt.show()
#     # # sys.exit()
#

    ''' Mass Distributions '''
    # f = plt.figure()
    # plt.hist(APK2['mass'],bins=np.linspace(0.7,2.0,50),label=r'APOKASC',color='grey',alpha=0.2,normed=True)
    # # plt.hist(a['mass'],bins=np.linspace(0.5,2.5,100),label=r'APOKASC - RGB',normed=True,histtype='step')
    # # plt.hist(b['mass'],bins=np.linspace(0.5,2.5,100),label=r'APOKASC - Clump',normed=True,histtype='step')
    # # plt.hist(TRI3['Mass'],bins=np.linspace(0.5,2.5,100),label=r'TRI',color='grey',alpha=0.2,normed=True)
    # # plt.hist(TRI3a['Mass'],bins=np.linspace(0.5,2.5,100),label=r'TRI - Thin',normed=True,histtype='step')
    # # plt.hist(TRI3b['Mass'],bins=np.linspace(0.5,2.5,100),label=r'TRI - Thick',normed=True,histtype='step')
    # plt.hist(APK2_alpha['mass'],bins=np.linspace(0.7,2.,50),histtype='step',label=r'$\alpha$-rich Kepler',color='r',normed=True,linewidth=2,linestyle='--')
    # # mdf.histo(AS,'mass',np.linspace(0.5,3.0,80),r'Radius [R$_{\odot}$]',0,r'K2 Spec.')
    # # mdf.histo(Luca,'rad',np.linspace(0,20,100),r'Radius [R$_{\odot}$]',0,r'K2')
    # plt.hist(Luca['mass'],bins=np.linspace(0.7,2.,50),histtype='step',label=r'K2 SM',normed=True,linewidth=2)
    # plt.hist(K2_alpha['mass'],bins=np.linspace(0.7,2.0,50),histtype='step',label=r'K2 $\alpha$-rich',normed=True,linewidth=2)
    # # mdf.histo(gold,'rad',np.linspace(0,20,100),r'Radius [R$_{\odot}$]',0,r'K2, Gold')
    # # plt.xlim(3.5,18)
    # plt.xlabel(r'Mass [M$_{\odot}$]',fontsize=15)
    # plt.yticks([])
    # plt.legend(prop={'size':10})
    # f.savefig('TRI_Mass_Thin_Thick.pdf', bbox_inches='tight')
    # plt.show()
    # sys.exit()
    # pdf.savefig(f)
#     #
#     #
#     # f, ax = plt.subplots()
#     # APK2 = APK2[APK2['feh'] > -0.25]
#     # # t = ax.scatter(APK2['mass'],APK2['Teff'],c=APK2['feh'],cmap=colormaps.parula,label=r'APOKASC')
#     # t = ax.scatter(TRI3['Mass'],TRI3['Teff'],c=TRI3['M_H'],cmap=colormaps.parula,label=r'TRILEGAL')
#     # # t = ax.scatter(AS['rad'],AS['feh'],c=AS['alpha'],cmap=colormaps.parula,label=r'K2 Spec.')
#     # ax.set_xlabel(r'Mass [M$_{\odot}$]',fontsize='15')
#     # ax.set_ylabel(r'$T_{\rm{eff}}$ [K]',fontsize='15')
#     # # ax.legend(prop={'size':15})
#     # cbar = fig.colorbar(t,ax=ax)
#     # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=15)
#     # # f.savefig('TRI_Mass_Teff_FeH.pdf', bbox_inches='tight')
#     # plt.show()
#     # sys.exit()
#

    ''' Mass vs logg scatter '''
    Luca = Luca[Luca['sig_age']/Luca['age']<0.35]
    # L1 = Luca[abs(Luca['Z']) < 0.5]
    # L2 = Luca[(abs(Luca['Z']) >= 0.5) & (abs(Luca['Z']) < 1.0)]
    # L3 = Luca[(abs(Luca['Z']) >= 1.0) & (abs(Luca['Z']) < 1.5)]
    # # L4 = Luca[(abs(Luca['Z']) >= 1.5) & (abs(Luca['Z']) < 2.0)]
    # # L5 = Luca[(abs(Luca['Z']) >= 2.0) & (abs(Luca['Z']) < 2.5)]
    L4 = Luca[(abs(Luca['Z']) >= 1.0) & (Luca['age']<7.0)]
    # print(L4.columns.values)
    # sys.exit()


    ''' AS Age cut '''
    # AS = AS[AS['sig_age']/AS['age']<0.35]
    # print(len(AS))
    # sys.exit()
    # AS = AS[(AS['age'] < 7) & (AS['alpha'] > 0.1)]
    APK2 = APK2[APK2['sig_age']/APK2['age']<0.35]
    L1 = AS[abs(AS['Z']) < 0.5]
    L2 = AS[(abs(AS['Z']) >= 0.5) & (abs(AS['Z']) < 1.0)]
    L3 = AS[(abs(AS['Z']) >= 1.0) & (abs(AS['Z']) < 1.5)]
    # L4 = AS[(abs(AS['Z']) >= 1.5) & (abs(AS['Z']) < 2.0)]
    L4 = AS[abs(AS['Z']) >= 1.5]
#     # L5 = AS[(abs(AS['Z']) >= 2.0) & (abs(AS['Z']) < 2.5)]
#     # L6 = AS[abs(AS['Z']) >= 2.5]
#
#     goldC = gold[(gold['rad'] < 9.5) | (gold['rad'] > 11.5)] # Crude clump
#     L1c = goldC[abs(goldC['Z']) < 0.5]
#     L2c = goldC[(abs(goldC['Z']) >= 0.5) & (abs(goldC['Z']) < 1.0)]
#     L3c = goldC[(abs(goldC['Z']) >= 1.0) & (abs(goldC['Z']) < 1.5)]
#     L4c = goldC[(abs(goldC['Z']) >= 1.5) & (abs(goldC['Z']) < 2.0)]
#     L5c = goldC[(abs(goldC['Z']) >= 2.0) & (abs(goldC['Z']) < 2.5)]
#     L6c = goldC[abs(goldC['Z']) >= 2.5]
#
    APK2 = APK2[APK2['age'] > 0.]
    APK2 = APK2[APK2['age'] < 19.9]
    # APK2 = APK2[APK2['sig_age']/APK2['age'] < .35]
    APK2 = APK2.reset_index(drop=True)
    AK1 = APK2[abs(APK2['Z']) < 0.5]
    AK2 = APK2[(abs(APK2['Z']) >= 0.5) & (abs(APK2['Z']) < 1.0)]
    AK3 = APK2[(abs(APK2['Z']) >= 1.0) & (abs(APK2['Z']) < 1.5)]
#     # print(len(L1),len(L2))
#     # print(len(AK1),len(AK2))

    L1 = L1.reset_index(drop=True)
    print(stats.ks_2samp(L3.age,AK3.age))
    # sys.exit()
    ''' Age as a function of Z '''
#     # fig, ((ax,ax1),(ax2,ax3),(ax4,ax5)) = plt.subplots(3,2,sharex='col',sharey=True) #
    fig, ((ax,ax1),(ax2,ax3)) = plt.subplots(2,2,sharex='col',sharey=True,figsize=(9,4))
    sns.despine(left=True)
    bins = np.linspace(0,20,40)
    sns.distplot(AK1['age'], color="grey", ax=ax,bins=bins,label=r'APOKASC; Z $< 0.5$ kpc')
    sns.distplot(L1['age'], color="orange", ax=ax, bins=bins,label=r'K2 Spec.; Z $< 0.5$ kpc')
    # ax.hist(AK1['age'],bins=bins,histtype='step',label='__nolegend__',normed=True,linewidth=2,color='grey',alpha=0.1)
    # ax.hist(L1['age'],bins=bins,histtype='step',label='__nolegend__',normed=True,linewidth=2,color='b',alpha=0.1)
    ax.set_yticks([])
    ax.set_xlim(0,20)
    ax.set_ylim(0,0.5)
    ax.set_xlabel('')
    ax.legend(loc=1)

    sns.distplot(AK2['age'], color="grey", ax=ax1,bins=bins,label='__nolegend__')
    sns.distplot(L2['age'], color="orange", ax=ax1,bins=bins,label=r'$0.5 <$ Z $< 1.0$ kpc')
    ax1.set_yticks([])
    ax1.set_xlim(0,20)
    ax1.set_xlabel('')
    ax1.legend(loc=1)

    sns.distplot(AK3['age'], color="grey", ax=ax2,bins=bins,label=r'__nolegend__')
    sns.distplot(L3['age'], color="orange", ax=ax2,bins=bins,label=r'$1.0 <$ Z $< 1.5$ kpc')
    ax2.set_yticks([])
    ax2.set_xlabel('')
    ax2.legend(loc=1)

    sns.distplot(L4['age'], color="orange", ax=ax3,bins=bins,label=r'Z $> 1.5$ kpc')
    # ax3.hist(L4['age'],bins=bins,histtype='step',label=r'__nolegend__',normed=True,linewidth=2,color='b',alpha=0.1)
    ax3.set_yticks([])
    ax3.set_xlabel('')
    ax3.legend(loc=1)
    # print(len(L3[(L3['age'] < 7) & (L3['alpha'] > 0.1) & (L3['age'] > 4)]),len(L3[(L3['age'] < 7)]),len(L3))
    # print(len(L4[(L4['age'] < 7) & (L4['alpha'] > 0.1) & (L4['age'] > 4)]),len(L4[(L4['age'] < 7)]),len(L4))
    # print(len(AS[(AS['age'] < 7) & (AS['alpha'] > 0.1) & (AS['age'] > 4)]),len(AS[(AS['age'] < 7)]),len(AS))
    # print(len(APK2[(APK2['age'] < 7) & (APK2['alpha'] > 0.1)]),len(APK2[(APK2['age'] < 7)]),len(APK2))
    # sns.distplot(L5['age'], color="b", ax=ax4,bins=bins,label=r'$2.0 <$ Z $< 2.5$ kpc')
    # ax4.set_yticks([])
    ax2.set_xlabel(r'Age [Gyr]')
    # ax4.legend()
    #
    # sns.distplot(L6['age'], color="b", ax=ax5,bins=bins,label=r'Z $> 2.5$ kpc')
    # ax5.set_yticks([])
    ax3.set_xlabel(r'Age [Gyr]')
    # ax5.legend()
    # plt.tight_layout()
    # fig.savefig('Age_Z_Spec_HQ.pdf', bbox_inches='tight')
    # pdf.savefig(fig)
    # pdf.close()
    # plt.show()
    plt.close()
    # sys.exit()
#
#
#     # plt.figure()
#     # plt.scatter(L1['mass'],L1['logg'],label=r'Z $< 0.5$',alpha=0.5)
#     # plt.scatter(L2['mass'],L2['logg'],label=r'$0.5 <=$ Z $< 1.0$',alpha=0.5)
#     # plt.scatter(L3['mass'],L3['logg'],label=r'$1.0 <=$ Z $< 1.5$',alpha=0.5)
#     # plt.scatter(L4['mass'],L4['logg'],label=r'$1.5 <=$ Z $< 2.0$',alpha=0.5)
#     # plt.scatter(L5['mass'],L5['logg'],label=r'$2.0 <=$ Z $< 2.5$',alpha=0.5)
#     # plt.scatter(L6['mass'],L6['logg'],label=r'Z $>= 2.5$',alpha=0.5)
#     # # Errorbar plot #
#     # # plt.errorbar(L1['mass'],L1['logg'],label=r'Z $< 0.5$',xerr=[abs(L1['mass_68L']-L1['mass']),abs(L1['mass']-L1['mass_68U'])],fmt='o',alpha=0.5)
#     # # plt.errorbar(L2['mass'],L2['logg'],label=r'$0.5 <=$ Z $< 1.0$',xerr=[abs(L2['mass_68L']-L2['mass']),abs(L2['mass']-L2['mass_68U'])],fmt='o',alpha=0.5)
#     # # plt.errorbar(L3['mass'],L3['logg'],label=r'$1.0 <=$ Z $< 1.5$',xerr=[abs(L3['mass_68L']-L3['mass']),abs(L3['mass']-L3['mass_68U'])],fmt='o',alpha=0.5)
#     # # plt.errorbar(L4['mass'],L4['logg'],label=r'$1.5 <=$ Z $< 2.0$',xerr=[abs(L4['mass_68L']-L4['mass']),abs(L4['mass']-L4['mass_68U'])],fmt='o',alpha=0.5)
#     # # plt.errorbar(L5['mass'],L5['logg'],label=r'$2.0 <=$ Z $< 2.5$',xerr=[abs(L5['mass_68L']-L5['mass']),abs(L5['mass']-L5['mass_68U'])],fmt='o',alpha=0.5)
#     # # plt.errorbar(L6['mass'],L6['logg'],label=r'Z $>= 2.5$',xerr=[abs(L6['mass_68L']-L6['mass']),abs(L6['mass']-L6['mass_68U'])],fmt='o',alpha=0.5)
#     # # plt.scatter(Luca['mass'],Luca['logg'],c=abs(Luca['Z']),label=r'Luca',alpha=0.75,cmap=colormaps.parula)
#     # # cbar = plt.colorbar()
#     # # cbar.set_label(r'Z$_{\rm{abs}}$ [Kpc]', rotation=270, fontsize=15, labelpad=15)
#     # plt.xlim(-0.07+min(Luca['mass']),max(Luca['mass'])+0.1)
#     # plt.ylim(-0.1+min(Luca['logg']),max(Luca['logg'])+0.1)
#     # plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
#     # plt.ylabel(r'logg', fontsize=15)
#     # plt.legend()
#     # plt.show()
#     # sys.exit()
#
    ''' KEPLER vs K2 simulations '''
#     # plt.figure()
#     # for i in range(100):
#     #     Kep_Sim2 = alt_sim_kep.sample(n=len(alt_sim))
#     #     Kep_Sim2 = Kep_Sim2.reset_index(drop=True)
#     #     plt.hist(Kep_Sim2['logAge'][np.isfinite(Kep_Sim2['logAge'])],bins=ranges[0],histtype='step',color='b',label=r'Kepler Sim',normed=True)
#     #     plt.hist(alt_sim['logAge'][np.isfinite(alt_sim['logAge'])],bins=ranges[0],histtype='step',color='orange',label=r'K2 Sim',normed=True)
#     # plt.legend([r'Kepler Sim',r'K2 Sim'],prop={'size':15},loc=2)
#     # plt.xlabel(r'log$_{10}$(Age)')
#     # # plt.title(r'$\forall$ R')
#     # plt.title(r'R $< 9$')
#     # plt.show()
#     # if save_out == 1:
#     #     plt.savefig(ext_fig+folder_loc+'Kep_multi_samp_K2_age_distr.png')
#     #
#     # thin_kep = Kep_Sim2[Kep_Sim2['#Gc'] == 1]
#     # thick_kep = Kep_Sim2[Kep_Sim2['#Gc'] == 2]
#     #
#     # # plt.figure()
#     # fig, axes = plt.subplots(2,1,sharex=True,sharey=True)
#     # ax0,ax1 = axes.flatten()
#     # ax0.hist([thick['logAge'],thin['logAge']],bins=ranges[0],stacked=True,label=[r'K2 Sim Thick',r'K2 Sim Thin'])
#     # # plt.xlabel(r'log$_{10}$(Age)')
#     # ax0.legend(prop={'size':15},loc=2)
#     # ax1.hist([thick_kep['logAge'],thin_kep['logAge']],bins=ranges[0],stacked=True,label=[r'Kepler Sim Thick',r'Kepler Sim Thin'])
#     # plt.xlabel(r'log$_{10}$(Age)')
#     # ax1.legend(prop={'size':15},loc=2)
#     # plt.show()
#     # if save_out == 1:
#     #     plt.savefig(ext_fig+folder_loc+'Kep_K2_age_distr.png')
#

    ''' Moving average for age-metallicity/alpha trends '''
#     # plt.figure()
#     # c_three = c_three.sort_values(['age'])
#     # c_three = c_three[c_three['feh'] > -5]
#     # c_three = c_three[c_three['alpha'] > -5]
#     # c_six = c_six.sort_values(['age'])
#     # c_six = c_six[c_six['feh'] > -5]
#     # c_six = c_six[c_six['alpha'] > -5]
#     #
#     # plt.scatter(c_three['age'],c_three['feh'])
#     # y_av = mdf.movingaverage(c_three['feh'],50)
#     # plt.plot(c_three['age'], y_av,'r')
#     # plt.xlabel(r'Age [Gyr]')
#     # plt.ylabel(r'[Fe/H]')
#     # plt.show()
#
#     ''' Kiel Diagram + Age KDE '''
#     APK2=APK2[APK2['mass']>0.]
#     a, b = 0.22, 0.79893 # Mosser(?)
#     # a, b = 0.263, 0.772 # Stello et al. 2009
#     mesa['logg'] = np.log10(27400 * (((mesa['dnu']/a)**(1/b))/3090) * np.sqrt(10**(mesa['logTe']/5777)))
# #
    mesa = pd.read_csv('MESA_track_example',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])
    mesa['logg']=4.434+np.log10(1.00 / ((10**mesa['logL'])/((10**mesa['logTe'])/5777)**4))
    mesa2 = pd.read_csv('MESA_track_example2',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])
    mesa2['logg']=4.434+np.log10(1.00 / ((10**mesa2['logL'])/((10**mesa2['logTe'])/5777)**4))
    mesa3 = pd.read_csv('MESA_track_example3',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])
    mesa3['logg']=4.434+np.log10(1.00 / ((10**mesa3['logL'])/((10**mesa3['logTe'])/5777)**4))
    mesa4 = pd.read_csv('MESA_track_example4',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])
    mesa4['logg']=4.434+np.log10(0.80 / ((10**mesa4['logL'])/((10**mesa4['logTe'])/5777)**4))
    mesa5 = pd.read_csv('MESA_track_example5',delimiter=r'\s+',skiprows=1,names=['age','logL','logTe','dnu','pi','mod','mod1'])
    mesa5['logg']=4.434+np.log10(0.80 / ((10**mesa5['logL'])/((10**mesa5['logTe'])/5777)**4))
#
#     young = C3_Luca[(C3_Luca['age'] < 2.0)]
#     APK2 = APK2[APK2['Z'] < 1.5]
#     # LucaZ1 = gold[abs(gold['Z']) < 1.]
#     # LucaZ2 = gold[abs(gold['Z']) > 1.]

    LucaZ1 = AS[abs(AS['Z']) < 1.]
    LucaZ2 = AS[abs(AS['Z']) > 1.]

    ''' Kiel Diagrams '''
    vmin = 0.7
    vmax = 2.25
    c = 'mass'
    # APK2 = APK2.sort_values(c)#,ascending=False)
    # AS = AS.sort_values(c,ascending=False)
    # Luca = Luca.sort_values(c,ascending=False)
    ### APOKASC, Z < 1. - significant sample lie below this value ###

    # f = plt.figure()
    # plt.scatter(APK2['Teff'],APK2['logg'],c=APK2[c],cmap=colormaps.parula,alpha=0.9,label='APOKASC',vmin=vmin, vmax=vmax,s=15)
    # plt.plot(10**mesa5['logTe'],mesa5['logg'],color='k',label=r'0.8 M$_{\odot}$; -0.5 dex',linestyle='--')
    # plt.plot(10**mesa['logTe'],mesa['logg'],color='k',label=r'1.0 M$_{\odot}$; -0.5 dex')
    # plt.plot(10**mesa4['logTe'],mesa4['logg'],color='m',label=r'0.8 M$_{\odot}$; -0.25 dex',linestyle='--')
    # plt.plot(10**mesa2['logTe'],mesa2['logg'],color='m',label=r'1.0 M$_{\odot}$; -0.25 dex')
    # cbar = plt.colorbar()
    # cbar.set_label(r'Mass [M$_{\odot}$]', rotation=270, fontsize=15, labelpad=25)
    # plt.arrow(4350, 2.8, 100, -0.08, head_width=0.15, head_length=0.2, fc='k', ec='grey')
    # # plt.gca().invert_xaxis()
    # # plt.gca().invert_yaxis()
    # plt.xlabel(r'$T_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log(g)', fontsize=15)
    # plt.legend()
    # plt.xlim(5600,4200)
    # plt.ylim(3.5,1.4)
    # plt.tight_layout()
    # f.savefig('APOKASC.pdf', bbox_inches='tight')
    # pdf.savefig(f)
    # plt.show()
    # sys.exit()

    ## Spectroscopic K2 sample ###
    # f = plt.figure()
    # plt.scatter(APK2['Teff'],APK2['logg'],alpha=0.05,label='__nolegend__',color='grey')
    # plt.scatter(AS['Teff'],AS['logg'],c=AS[c],cmap=colormaps.parula,alpha=0.9,label='K2 Spec.',vmin=vmin, vmax=vmax,s=15)
    # plt.plot(10**mesa5['logTe'],mesa5['logg'],color='k',label=r'0.8 M$_{\odot}$; -0.5 dex',linestyle='--')
    # plt.plot(10**mesa['logTe'],mesa['logg'],color='k',label=r'1.0 M$_{\odot}$; -0.5 dex')
    # plt.plot(10**mesa4['logTe'],mesa4['logg'],color='m',label=r'0.8 M$_{\odot}$; -0.25 dex',linestyle='--')
    # plt.plot(10**mesa2['logTe'],mesa2['logg'],color='m',label=r'1.0 M$_{\odot}$; -0.25 dex')
    # cbar = plt.colorbar()
    # cbar.set_label(r'Mass [M$_{\odot}$]', rotation=270, fontsize=15, labelpad=25)
    # plt.arrow(4350, 2.8, 100, -0.08, head_width=0.15, head_length=0.2, fc='k', ec='grey')
    # # plt.gca().invert_xaxis()
    # # plt.gca().invert_yaxis()
    # plt.xlabel(r'$T_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log(g)', fontsize=15)
    # plt.legend()
    # plt.xlim(5600,4200)
    # plt.ylim(3.5,1.4)
    # plt.tight_layout()
    # f.savefig('All_Spec.pdf', bbox_inches='tight')
    # plt.show()
    # sys.exit()
    # pdf.savefig(f)
#
    ### Luca photom, Z < 1. ###
    # f = plt.figure()
    # plt.scatter(LucaZ1['teff'],LucaZ1['logg'],c=LucaZ1['feh'],cmap=colormaps.parula,alpha=0.84,label=r'K2 $< 1$kpc')
    # plt.plot(10**mesa5['logTe'],mesa5['logg'],color='k',label=r'0.8 M$_{\odot}$; -0.5 dex',linestyle='--')
    # plt.plot(10**mesa['logTe'],mesa['logg'],color='k',label=r'1.0 M$_{\odot}$; -0.5 dex')
    # plt.plot(10**mesa4['logTe'],mesa4['logg'],color='orange',label=r'0.8 M$_{\odot}$; -0.25 dex',linestyle='--')
    # plt.plot(10**mesa2['logTe'],mesa2['logg'],color='orange',label=r'1.0 M$_{\odot}$; -0.25 dex')
    # cbar = plt.colorbar()
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # plt.xlabel(r'T$_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log$_{10}$(g)', fontsize=15)
    # plt.legend()
    # plt.xlim(5500,4200)
    # plt.ylim(3.5,1.8)
    # plt.tight_layout()
    # pdf.savefig(f)
    # f.savefig('K2_Zlt1.pdf', bbox_inches='tight')

    ### Luca photom, Z > 1. ###
    # f = plt.figure()
    # plt.scatter(LucaZ2['teff'],LucaZ2['logg'],c=LucaZ2['feh'],cmap=colormaps.parula,alpha=0.84,label=r'K2 $> 1$kpc')
    # plt.plot(10**mesa5['logTe'],mesa5['logg'],color='k',label=r'0.8 M$_{\odot}$; -0.5 dex',linestyle='--')
    # plt.plot(10**mesa['logTe'],mesa['logg'],color='k',label=r'1.0 M$_{\odot}$; -0.5 dex')
    # plt.plot(10**mesa4['logTe'],mesa4['logg'],color='orange',label=r'0.8 M$_{\odot}$; -0.25 dex',linestyle='--')
    # plt.plot(10**mesa2['logTe'],mesa2['logg'],color='orange',label=r'1.0 M$_{\odot}$; -0.25 dex')
    # cbar = plt.colorbar()
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # plt.xlabel(r'T$_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log$_{10}$(g)', fontsize=15)
    # plt.legend()
    # plt.xlim(5500,4200)
    # plt.ylim(3.5,1.8)
    # plt.tight_layout()
    # pdf.savefig(f)
    # f.savefig('K2_Zlt1.pdf', bbox_inches='tight')

    ### Luca photom, gold ###
    # f = plt.figure()
    # plt.scatter(APK2['Teff'],APK2['logg'],alpha=0.05,label='APOKASC',color='grey')
    # plt.scatter(gold['teff'],gold['logg'],c=gold[c],cmap=colormaps.parula,alpha=0.84,label=r'K2, Gold',vmin=vmin, vmax=vmax)
    # plt.plot(10**mesa5['logTe'],mesa5['logg'],color='k',label=r'0.8 M$_{\odot}$; -0.5 dex',linestyle='--')
    # plt.plot(10**mesa['logTe'],mesa['logg'],color='k',label=r'1.0 M$_{\odot}$; -0.5 dex')
    # plt.plot(10**mesa4['logTe'],mesa4['logg'],color='m',label=r'0.8 M$_{\odot}$; -0.25 dex',linestyle='--')
    # plt.plot(10**mesa2['logTe'],mesa2['logg'],color='m',label=r'1.0 M$_{\odot}$; -0.25 dex')
    # cbar = plt.colorbar()
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # plt.xlabel(r'T$_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log$_{10}$(g)', fontsize=15)
    # plt.legend()
    # plt.xlim(5500,4200)
    # plt.ylim(3.5,1.8)
    # plt.tight_layout()
    # f.savefig('Gold_feh.pdf', bbox_inches='tight')
    # pdf.savefig(f)
    # plt.show()
    # sys.exit()

    ### Luca photom ###
    # f = plt.figure()
    # plt.scatter(APK2['Teff'],APK2['logg'],alpha=0.05,label='__nolegend__',color='grey')
    # plt.scatter(Luca['teff'],Luca['logg'],c=Luca[c],cmap=colormaps.parula,alpha=0.9,label=r'K2 SM',vmin=vmin, vmax=vmax,s=15)
    # plt.plot(10**mesa5['logTe'],mesa5['logg'],color='k',label=r'0.8 M$_{\odot}$; -0.5 dex',linestyle='--')
    # plt.plot(10**mesa['logTe'],mesa['logg'],color='k',label=r'1.0 M$_{\odot}$; -0.5 dex')
    # plt.plot(10**mesa4['logTe'],mesa4['logg'],color='m',label=r'0.8 M$_{\odot}$; -0.25 dex',linestyle='--')
    # plt.plot(10**mesa2['logTe'],mesa2['logg'],color='m',label=r'1.0 M$_{\odot}$; -0.25 dex')
    # cbar = plt.colorbar()
    # cbar.set_label(r'Mass [M$_{\odot}$]', rotation=270, fontsize=15, labelpad=25)
    # plt.arrow(4350, 2.8, 100, -0.08, head_width=0.15, head_length=0.2, fc='k', ec='grey')
    # # plt.gca().invert_xaxis()
    # # plt.gca().invert_yaxis()
    # plt.xlabel(r'$T_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log(g)', fontsize=15)
    # plt.legend()
    # plt.xlim(5600,4200)
    # plt.ylim(3.5,1.4)
    # plt.tight_layout()
    # # f.savefig('Luca.pdf', bbox_inches='tight')
    # # pdf.savefig(f)
    # plt.show()
    # sys.exit()

    ''' Age Distribution KDEs '''
    APK2 = APK2[APK2['sig_age']/APK2['age'] < 0.35]
    APK2 = APK2[APK2['age'] > 0.]
    APK2_v2 = APK2[APK2['alpha'] < 0.1]
    APK2_alpha = APK2_alpha[APK2_alpha['sig_age']/APK2_alpha['age'] < 0.35]
    AS = AS[AS['sig_age']/AS['age'] < 0.35]
    Luca = Luca[Luca['sig_age']/Luca['age'] < 0.35]
    LucaZ1 = LucaZ1[LucaZ1['sig_age']/LucaZ1['age'] < 0.35]
    LucaZ2 = LucaZ2[LucaZ2['sig_age']/LucaZ2['age'] < 0.35]

    Luca3en['sig_age'] = abs(Luca3en['age_68U'] - Luca3en['age_68L'])/2
    Luca6en['sig_age'] = abs(Luca6en['age_68U'] - Luca6en['age_68L'])/2
    Luca3en = Luca3en[Luca3en['sig_age']/Luca3en['age'] < 0.35]
    Luca6en = Luca6en[Luca6en['sig_age']/Luca6en['age'] < 0.35]
    en = pd.concat([Luca3en,Luca6en],ignore_index=True)

    f1, (ax,ax1,ax2) = plt.subplots(1,3,figsize=(15,5))
    # f1, (ax,ax2) = plt.subplots(1,2,figsize=(10,5))
    d1 = kde.KDE1D(APK2['age'])
    x1 = np.r_[min(APK2['age']):max(APK2['age']):1024j]
    ax.plot(x1,d1(x1),linewidth=2,label=r'Full',color='grey')
    d1 = kde.KDE1D(APK2_alpha['age'])
    x1 = np.r_[min(APK2_alpha['age']):max(APK2_alpha['age']):1024j]
    ax.plot(x1,d1(x1),linewidth=2,label=r'[$\alpha$/Fe] $> 0.1$',color='r')
    d1 = kde.KDE1D(APK2_v2['age'])
    x1 = np.r_[min(APK2_v2['age']):max(APK2_v2['age']):1024j]
    ax.plot(x1,d1(x1),linewidth=2,label=r'[$\alpha$/Fe] $< 0.1$',linestyle='--')
    ax.legend(loc=8,prop={'size':8})

    d1 = kde.KDE1D(Luca['age'])
    x1 = np.r_[min(Luca['age']):max(Luca['age']):1024j]
    ax1.plot(x1,d1(x1),linewidth=2,label=r'Full',color='b')
    # d1 = kde.KDE1D(LucaZ1['age'])
    # x1 = np.r_[min(LucaZ1['age']):max(LucaZ1['age']):1024j]
    d1 = kde.KDE1D(en['age'])
    x1 = np.r_[min(en['age']):max(Luca3en['age']):1024j]
    ax1.plot(x1,d1(x1),linewidth=2,label=r'Z $< 1.0$',color='g')
    # d1 = kde.KDE1D(LucaZ2['age'])
    # x1 = np.r_[min(LucaZ2['age']):max(LucaZ2['age']):1024j]
    # ax1.plot(x1,d1(x1),linewidth=2,label=r'Z $> 1.0$',color='m')
    ax1.legend(loc=8,ncol=1,prop={'size':8})

    d1 = kde.KDE1D(AS['age'])
    x1 = np.r_[min(AS['age']):max(AS['age']):1024j]
    ax2.plot(x1,d1(x1),linewidth=2,label=r'Full',color='orange')
    for i in range(len(AS)):
        for j in range(len(AP3)):
            if AS['#Id'].iloc[i] == AP3['#Id'].iloc[j]:
                AS['APO'].iloc[i] = 1
    for i in range(len(AS)):
        for j in range(len(AP6)):
            if AS['#Id'].iloc[i] == AP6['#Id'].iloc[j]:
                AS['APO'].iloc[i] = 1
    # AS = AS[AS['APO'] == 1]
    ar = AS[AS['alpha']>0.1]
    ap = AS[AS['alpha']<=0.1]
    d1 = kde.KDE1D(ar['age'])
    x1 = np.r_[min(ar['age']):max(ar['age']):1024j]
    ax2.plot(x1,d1(x1),linewidth=2,label=r'[$\alpha$/Fe] $> 0.1$',color='r')
    d1 = kde.KDE1D(ap['age'])
    x1 = np.r_[min(ap['age']):max(ap['age']):1024j]
    ax2.plot(x1,d1(x1),linewidth=2,label=r'[$\alpha$/Fe] $< 0.1$',linestyle='--')
    ax2.legend(loc=8,ncol=1,prop={'size':8})

    # d1 = kde.KDE1D(TRI3['age'])
    # x1 = np.r_[min(TRI3['age']):max(TRI3['age']):1024j]
    # ax2.plot(x1,d1(x1),linewidth=2,label=r'TRILEGAL')#,color='m')
    # d1 = kde.KDE1D(TRI3a['age'])
    # x1 = np.r_[min(TRI3a['age']):max(TRI3a['age']):1024j]
    # ax2.plot(x1,d1(x1),linewidth=2,label=r'TRILEGAL (Thin)')#,color='m')
    # d1 = kde.KDE1D(TRI3b['age'])
    # x1 = np.r_[min(TRI3b['age']):max(TRI3b['age']):1024j]
    # ax2.plot(x1,d1(x1),linewidth=2,label=r'TRILEGAL (Thick)')#,color='m')
    # d1 = kde.KDE1D(TRI3c['age'])
    # x1 = np.r_[min(TRI3c['age']):max(TRI3c['age']):1024j]
    # ax2.plot(x1,d1(x1),linewidth=2,label=r'TRILEGAL ()')#,color='m')
    # ax2.legend()

    ax.set_yticks([])
    # ax1.set_yticks([])
    ax2.set_yticks([])
    ax.set_title(r'APOKASC')
    # ax1.set_title(r'K2 HQ')
    ax2.set_title(r'K2 Spec.$_{\rm{HQ}}$')
    ax.set_xlim(0,20)
    # ax1.set_xlim(0,20)
    ax2.set_xlim(0,20)
    ax.set_xlabel(r'Age [Gyr]',fontsize=15)
    # ax1.set_xlabel(r'Age [Gyr]',fontsize=15)
    ax2.set_xlabel(r'Age [Gyr]',fontsize=15)
    plt.subplots_adjust(hspace=0.07, wspace=0.1,top=0.91,bottom=0.15)
    # f1.savefig('Age_KDE_comb.pdf', bbox_inches='tight')
    plt.show()
    # pdf.savefig(f1)
    sys.exit()

    ''' Replication of Andrea's plot (10/09/2018) using spectroscopic K2 data (mass/age vs Z with [Fe/H] colour bar - date split by alpha) '''
    # fig, ((ax,ax1),(ax2,ax3)) = plt.subplots(2,2)
    # AS = AS[AS['age'] < 19.9]
    # AS = AS.reset_index(drop=True)
    # AS1 = AS[AS['alpha'] <= 0.1]
    # AS2 = AS[AS['alpha'] > 0.1]
    #
    # ax.grid(color='grey', linestyle='-', linewidth=1, alpha=0.25)
    # a = ax.scatter(AS1['mass'],abs(AS1['Z']),c=AS1['feh'],cmap=colormaps.parula,vmin=-2.5, vmax=0.5)
    # ax.title.set_text(r'[$\alpha$/Fe] $<$ 0.1')
    # ax.set_xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    # ax.set_ylabel(r'Z [kpc]', fontsize=15)
    # cbar = fig.colorbar(a, ax=ax)
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    #
    # ax1.grid(color='grey', linestyle='-', linewidth=1, alpha=0.25)
    # a = ax1.scatter(AS2['mass'],abs(AS2['Z']),c=AS2['feh'],cmap=colormaps.parula,vmin=-2.5, vmax=0.5)
    # ax1.title.set_text(r'[$\alpha$/Fe] $>$ 0.1')
    # ax1.set_xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    # ax1.set_ylabel(r'Z [kpc]', fontsize=15)
    # cbar = fig.colorbar(a, ax=ax1)
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    #
    # ax2.grid(color='grey', linestyle='-', linewidth=1, alpha=0.25)
    # ax2.grid(which='minor', color='grey', linestyle='--',alpha=0.25)
    # a = ax2.scatter(AS1['age'],abs(AS1['Z']),c=AS1['feh'],cmap=colormaps.parula,vmin=-2.5, vmax=0.5)
    # ax2.set_xlabel(r'Age [Gyr]', fontsize=15)
    # ax2.set_ylabel(r'Z [kpc]', fontsize=15)
    # ax2.set_xscale('log')
    # cbar = fig.colorbar(a, ax=ax2)
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    #
    # ax3.grid(color='grey', linestyle='-', linewidth=1, alpha=0.25)
    # ax3.grid(which='minor', color='grey', linestyle='--',alpha=0.25)
    # a = ax3.scatter(AS2['age'],abs(AS2['Z']),c=AS2['feh'],cmap=colormaps.parula,vmin=-2.5, vmax=0.5)
    # ax3.set_xlabel(r'Age [Gyr]', fontsize=15)
    # ax3.set_ylabel(r'Z [kpc]', fontsize=15)
    # ax3.set_xscale('log')
    # cbar = fig.colorbar(a, ax=ax3)
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    #
    # plt.tight_layout()
    # fig.savefig('Spec_Alpha.pdf', bbox_inches='tight')
    # plt.show()
    # sys.exit()

    ''' Multiple Plots for different age and Z combinations '''
    # sys.exit()
    # plt.figure()
    # young = C3_Luca[(C3_Luca['age'] < 2.0) & (abs(C3_Luca['Z']) > 1.0)]
    # # plt.scatter(APK2['Teff'],APK2['Av'],alpha=0.4,label=r'APOKASC')
    # plt.scatter(C3_Luca['teff'],C3_Luca['logg'],alpha=0.4,label=r'K2 (Full)')
    # plt.scatter(young['teff'],young['logg'],color='k')
    # # cbar = plt.colorbar()
    # # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # plt.xlabel(r'T$_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log$_{10}$(g)', fontsize=15)
    # plt.legend()
    # plt.tight_layout()
    #
    # plt.figure()
    # young = C3_Luca[(C3_Luca['age'] < 2.0) & (abs(C3_Luca['Z']) > 1.5)]
    # # plt.scatter(APK2['Teff'],APK2['Av'],alpha=0.4,label=r'APOKASC')
    # plt.scatter(C3_Luca['teff'],C3_Luca['logg'],alpha=0.4,label=r'K2 (Full)')
    # plt.scatter(young['teff'],young['logg'],color='k')
    # # cbar = plt.colorbar()
    # # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # plt.xlabel(r'T$_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log$_{10}$(g)', fontsize=15)
    # plt.legend()
    # plt.tight_layout()
    #
    # plt.figure()
    # young = C3_Luca[(C3_Luca['age'] < 2.0) & (abs(C3_Luca['Z']) > 2.0)]
    # # plt.scatter(APK2['Teff'],APK2['Av'],alpha=0.4,label=r'APOKASC')
    # plt.scatter(C3_Luca['teff'],C3_Luca['logg'],alpha=0.4,label=r'K2 (Full)')
    # plt.scatter(young['teff'],young['logg'],color='k')
    # # cbar = plt.colorbar()
    # # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # plt.xlabel(r'T$_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log$_{10}$(g)', fontsize=15)
    # plt.legend()
    # plt.tight_layout()
    #
    # plt.show()
    #
    # sys.exit()
    # plt.savefig(ext_fig+'Teff_logg.png')

    ''' Age vs Teff '''
    # plt.figure()
    # plt.scatter(AS['Teff'],AS['age'],c=AS['feh'],cmap=colormaps.magma)
    # cbar=plt.colorbar()
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.xlabel(r'T$_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'Age [Gyr]', fontsize=15)
    # plt.tight_layout()
    # plt.show()
    # plt.savefig(ext_fig+'Teff_Age_spectro.png')

    ''' Radius vs Scaling Radius '''
    # C3['Rs'] = (C3['nmx']/3090) * (C3['dnu']/135.1)**-2 * (C3['Teff']/5777)**0.5
    # C6['Rs'] = (C6['nmx']/3090) * (C6['dnu']/135.1)**-2 * (C6['Teff']/5777)**0.5
    # K2_New['Rs'] = (K2['nmx']/3090) * (K2['dnu']/135.1)**-2 * (K2['Teff']/5777)**0.5
    # # AS['Rs'] = (AS['nmx']/3090) * (AS['dnu']/135.1)**-2 * (AS['Teff']/5777)**0.5
    # # K2_AS['Rs'] = (K2_AS['nmx']/3090) * (K2_AS['dnu']/135.1)**-2 * (K2_AS['Teff']/5777)**0.5
    # K2_New['drad'] = (K2_New['rad']-K2_New['Rs'])/K2_New['rad']
    # # K2_New = K2_New[K2_New['drad'] > -0.5]
    # # K2_New = K2_New[K2_New['drad'] < 0.5]
    # plt.figure()
    # plt.scatter(K2_New['rad'],K2_New['drad'])
    # plt.plot([0,max(K2_New['rad'])+0.1],[0,0])
    # plt.xlim(0,max(K2_New['rad'])+0.1)
    # plt.plot([0,max(K2_New['rad'])+0.1],[0.5,0.5],color='r',linewidth=3)
    # plt.plot([0,max(K2_New['rad'])+0.1],[-0.5,-0.5],color='r',linewidth=3)
    # plt.xlabel(r'Radius [R$_{\odot}$]', fontsize=15)
    # plt.ylabel(r'R - R$_{sr}$', fontsize=15)
    # plt.tight_layout()
    # # plt.show()
    #
    # K2_New['Ms'] = (K2_New['nmx']/3090)**3 * (K2_New['dnu']/135.1)**-4 * (K2_New['Teff']/5777)**1.5
    # plt.figure()
    # plt.scatter(K2_New['mass'],(K2_New['mass']-K2_New['Ms'])/K2_New['mass'])
    # plt.plot([0,max(K2_New['mass'])+0.1],[0,0])
    # plt.xlim(0,max(K2_New['mass'])+0.1)
    # plt.plot([0,max(K2_New['mass'])+0.1],[0.5,0.5],color='r',linewidth=3)
    # plt.plot([0,max(K2_New['mass'])+0.1],[-0.5,-0.5],color='r',linewidth=3)
    # plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    # plt.ylabel(r'M - M$_{sr}$', fontsize=15)
    # plt.tight_layout()
    # # plt.show()

    # plt.figure()
    # hist, bins, patches = plt.hist(C3['rad'],bins=50,histtype='step',label=r'PARAM R',normed=True,linewidth=2)
    # plt.hist(C3['Rs'],bins=bins,histtype='step',label=r'Scaling R',normed=True,linewidth=2)
    # # plt.hist(alt_sim['radius'],bins=bins,histtype='step',label=r'TRI R',normed=True,linewidth=2,alpha=0.5)
    # plt.xlabel(r'Radius [R$_{\odot}$]', fontsize=15)
    # plt.legend()
    # plt.tight_layout()
    # plt.show()
    #
    # plt.figure()
    # hist, bins, patches = plt.hist(C6['rad'],bins=50,histtype='step',label=r'PARAM R',normed=True,linewidth=2)
    # plt.hist(C6['Rs'],bins=bins,histtype='step',label=r'Scaling R',normed=True,linewidth=2)
    # # plt.hist(alt_sim['radius'],bins=bins,histtype='step',label=r'TRI R',normed=True,linewidth=2,alpha=0.5)
    # plt.xlabel(r'Radius [R$_{\odot}$]', fontsize=15)
    # plt.tight_layout()
    # plt.legend()
    # plt.show()
    #
    # plt.figure()
    # hist, bins, patches = plt.hist(K2_AS['Rs'],bins=50,histtype='step',label=r'Photometry',normed=True,linewidth=2)
    # plt.hist(AS['Rs'],bins=bins,histtype='step',label=r'Spectroscopy',normed=True,linewidth=2)
    # # plt.xlabel(r'Scaling Relation Radius [R$_{\odot}$]', fontsize=15)
    # plt.legend()
    # plt.tight_layout()
    # plt.show()
    #
    # plt.figure()
    # hist, bins, patches = plt.hist(K2_AS['Teff'],bins=50,histtype='step',label=r'Photometry',normed=True,linewidth=2)
    # plt.hist(AS['Teff'],bins=bins,histtype='step',label=r'Spectroscopy',normed=True,linewidth=2)
    # plt.xlabel(r'Teff [K]', fontsize=15)
    # plt.legend()
    # plt.tight_layout()
    # plt.show()
    #
    # plt.figure()
    # plt.scatter(AS['logAge'],AS['rad'],c=AS['Teff'],cmap=colormaps.parula)
    # cbar = plt.colorbar()
    # cbar.set_label(r'T$_{\rm{eff}}$ [K]', rotation=270, fontsize=15, labelpad=25)
    # plt.xlabel(r'log Age', fontsize=15)
    # plt.ylabel(r'Radius [R$_{\odot}$]', fontsize=15)
    # plt.tight_layout()
    # plt.show()

    ''' Mass vs Radius ([Fe/H] colourbar) '''
    # plt.figure()
    # plt.scatter(K2['mass'],K2['rad'],c=K2['feh'],cmap=colormaps.magma)
    # cbar = plt.colorbar()
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    # plt.ylabel(r'Radius [R$_{\odot}$]', fontsize=15)
    # plt.tight_layout()
    # plt.show()

    ''' Mass vs Age (scaling mass colourbar) '''
    # K2['Ms'] = (K2['nmx']/3090)**3 * (K2['dnu']/135.1)**-4 * (K2['Teff']/5777)**1.5
    # plt.figure()
    # plt.scatter(K2['mass'],K2['logAge'],c=K2['Ms'],cmap=colormaps.parula)
    # cbar = plt.colorbar()
    # cbar.set_label(r'Mass - Scaling Relation', rotation=270, fontsize=15, labelpad=25)
    # plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    # plt.ylabel(r'log$_{10}$(Age)', fontsize=15)
    # plt.tight_layout()
    # plt.show()

    ''' Space Density plots '''
    # import mass_distr_functs as mdf
    # mdf.space_density2(C6_New)
    # mdf.space_density2(C3_New)
    # mdf.space_density2(K2_New)
    # mdf.space_density2(TRI3)
    # mdf.space_density2(TRI6)
    # # mdf.space_density2(AS)
    # mdf.sp3(C6)
    # plt.show()
    # sys.exit()

    ''' Mass vs Z scatter plots with trend lines '''
    # K2['mass_err'] = np.sqrt(((K2['mass_68U']-K2['mass'])**2 + (K2['mass']-K2['mass_68L'])**2)/4)
    # K2['Z_err'] = ((abs(K2['dist_68U'])-abs(K2['dist'])) + (abs(K2['dist'])-abs(K2['dist_68L'])))/2
    # K2['rad_err'] = 0.5*(K2['rad_68U']-K2['rad_68L'])
    # m_err = np.median(K2['mass_err'])
    # m_per = np.median(K2['mass_err']/K2['mass'])
    # z_err = np.median(K2['Z_err'])
    # z_per = np.median(K2['Z_err']*1e-3/K2['Z'])
    # K2['age'] = 10**K2['logAge'] * 1e-9
    # bins = np.linspace(0.75,2.25,10)
    # mass_Z, edges, number = scipy.stats.binned_statistic(C3['mass'],C3['feh'],statistic='median',bins=bins)
    # mass_Z6, edges6, number6 = scipy.stats.binned_statistic(C6['mass'],C6['feh'],statistic='median',bins=bins)
    #
    Luca3 = AS[AS['Z']<0]
    f = plt.figure()
    # plt.scatter(C6_New['mass'],C6_New['Z'],label=r'Full Sample')
    x = np.linspace(0,max(C6_New['mass'])+0.07)
    # plt.fill_between(x, 0.1, 1.5, facecolor='gray', alpha=0.2, interpolate=True,label=r'Kepler Z range')
    plt.scatter(APK2['mass'],APK2['Z'],c=APK2['feh'],cmap=colormaps.parula,label=None,vmin=-2.5,vmax=0.75)
    plt.scatter(Luca3['mass'],Luca3['Z'],c=Luca3['feh'],cmap=colormaps.parula,label=None,vmin=-2.5,vmax=0.75)
    # plt.errorbar(2.25, -3, xerr=m_err, yerr=z_err*1e-3,color='k')
    cbar = plt.colorbar()
    cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.plot(edges[:-1],mass_Z,color='k',linewidth=2)
    # plt.plot(edges6[:-1],abs(mass_Z6),color='k',linewidth=2)
    plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    plt.ylabel(r'Z [kpc]', fontsize=15)
    plt.xlim(min(C6_New['mass'])-0.07,max(C6_New['mass'])+0.07)
    # f.set_rasterized(True)
    # plt.title(r'Exploring Spectroscopic Peak Age Properties')
    # plt.legend()
    plt.tight_layout()
    # f.savefig('/home/bmr135/Dropbox/GES-K2/Ages/figure_1.eps',rasterized=True,dpi=400)


    fig, ax = plt.subplots()
    x = pd.DataFrame()
    x['V'] = [9.5,10.5,11.5,12.5,13.5,14.5]
    x['n'] = [18,37,66,115,97,22]
    ax.bar(x.V,x.n)
    ax.set_xlabel(r'V',fontsize=15)
    plt.close()
    # plt.show()
    # sys.exit()

    # mu_m = C3_New['sig_mass'].mean()
    # mu_Z = C3_New['sig_Z'].mean()
    # plt.figure()
    # x = np.linspace(0,max(C6_New['mass'])+0.07)
    # # plt.fill_between(x, 0.1, 1.5, facecolor='gray', alpha=0.2, interpolate=True,label=r'Kepler Z range')
    # plt.scatter(APK2['mass'],APK2['Z'],c=APK2['feh'],cmap=colormaps.parula,label=None,vmin=-2.5,vmax=0.5)
    # plt.scatter(C3_New['mass'],C3_New['Z'],c=C3_New['[Fe/H]'],cmap=colormaps.parula,label=None,vmin=-2.5,vmax=0.5)
    # cbar = plt.colorbar(spacing='uniform')
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # # cbar.set_clim(-2.5,1.25)
    # # cbar.set_ticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0],update_ticks=True)
    # # cbar.set_ticklabels([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0],update_ticks=True)
    # plt.plot([2.1,2.1],[-4.25-mu_Z,-4.25+mu_Z],color='k',linewidth=2,alpha=0.7,label=None)
    # plt.plot([2.1-mu_m,2.1+mu_m],[-4.25,-4.25],color='k',linewidth=2,alpha=0.7,label=None)
    # plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    # plt.ylabel(r'Z [kpc]', fontsize=15)
    # plt.xlim(min(C6_New['mass'])-0.07,max(C6_New['mass'])+0.07)
    # # plt.ylim(-6,8)
    # # plt.legend()
    # plt.tight_layout()
    # print(max(C3_New['[Fe/H]']),max(APK2['feh']))
    # plt.show()
    #
    # mu_m = AS['sig_mass'].mean()
    # mu_Z = AS['sig_Z'].mean()
    # plt.figure()
    # x = np.linspace(0,max(C6_New['mass'])+0.07)
    # plt.fill_between(x, 0.1, 1.5, facecolor='gray', alpha=0.2, interpolate=True,label=r'Kepler Z range')
    # plt.scatter(AS['mass'],AS['Z'],c=AS['feh'],cmap=colormaps.parula,label=None)
    # cbar = plt.colorbar(spacing='uniform')
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # cbar.set_clim(-2.5,1.25)
    # cbar.set_ticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0],update_ticks=True)
    # cbar.set_ticklabels([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0],update_ticks=True)
    # cbar.update_ticks()
    # plt.plot([2.1,2.1],[-4.25-mu_Z,-4.25+mu_Z],color='k',linewidth=2,alpha=0.7,label=None)
    # plt.plot([2.1-mu_m,2.1+mu_m],[-4.25,-4.25],color='k',linewidth=2,alpha=0.7,label=None)
    # plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    # plt.ylabel(r'Z [kpc]', fontsize=15)
    # plt.xlim(min(C6_New['mass'])-0.07,max(C6_New['mass'])+0.07)
    # plt.ylim(-6,8)
    # plt.legend()
    # plt.tight_layout()
    # plt.show()

    ''' Joint Photom and Spec MZFeH '''
    mu_m = Luca['sig_mass'].mean()
    mu_Z = Luca['sig_Z'].mean()

    # print(Luca['rad'])
    # sys.exit()
    # Lucaen = pd.merge(Lucaen,Luca[['#Id']],how='inner',on=['#Id'])
    # Lucaen = Lucaen.reset_index(drop=True)
    # ASen = pd.merge(ASen,AS[['#Id']],how='inner',on=['#Id'])
    # ASen = ASen.reset_index(drop=True)
    # # Luca = Luca[(Luca['rad'] <= 9.5) | (Luca['rad'] >= 11.5)]
    # # Luca=Luca.reset_index(drop=True)
    # # AS = AS[(AS['rad'] <= 9.5) | (AS['rad'] >= 11.5)]
    # # AS=AS.reset_index(drop=True)
    # # APK2 = APK2[(APK2['rad'] <= 9.5) | (APK2['rad'] >= 11.5)]
    # # APK2=APK2.reset_index(drop=True)
    # # Luca = Luca[Luca['age'] <= 9.]
    # # Luca=Luca.reset_index(drop=True)
    # # AS = AS[AS['age'] <= 9.]
    # # AS=AS.reset_index(drop=True)
    # # APK2 = APK2[APK2['age'] <= 9.]
    # # APK2=APK2.reset_index(drop=True)
    #
    # fig, (ax1,ax2,cax) = plt.subplots(ncols=3, gridspec_kw={"width_ratios" : [5,5,0.2]})
    # x = np.linspace(0,max(Luca['mass'])+0.07)
    # # ax1.fill_between(x, 0.1, 1.5, facecolor='gray', alpha=0.1, interpolate=True,label=r'Kepler $Z$ range')
    # ax1.scatter(APK2['mass'],APK2['Z'],alpha=0.05,cmap=colormaps.parula,vmin=-2.0, vmax=1.0,s=10,color='gray')
    # one = ax1.scatter(Luca['m_seismo'],Luca['Z'],c=Luca['feh'],cmap=colormaps.parula,label=None,vmin=-2.0, vmax=1.0,alpha=0.75,s=10)
    # # ax1.plot([2.1,2.1],[-4.25-mu_Z,-4.25+mu_Z],color='k',linewidth=2,alpha=0.7,label=None)
    # # ax1.plot([2.1-mu_m,2.1+mu_m],[-4.25,-4.25],color='k',linewidth=2,alpha=0.7,label=None)
    # ax1.set_xlabel(r'Mass [M$_{\odot}$], K2 SM', fontsize=15)
    # ax1.set_ylabel(r'Z [kpc]', fontsize=15)
    # ax1.set_xlim(min(Luca['mass'])-0.07,2.5)#max(Luca['mass'])+0.07)
    # ax1.set_ylim(-5,5)
    # # ax1.legend()
    #
    # mu_m = AS['sig_mass'].mean()
    # mu_Z = AS['sig_Z'].mean()
    # x = np.linspace(0,max(AS['mass'])+0.07)
    # # ax2.fill_between(x, 0.1, 1.5, facecolor='gray', alpha=0.2, interpolate=True,label=r'Kepler Z range')
    # # ax2.fill_between(x, 0.1, 1.5, facecolor='gray', alpha=0.1, interpolate=True,label='__nolegend__')
    # ax2.scatter(APK2['mass'],APK2['Z'],alpha=0.05,cmap=colormaps.parula,vmin=-2.0, vmax=1.0,s=10,color='gray')
    # two = ax2.scatter(AS['m_seismo'],AS['Z'],c=AS['feh'],cmap=colormaps.parula,label=None,vmin=-2.0, vmax=1.0,alpha=0.75,s=10)
    # ax2.text(0.9, 0.95, '(B)', horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes,fontsize=15)
    # # ax2.plot([2.1,2.1],[-4.25-mu_Z,-4.25+mu_Z],color='k',linewidth=2,alpha=0.7,label=None)
    # # ax2.plot([2.1-mu_m,2.1+mu_m],[-4.25,-4.25],color='k',linewidth=2,alpha=0.7,label=None)
    # ax2.set_xlabel(r'Mass [M$_{\odot}$], K2 Spec.', fontsize=15)
    # # ax2.set_ylabel(r'Z [kpc]', fontsize=15)
    # ax2.set_xlim(min(Luca['mass'])-0.07,2.5)#max(Luca['mass'])+0.07)
    # ax2.set_ylim(-5,5)
    # ax2.set_yticklabels([])
    # # ax2.legend()
    #
    # # cax.set_yticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0])
    # # cax.set_yticklabels(['-2.0','-1.5','-1.0','-0.5','0.0','0.5','1.0'])
    # cbar = fig.colorbar(two,cax=cax)
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=15)
    # cbar.set_clim(-2.0,1.0)
    # cbar.ax.set_yticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0])#,update_ticks=True)
    # cbar.ax.set_yticklabels(['-2.0','-1.5','-1.0','-0.5','0.0','0.5','1.0'])#,update_ticks=True)
    # # cbar.update_ticks()
    # # plt.tight_layout()
    # # fig.savefig('MZ_feh_massSR.pdf', bbox_inches='tight')
    # # pdf.savefig(fig)
    # plt.show()
    # sys.exit()

    ''' M vs Z - no colourbar '''
    # # mu_m = (C3_New['sig_mass'].mean() + C6_New['sig_mass'].mean())/2
    # # mu_Z = (C3_New['sig_Z'].mean() + C6_New['sig_Z'].mean())/2
    # # C6_New = C6_New[abs((C6_New['m_seismo']-C6_New['mass'])/C6_New['m_seismo']) < 0.1]
    # # C3_New = C3_New[abs((C3_New['m_seismo']-C3_New['mass'])/C3_New['m_seismo']) < 0.1]
    # f = plt.figure()
    # plt.scatter(APK2['mass'],APK2['Z'],c=APK2['feh'],alpha=0.7,cmap=colormaps.parula,vmin=-2.0, vmax=1.0)
    # # plt.scatter(C6_New['mass'],C6_New['Z'],label=r'C6',alpha=0.7)
    # # plt.scatter(C3_New['mass'],C3_New['Z'],label=r'C3',alpha=0.7)
    # cbar = plt.colorbar()
    # cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=15)
    # cbar.set_clim(-2.0,1.0)
    # cbar.ax.set_yticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0])#,update_ticks=True)
    # cbar.ax.set_yticklabels(['-2.0','-1.5','-1.0','-0.5','0.0','0.5','1.0'])#,update_ticks=True)
    # # plt.plot([2.1,2.1],[-4.25-mu_Z,-4.25+mu_Z],color='k',linewidth=2,alpha=0.7,label=None)
    # # plt.plot([2.1-mu_m,2.1+mu_m],[-4.25,-4.25],color='k',linewidth=2,alpha=0.7,label=None)
    # plt.xlabel(r'Mass [M$_{\odot}$], APOKASC',fontsize=15)
    # plt.ylabel(r'Z [kpc]',fontsize=15)
    # plt.xlim(min(C6_New['mass'])-0.07,max(C6_New['mass'])+0.07)
    # # plt.legend()
    # plt.tight_layout()
    # # f.savefig('MZ_feh_APOKASC.pdf', bbox_inches='tight')
    # # plt.show()
    # # sys.exit()

    ''' Age vs Z '''
    # plt.figure()
    # plt.scatter(K2['age'],K2['Z'])
    # plt.xlabel(r'Age', fontsize=15)
    # plt.ylabel(r'Z [kpc]', fontsize=15)
    # plt.tight_layout()
    # # plt.show()

    # plt.figure()
    # K2_RGB = K2[K2['rad'] < 10.0]
    # AS_RGB = AS[AS['rad'] < 10.0]
    # plt.hist(C6_New['logAge'],bins=np.linspace(8.5,10.5,40),histtype='step',normed=True,label=r'C6 Photometry',linewidth=2)
    # plt.hist(RC6['logAge'],bins=np.linspace(8.5,10.5,40),histtype='step',normed=True,label=r'RAVE C6',linewidth=2,alpha=0.65)
    # # plt.hist(np.log10(APK2['age']*10**9),bins=np.linspace(8.5,10.5,75),histtype='step')#,normed=True)
    # plt.xlabel(r'log$_{10}$(Age)')
    # cur_axes = plt.gca()
    # cur_axes.axes.get_yaxis().set_ticklabels([])
    # cur_axes.axes.get_yaxis().set_ticks([])
    # # plt.title(r'Clump cut: R $< 10.0$')
    # plt.legend()
    # plt.tight_layout()
    # plt.show()

    # APK2.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/APK2',index=False)
    # RC3.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/RC3',index=False)
    # GES.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/GES',index=False)
    # C6_New.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/C6_New',index=False)
    # RC6.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/RC6',index=False)

    # C3_New['sig_age'] = ((C3_New['age_68U']-C3_New['age']) + (C3_New['age']-C3_New['age_68L']))/2
    # C3_New.to_csv('/home/bmr135/K2_BG/C3_New',index=False)
    # C6_New['sig_age'] = ((C6_New['age_68U']-C6_New['age']) + (C6_New['age']-C6_New['age_68L']))/2
    # C6_New.to_csv('/home/bmr135/K2_BG/C6_New',index=False)

    ''' [Fe/H] vs [Alpha/Fe] '''
    AS = AS[AS['alpha'] > -4]
    RG = pd.merge(RC3,GES,how='inner',on=['#Id'])
    GES_rev = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Gaia_ESO/Alpha_Abundances.dat',delimiter=r'\s+')
    GES_rev = GES_rev.rename(columns={'EPIC':'#Id'})
    GES_rev = pd.merge(GES_rev,AS[['#Id']],how='inner',on=['#Id'])
    AP3_rev = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/APOGEE_DR14_C3_180418')
    AP3_rev = AP3_rev.rename(columns={'EPIC':'#Id'})
    AP3_rev = pd.merge(AP3_rev,AS[['#Id']],how='inner',on=['#Id'])
    AP6_rev = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/APOGEE_DR14_C6_180418')
    AP6_rev = AP6_rev.rename(columns={'EPIC':'#Id'})
    AP6_rev = pd.merge(AP6_rev,AS[['#Id']],how='inner',on=['#Id'])
    data = fits.getdata('/home/bmr135/Downloads/APOKASC4BEN_joined_allStar_l31c2.fits', 1)
    t = Table(data)
    cols = ['NINST','STABLERV_CHI2','STABLERV_RCHI2','CHI2_THRESHOLD','STABLERV_CHI2_PROB', \
            'PARAM','FPARAM','PARAM_COV','FPARAM_COV','PARAMFLAG','FELEM','FELEM_ERR', \
            'X_H','X_H_ERR','X_M','X_M_ERR','ELEM_CHI2','ELEMFLAG','ALL_VISIT_PK', \
            'VISIT_PK','FPARAM_CLASS','CHI2_CLASS']
    cols_corogee = ['PARAM_COV_2','FPARAM_COV_2','FELEM_2','FELEM_ERR_2','X_H','X_H_ERR', \
                    'X_M','X_M_ERR','ELEM_CHI2_2','ELEMFLAG_2']
    t.remove_columns(cols)
    flag = t['ASPCAPFLAG'] <= 5
    APK2_rev = t[flag].to_pandas()
    APK2_rev = APK2_rev[APK2_rev['FE_H'] > -2]
    APK2_rev = APK2_rev[APK2_rev['MG_FE'] > -5]
    # sys.exit()
    f2 = plt.figure()
    # print(len(APK2),len(APK2[APK2['alpha']>0.1]))
    APK2 = APK2[APK2['alpha'] > -5]
    RC3 = RC3[RC3['alpha'] > -5]
    RC6 = RC6[RC6['alpha'] > -5]
    AP3 = AP3[AP3['alpha'] > -5]
    AP6 = AP6[AP6['alpha'] > -5]
    GES = GES[GES['alpha'] > -5]
    # plt.scatter(AS['feh'],AS['ALPHA'],marker='<',label=r'K2 Spec.')#,c=AS['age'],cmap=colormaps.parula)
    plt.scatter(APK2_rev.FE_H,APK2_rev['MG_FE'],label=r'APOKASC',alpha=0.025,color='gray',s=10)
    # plt.scatter(AS['feh'],AS['alpha'],label=r'K2 Spec.')
    # plt.scatter(RC3['feh'],RC3['alpha'],label=r'RAVE C3',s=20,alpha=0.75)
    # plt.scatter(RC6['feh'],RC6['alpha'],label=r'RAVE C6',s=20,alpha=0.75)
    plt.scatter(AP6_rev['FE_H'],AP6_rev['MG_FE'],label=r'APOGEE C6',s=20,alpha=0.75)
    plt.scatter(AP3_rev['FE_H'],AP3_rev['MG_FE'],label=r'APOGEE C3',s=20,alpha=0.75)
    plt.scatter(GES_rev['FEH'],GES_rev['MG1'],label=r'Gaia-ESO C3',s=20,alpha=0.75,marker='d')
    # plt.scatter(GES['feh'],GES['alpha'],label=r'Gaia-ESO C3',s=20,alpha=0.75)
    # plt.scatter(AS['feh'],AS['alpha_y'],marker='>',label=r'Gaia-ESO')#,c=AS['age'],cmap=colormaps.parula)
    # cbar = plt.colorbar()
    # cbar.set_label(r'Age [Gyr]', rotation=270, fontsize=15, labelpad=25)
    plt.xlabel(r'[Fe/H]', fontsize=15)
    plt.ylabel(r'[Mg/Fe]', fontsize=15)
    # plt.title(r'RAVE Gaia-ESO cross match')
    plt.legend(loc=3)
    plt.tight_layout()
    plt.close()
    # f2.savefig('Mg_fe.pdf',bbox_inches='tight')
    # plt.show(f2)

    # fig, ((ax,ax1),(rax,rax1)) = plt.subplots(2,2,sharex='col',gridspec_kw={"height_ratios" : [5,1]})
    # ax.plot([-0.9,0.4],[-0.9,0.4],linestyle='--',color='k')
    # ax.scatter(RG['feh_x'],RG['feh_y'])
    # ax.set_ylabel(r'[Fe/H] - Gaia-ESO', fontsize=15)
    # ax.set_xlim(-0.9,0.4)
    #
    # ax1.plot([-0.1,0.25],[-0.1,0.25],linestyle='--',color='k')
    # ax1.scatter(RG['ALPHA'],RG['alpha_y'])
    # ax1.set_ylabel(r'[$\alpha$/Fe] - Gaia-ESO', fontsize=15)
    # ax1.set_xlim(-0.1,0.25)
    #
    # rax.plot([-0.9,0.4],[0,0],color='k')
    # rax.scatter(RG['feh_x'],(RG['feh_x']-RG['feh_y'])/RG['feh_x'])
    # rax.set_xlabel(r'[Fe/H] - RAVE')
    # rax.set_ylabel(r'Residual [%]')
    # rax1.plot([-0.1,0.25],[0,0],color='k')
    # rax1.scatter(RG['ALPHA'],(RG['ALPHA']-RG['alpha_y'])/RG['ALPHA'])
    # rax1.set_xlabel(r'[$\alpha$/Fe] - RAVE')
    # rax1.set_ylabel(r'Residual [%]')
    # plt.tight_layout()
    # plt.show()
    # sys.exit()

    ''' Mag vs Z with trend lines '''
    # bins = np.linspace(-8,0,16)
    # bins6 = np.linspace(0,8,16)
    # Kp_Z, edges, number = scipy.stats.binned_statistic(C3['Z'],C3['Kepler'],statistic='median',bins=bins)
    # Kp_Z6, edges6, number6 = scipy.stats.binned_statistic(C6['Z'],C6['Kepler'],statistic='median',bins=bins6)
    # plt.figure()
    # plt.scatter(K2['Z'],K2['Kepler'],c=K2['mass'],cmap=colormaps.parula)
    # cbar = plt.colorbar()
    # cbar.set_label(r'Mass [M$_{\odot}$]', rotation=270, fontsize=15, labelpad=25)
    # plt.plot(edges[:-1],Kp_Z,color='k')
    # plt.plot(edges6[1:],Kp_Z6,color='k')
    # plt.xlabel(r'Z [kpc]', fontsize=15)
    # plt.ylabel(r'Kp', fontsize=15)
    # plt.tight_layout()
    # plt.show()
    #
    # Kp_Z, edges, number = scipy.stats.binned_statistic(C3_40_06FeH['Z'],C3_40_06FeH['Kepler'],statistic='count',bins=bins)
    # Kp_Z6, edges6, number6 = scipy.stats.binned_statistic(C6_40_06FeH['Z'],C6_40_06FeH['Kepler'],statistic='count',bins=bins6)

    ''' Distibution of stars by Galactic Radius and Z '''
    # print(gold.columns.values)
    mu_GR = (C3_New['sig_Gal_Rad'].mean() + C6_New['sig_Gal_Rad'].mean())/2
    mu_Z = (C3_New['sig_Z'].mean() + C6_New['sig_Z'].mean())/2
    AS_clump = AS[(AS['rad'] > 10.0) & (AS['rad'] < 12.0)]
    Yu1 = pd.read_csv('/home/bmr135/K2_Poles/Yu_2018_T1')
    Yu2 = pd.read_csv('/home/bmr135/K2_Poles/Yu_2018_T2')
    Yu = pd.merge(Yu1,Yu2,how='inner',on=['#Id'])
    id = pd.read_csv('/home/bmr135/Downloads/kepler_search(2).txt')
    id = id.drop_duplicates(subset=['#Id'])
    Yu = pd.merge(Yu,id,how='inner',on=['#Id'])
    # print(max(Yu.Glat),max(C6_Luca.Glat))
    # sys.exit()


    # cols = ['#Id', 'logg', 'FeH', 'Teff', 'EBV']
    # Yu['EBV'] = 0
    # g2 = Yu[cols]
    # g2.to_csv('input.sample.all', header=None, sep=r' ', index=False)
    # p = subprocess.Popen(['./bcall'])
    # p.wait()
    # p = subprocess.Popen(["mv", "output.file.all", "/home/bmr135/K2_Poles/Yu_BCs.csv"])
    # p.wait()
    # BC = pd.read_table('/home/bmr135/K2_Poles/Yu_BCs.csv', sep=r'\s+', header=0)
    # BC = BC[['ID', 'BC_1']]
    # BC = BC.rename(columns={'ID':'#Id','BC_1':'BC_K'})
    # Yu = pd.merge(Yu, BC, on=['#Id'], how='inner')
    # # print(Yu.columns.values)
    # Yu['A_K'] = 0
    # MbSun = 4.75 # Casagrande et al. 2018
    # Yu['dist'] = 10 * 10**((Yu.K + 2.5*np.log10((Yu.R_RGB**2) * (Yu.Teff/5777)**4) - 4.75 - Yu.A_K + Yu.BC_K)/5)
    # Yu = mdf.vert_dist(Yu)


    # sys.exit()
    # AS = AS[AS['age'] <= 7.0]
    # AS = AS[AS['alpha'] < -10]
    APK2 = APK2[APK2['Z'] < 1.5]
    # print(len(AS[AS['alpha']>0.1]))
    K2_ar = AS[(AS['alpha'] > 0.1) & ((10**AS['logAge']/1e9) < 7)]
    # K2_ar = K2_ar[abs((K2_ar['mass']-K2_ar['m_seismo'])/K2_ar['mass']) < 0.5]
    K2_ar = K2_ar.reset_index(drop=True)
    f = plt.figure()
    # print(K2_New['feh'])
    plt.plot([min(Luca['Gal_Rad'])-0.2,max(APK2['Gal_Rad'])+0.2],[0,0],color='r',linewidth=2,alpha=0.7,linestyle='--',label=None)
    # plt.scatter(Yu['Gal_Rad'],Yu['Z'],color='g',alpha=0.05,label=r'Yu et al. 2018',marker='d')
    plt.scatter(APK2['Gal_Rad'],APK2['Z'],color='k',alpha=0.15,label=r'APOKASC')
    # plt.scatter(Luca['Gal_Rad'],Luca['Z'],label=r'K2 Spec.$_{\rm{HQ}}$'),alpha=1.,c=AS['age'],cmap=colormaps.parula,vmin=0,vmax=14,edgecolor='grey')
    # plt.scatter(C6_Luca['Gal_Rad'],C6_Luca['Z'],label=r'C6',alpha=0.4)
    # plt.scatter(C3_Luca['Gal_Rad'],C3_Luca['Z'],label=r'C3',alpha=0.4)
    # plt.scatter(K2_ar['Gal_Rad'],K2_ar['Z'],label=r'Young $\alpha$-rich',alpha=0.7,color='y',edgecolor='k')
    # plt.scatter(gold['Gal_Rad'],gold['Z'],label=r'Gold',edgecolors='g',facecolor='none')
    plt.scatter(Luca['Gal_Rad'],Luca['Z'],c=Luca['age'],cmap=colormaps.parula,label=r'K2 SkyMapper',edgecolor='k')
    # plt.scatter(gold['Gal_Rad'],gold['Z'],c=gold['age'],label=r'Gold',marker='^',edgecolors='r')
    plt.plot([14.5,14.5],[-2-mu_Z,-2+mu_Z],color='k',linewidth=2,alpha=0.7,label=None)
    plt.plot([14.5-mu_GR,14.5+mu_GR],[-2,-2],color='k',linewidth=2,alpha=0.7,label=None)
    cbar = plt.colorbar()
    cbar.set_label(r'Age [Gyr]', rotation=270, fontsize=15, labelpad=25)
    plt.xlabel(r'Galactic Radius [kpc]',fontsize=20)
    # plt.xlim(min(K2_New['Gal_Rad'])-0.2,15)#max(APK2['Gal_Rad'])+0.2)
    plt.ylabel(r'Z [kpc]',fontsize=20)
    plt.tick_params(labelsize=15)
    plt.subplots_adjust(left=0.13,top=0.93,bottom=0.15)
    plt.legend()
    # plt.ylim(-4,4)
    plt.xlim(min(Luca['Gal_Rad'])-0.2,max(APK2['Gal_Rad'])+0.2)
    # plt.show()
    f.savefig('GalR_Z_Paula.pdf', bbox_inches='tight')
    sys.exit()
    # pdf.savefig(f)
    # pdf.close()


    # K2_ar = K2_ar[(K2_ar['mass'] > 1.8) & (K2_ar['mass'] < 2.2)]
    # fig, ax = plt.subplots()
    # ax.scatter(AS.mass, AS.nmx)
    # ax.scatter(K2_ar.mass,K2_ar.nmx)
    #
    # fig, ax = plt.subplots()
    # # print(APK2.columns.values)
    # # ax.scatter(APK2['Teff'],np.log10(APK2['lum']),color='gray',alpha=0.2)
    # ax.scatter(AS['Teff'],np.log10(AS['gLumo']))
    # ax.scatter(K2_ar['Teff'],np.log10(K2_ar['gLumo']))
    # ax.invert_xaxis()
    # print(len(K2_ar)/len(AS),len(AS),len(K2_ar))
    # plt.show()
    # sys.exit()

    # APK2['gal_rad2'] = APK2['dist']*np.cos(APK2['Glat']*np.pi/180)*np.tan(APK2['Glon']*np.pi/180)*1e-3
    APK2['gal_rad3'] = np.sqrt(8250**2 + APK2['dist']**2 - 2*8250*APK2['dist']*np.cos(APK2['Glon']*np.pi/180))*1e-3
    APK2['Z2'] = APK2['dist']*np.sin(APK2['Glat'])*1e-3

    ''' Distance check plots '''
    # plt.figure()
    # plt.scatter(APK2['dist']*1e-3,APK2['gal_rad3'], label=r'Casagrande 2016')
    # plt.scatter(APK2['dist']*1e-3,APK2['Gal_Rad'],label=r'sqrt(X**2+Y**2)')
    # # plt.scatter(APK2['dist']*1e-3,APK2['gal_rad2'])
    # plt.xlabel(r'Distance [kpc]')
    # plt.ylabel(r'R$_{\rm{Gal}}$ [kpc]')
    # plt.legend()
    #
    # plt.figure()
    # plt.scatter(APK2['dist']*1e-3,APK2['dist']*np.sin(APK2['Glat']*np.pi/180)*1e-3,label=r'Casagrande 2016')
    # plt.scatter(APK2['dist']*1e-3,APK2['Z'], label=r'Original')
    # # plt.scatter(AS['dist']*1e-3,AS['gal_rad2'])
    # plt.xlabel(r'Distance [kpc]')
    # plt.ylabel(r'Z [kpc]')
    #
    # plt.show()

    ''' Alpha-rich young stars '''
    # K2_ar = AS[(AS['alpha'] > 0.1) & ((10**AS['logAge']/1e9) < 7)]
    # K2_ar = K2_ar[abs((K2_ar['mass']-K2_ar['m_seismo'])/K2_ar['mass']) < 0.5]
    # K2_ar = K2_ar.reset_index(drop=True)
    # # K2_ar.to_csv('/home/ben/K2_Poles/Mass_Distr_In/Young_Alpha_Rich',index=False)
    # plt.figure()
    # plt.scatter(K2_ar['Gal_Rad'],K2_ar['Z'])
    # plt.xlabel(r'Galactic Radius [kpc]',fontsize=20)
    # plt.ylabel(r'Z [kpc]',fontsize=20)
    # plt.tick_params(labelsize=15)
    # plt.subplots_adjust(left=0.13,top=0.93,bottom=0.15)
    # # plt.scatter(K2_ar['mass'],((K2_ar['mass']-K2_ar['m_seismo'])/K2_ar['mass']))
    # # plt.xlabel(r'Mass')
    # # plt.ylabel(r'(Param - Scaling)/Param')
    # plt.show()
    #
    # plt.figure()
    # plt.scatter(K2_New['nmx'],K2_New['dnu'],alpha=0.1,label=r'K2')
    # plt.scatter(K2_ar['nmx'],K2_ar['dnu'],alpha=0.5,label=r'K2 $\alpha$-RY')
    # plt.xlabel(r'$\nu_{\rm{max}}$ [$\mu$Hz]')
    # plt.ylabel(r'$\Delta\nu$ [$\mu$Hz]')
    # plt.legend()
    # plt.show()

    ''' Radial MDF - 1kpc bins '''
    # low3 = c_three[c_three['Gal_Rad'] < 7]
    # high3 = c_three[c_three['Gal_Rad'] > 7]
    # low6 = c_six[c_six['Gal_Rad'] < 7]
    # high6 = c_six[c_six['Gal_Rad'] > 7]
    # fig, (ax1,ax2) = plt.subplots(2)
    # ax1.hist(low3['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'R $<$ 7kpc',normed=1,linewidth=2)
    # ax1.hist(high3['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'R $>$ 7kpc',normed=1,linewidth=2)
    # ax1.legend()
    # ax1.set_xlabel(r'[Fe/H] - C3')
    # ax2.hist(low6['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'R $<$ 7kpc',normed=1,linewidth=2)
    # ax2.hist(high6['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'R $>$ 7kpc',normed=1,linewidth=2)
    # ax2.legend()
    # ax2.set_xlabel(r'[Fe/H] - C6')
    # plt.tight_layout()
    # plt.show()
