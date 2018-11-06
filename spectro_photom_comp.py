''' Program for determining offsets between spectroscopic and photometric values '''

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.optimize as op
import scipy.odr.odrpack as odrpack
from sklearn import linear_model, datasets
import sys
import os

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

def truncate(a,b):
    ''' Truncate datasets to find overlapping stars.

        a = dataset 1
        b = dataset to find values of a in
        x = placeholder for EPICS in a to find in b
        y = EPICS in turncated b to act as indexing tool for selecting values
            from a
    '''
    x = pd.DataFrame()
    x['EPIC'] = a['EPIC']
    b = pd.merge(x,b,how='inner',on=['EPIC'])
    b.reset_index(drop=True)
    y = pd.DataFrame()
    y['EPIC'] = b['EPIC']
    a = pd.merge(y,a,how='inner',on=['EPIC'])
    a.reset_index(drop=True)
    # print(len(a),len(b))

    return a,b

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
    # df['sig'] = np.sqrt(df[uncert[0]]**2 + df1[uncert[1]]**2)
    X = (df[param[0]])
    X = X.values.reshape(-1,1)
    y = (df['Diff'])

    ''' Fit line using all data '''
    lr = linear_model.LinearRegression()
    lr.fit(X, y)

    ''' Robustly fit linear model with RANSAC algorithm '''
    a = np.zeros((100,2))
    b = np.zeros((100,np.shape(df)[0])) # Size of inlier mask
    for i in range(100):
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
    # df1['Teff_apo_corr'] = 0 - (a[medIdx,0]*df[param[0]] + a[medIdx,1])
    df1['corr'] = df1[param[1]] + (a[medIdx,0]*df[param[0]] + a[medIdx,1])

    plt.figure()
    # plt.scatter(df[param[0]],df1['corr'])
    plt.scatter(df[param[0]],df1['corr']-df[param[0]])

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
    print(df['outlier_flag'])
    # plt.savefig('/home/bmr135/spectro_comp/C3_RAVE_GES_'+param[1]+'.png')

if __name__ == "__main__":
    ''' Read in Data '''
    ext = '/media/bmr135/SAMSUNG/GA/K2Poles'
    A3 = pd.read_csv(ext+'/APO_LAMOST/APOGEE_full_C3.csv')
    A3 = A3.dropna(subset=['TEFF','LOGG','FE_H'])
    A6 = pd.read_csv(ext+'/APO_LAMOST/APOGEE_full_C6.csv')
    A6 = A6.dropna(subset=['TEFF','LOGG','FE_H'])
    G = pd.read_csv(ext+'/Gaia_ESO/GES_full.csv')
    G = G.dropna(subset=['TEFF','LOGG','FEH'])
    R3 = pd.read_csv(ext+'/RAVE_C3.csv')
    R3 = R3.dropna(subset=['Teff_RAVE','logg_RAVE','[Fe/H]_RAVE'])
    R3 = R3[R3['logg_RAVE'] > 0]
    R6 = pd.read_csv(ext+'/RAVE_C6.csv')
    R6 = R6.dropna(subset=['Teff_RAVE','logg_RAVE','[Fe/H]_RAVE'])
    R6 = R6[R6['logg_RAVE'] > 0]
    L6 = pd.read_csv(ext+'/APO_LAMOST/LAMOST_full_C6.csv')
    L6 = L6.dropna(subset=['teff_L','logg_L','feh_L'])
    C3 = pd.read_csv(ext+'/matlab_in/C3_02112018.csv')
    C3 = C3.dropna(subset=['Teff','logg','[Fe/H]'])
    C6 = pd.read_csv(ext+'/matlab_in/C6_02112018.csv')
    C6 = C6.dropna(subset=['Teff','logg','[Fe/H]'])

    ''' Combine datasets prior to comparisons '''
    a3, ac3 = truncate(A3,C3)
    a6, ac6 = truncate(A6,C6)
    g, gc = truncate(G,C3)
    r3, rc3 = truncate(R3,C3)
    r6, rc6 = truncate(R6,C6)
    l6, lc6 = truncate(L6,C6)
    ar3, ra3 = truncate(A3,R3)
    ag, ga = truncate(A3,G)
    gr, rg = truncate(G,R3)
    ar, ra = truncate(A6,R6)
    al, la = truncate(A6,L6)

    # print(R3.columns.values)

    ''' Compare Teff, [Fe/H], logg of given datasets '''

    ''' APOGEE vs EPIC - C3 '''
    # ransac_fit(a3,ac3,['TEFF','Teff'],[r'$T_{\rm{eff}}$ - APOGEE, C3',r'$\Delta(T_{\rm{eff}})_{APO-EPIC}$'],10)#,['TEFF_ERR','sig_Teff'])
    # ransac_fit(a3,ac3,['LOGG','logg'],[r'log$_{10}$(g) - APOGEE, C3',r'$\Delta($log$_{10}$(g)$)_{APO-EPIC}$'],0.1)#,['LOGG_ERR','sig_logg'])
    # ransac_fit(a3,ac3,['FE_H','[Fe/H]'],[r'[Fe/H] - APOGEE, C3',r'$\Delta($[Fe/H]$)_{APO-EPIC}$'],0.1)#,['FE_H_ERR','sig_feh'])
    # plt.show()

    ''' APOGEE vs EPIC - C6 '''
    # ransac_fit(a6,ac6,['TEFF','Teff'],[r'$T_{\rm{eff}}$ - APOGEE, C6',r'$\Delta(T_{\rm{eff}})_{APO-EPIC}$'],10)
    # ransac_fit(a6,ac6,['LOGG','logg'],[r'log$_{10}$(g) - APOGEE, C6',r'$\Delta($log$_{10}$(g)$)_{APO-EPIC}$'],0.1)
    # ransac_fit(a6,ac6,['FE_H','[Fe/H]'],[r'[Fe/H] - APOGEE, C6',r'$\Delta($[Fe/H]$)_{APO-EPIC}$'],0.1)
    # plt.show()

    ''' Gaia-ESO vs EPIC - C3 '''
    # ransac_fit(g,gc,['TEFF','Teff'],[r'$T_{\rm{eff}}$ - Gaia-ESO, C3',r'$\Delta(T_{\rm{eff}})_{GES-EPIC}$'],10)
    # ransac_fit(g,gc,['LOGG','logg'],[r'log$_{10}$(g) - Gaia-ESO, C3',r'$\Delta($log$_{10}$(g)$)_{GES-EPIC}$'],0.1)
    # ransac_fit(g,gc,['FEH','[Fe/H]'],[r'[Fe/H] - Gaia-ESO, C3',r'$\Delta($[Fe/H]$)_{GES-EPIC}$'],0.1)
    # plt.show()

    ''' RAVE vs EPIC - C3 '''
    # ransac_fit(r3,rc3,['Teff_RAVE','Teff'],[r'$T_{\rm{eff}}$ - RAVE, C3',r'$\Delta(T_{\rm{eff}})_{RAVE-EPIC}$'],10)
    # ransac_fit(r3,rc3,['logg_RAVE','logg'],[r'log$_{10}$(g) - RAVE, C3',r'$\Delta($log$_{10}$(g)$)_{RAVE-EPIC}$'],0.1)
    # ransac_fit(r3,rc3,['[Fe/H]_RAVE','[Fe/H]'],[r'[Fe/H] - RAVE, C3',r'$\Delta($[Fe/H]$)_{RAVE-EPIC}$'],0.1)
    # plt.show()

    ''' RAVE vs EPIC - C6 '''
    # ransac_fit(r6,rc6,['Teff_RAVE','Teff'],[r'$T_{\rm{eff}}$ - RAVE, C6',r'$\Delta(T_{\rm{eff}})_{RAVE-EPIC}$'],10)
    # ransac_fit(r6,rc6,['logg_RAVE','logg'],[r'log$_{10}$(g) - RAVE, C6',r'$\Delta($log$_{10}$(g)$)_{RAVE-EPIC}$'],0.1)
    # ransac_fit(r6,rc6,['[Fe/H]_RAVE','[Fe/H]'],[r'[Fe/H] - RAVE, C6',r'$\Delta($[Fe/H]$)_{RAVE-EPIC}$'],0.1)
    # plt.show()

    ''' LAMOST vs EPIC - C6 '''
    # ransac_fit(l6,lc6,['teff_L','Teff'],[r'$T_{\rm{eff}}$ - LAMOST, C6',r'$\Delta(T_{\rm{eff}})_{LAMOST-EPIC}$'],10)
    # ransac_fit(l6,lc6,['logg_L','logg'],[r'log$_{10}$(g) - LAMOST, C6',r'$\Delta($log$_{10}$(g)$)_{LAMOST-EPIC}$'],0.1)
    # ransac_fit(l6,lc6,['feh_L','[Fe/H]'],[r'[Fe/H] - LAMOST, C6',r'$\Delta($[Fe/H]$)_{LAMOST-EPIC}$'],0.1)
    # plt.show()

    ''' APOGEE vs RAVE - C3 '''
    # if len(ar3) > 2:
    #     ransac_fit(ar3,ra3,['TEFF','Teff_RAVE'],[r'$T_{\rm{eff}}$ - APOGEE, C3',r'$\Delta(T_{\rm{eff}})_{APO-RAVE}$'],10)
    #     ransac_fit(ar3,ra3,['LOGG','logg_RAVE'],[r'log$_{10}$(g) - APOGEE, C3',r'$\Delta($log$_{10}$(g)$)_{APO-RAVE}$'],0.1)
    #     ransac_fit(ar3,ra3,['FE_H','[Fe/H]_RAVE'],[r'[Fe/H] - APOGEE, C3',r'$\Delta($[Fe/H]$)_{APO-RAVE}$'],0.1)
    # plt.show()

    ''' APOGEE vs RAVE - C6 '''
    if len(ar) > 2:
        ransac_fit(ar,ra,['TEFF','Teff_RAVE'],[r'$T_{\rm{eff}}$ - APOGEE, C6',r'$\Delta(T_{\rm{eff}})_{APO-RAVE}$'],10)
        plt.show()
        ransac_fit(ar,ra,['LOGG','logg_RAVE'],[r'log$_{10}$(g) - APOGEE, C6',r'$\Delta($log$_{10}$(g)$)_{APO-RAVE}$'],0.1)
        plt.show()
        ransac_fit(ar,ra,['FE_H','[Fe/H]_RAVE'],[r'[Fe/H] - APOGEE, C6',r'$\Delta($[Fe/H]$)_{APO-RAVE}$'],0.1)
        plt.show()

    ''' APOGEE vs Gaia-ESO - C3 '''
    if len(ag) > 2:
        ransac_fit(ag,ga,['TEFF','TEFF'],[r'$T_{\rm{eff}}$ - APOGEE, C3',r'$\Delta(T_{\rm{eff}})_{APO-GES}$'],10)
        plt.show()
        ransac_fit(ag,ga,['LOGG','LOGG'],[r'log$_{10}$(g) - APOGEE, C3',r'$\Delta($log$_{10}$(g)$)_{APO-GES}$'],0.1)
        plt.show()
        ransac_fit(ag,ga,['FE_H','FEH'],[r'[Fe/H] - APOGEE, C3',r'$\Delta($[Fe/H]$)_{APO-GES}$'],0.1)
        plt.show()

    ''' RAVE vs Gaia-ESO - C3 '''
    # if len(rg) > 2:
    #     ransac_fit(rg,gr,['Teff_RAVE','TEFF'],[r'$T_{\rm{eff}}$ - RAVE, C3',r'$\Delta(T_{\rm{eff}})_{RAVE-GES}$'],10)
    #     ransac_fit(rg,gr,['logg_RAVE','LOGG'],[r'log$_{10}$(g) - RAVE, C3',r'$\Delta($log$_{10}$(g)$)_{RAVE-GES}$'],0.1)
    #     ransac_fit(rg,gr,['[Fe/H]_RAVE','FEH'],[r'[Fe/H] - RAVE, C3',r'$\Delta($[Fe/H]$)_{RAVE-GES}$'],0.1)
    # plt.show()

    ''' APOGEE vs LAMOST - C6 '''
    # if len(al) > 2:
    #     ransac_fit(al,la,['TEFF','teff_L'],[r'$T_{\rm{eff}}$ - APOGEE, C6',r'$\Delta(T_{\rm{eff}})_{APO-LAMOST}$'],10)
    #     ransac_fit(al,la,['LOGG','logg_L'],[r'log$_{10}$(g) - APOGEE, C6',r'$\Delta($log$_{10}$(g)$)_{APO-LAMOST}$'],0.1)
    #     ransac_fit(al,la,['FE_H','feh_L'],[r'[Fe/H] - APOGEE, C6',r'$\Delta($[Fe/H]$)_{APO-LAMOST}$'],0.1)
    # plt.show()
