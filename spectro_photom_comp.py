''' Program for determining offsets between spectroscopic and photometric values '''

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import scipy.optimize as op
import scipy.odr.odrpack as odrpack
from sklearn import linear_model, datasets
import sys
import os

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams["font.family"] = "serif"

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

def ransac_fit(df,df1,param,label,f,ax,uncert):
    ''' Using RANSAC fitting algorithm to fit the trends observed between different
        spectroscopic pipelines.

    Data in:
        - DataFrames containing the parameters to be compared (should already be
          the same length)
        - Parameters to be compared (string)
        - Uncertainty if applicable (string)
        - Axes labels (string)
        - f, factor to extend plot limits to in the x-axis (float)
        - ax, subplot axes for plot to be associated with
    '''
    df['Diff'] = df[param[0]] - df1[param[1]]
    # df['sig'] = np.sqrt(df[uncert[0]]**2 + df1[uncert[1]]**2)
    x = (df1[param[1]])
    # X = X.values.reshape(-1,1)
    # y = (df['Diff'])
    Y = (df[param[0]])
    # print(df1.columns.values)
    # sys.exit()
    ''' Fit line using all data '''
    # lr = linear_model.LinearRegression()
    # lr.fit(X, y)

    ''' Robustly fit linear model with RANSAC algorithm '''
    # plt.figure()
    ax.set_xlim(df1[param[1]].min()-f, df1[param[1]].max()+f)
    # ax.set_xlabel(label[1])
    # ax.set_ylabel(label[0])
    # ax.scatter(X, y, color='yellowgreen', label='__nolegend__')
    lw = 2
    a = np.zeros((100,2))
    b = np.zeros((100,np.shape(df1)[0])) # Size of inlier mask
    df['pert'] = 0.
    df1['pert'] = 0.
    for i in range(100):
        for j in range(len(df1)):
            gauss1 = np.random.normal(df1[param[1]].iloc[j],df1[uncert[1]].iloc[j],50)
            value = gauss1[np.random.randint(len(gauss1),size=1)][0]
            df1['pert'].iloc[j] = value
            gauss2 = np.random.normal(df[param[0]].iloc[j],df[uncert[0]].iloc[j],50)
            value2 = gauss2[np.random.randint(len(gauss2),size=1)][0]
            df['pert'].iloc[j] = value2

        X = (df1['pert'])
        X = X.values.reshape(-1,1)
        y = (df['pert'])

        ransac = linear_model.RANSACRegressor(min_samples=7)
        ransac.fit(X, y)
        inlier_mask = ransac.inlier_mask_
        outlier_mask = np.logical_not(inlier_mask)

        ''' Predict data of estimated models '''
        # z = np.arange(X.min()-f, X.max()+f, 0.01)
        line_X = np.arange(x.min()-f, x.max()+f, 0.01)#[:, np.newaxis]
        line_X = line_X.reshape(-1,1)
        # line_y = lr.predict(line_X)
        line_y_ransac = ransac.predict(line_X)
        ''' Compare estimated coefficients '''
        a[i,0], a[i,1] = ransac.estimator_.coef_,ransac.estimator_.intercept_
        b[i,:] = inlier_mask
        # if i%10 == 0:
        ax.plot(line_X, a[i,0]*line_X + a[i,1], color='k', linewidth=3, label='__nolegend__', alpha=0.01)
            # plt.draw()
            # plt.waitforbuttonpress(0.1)

    ''' Make inliers boolean and find optimal fit using the median of the data '''
    b = b.astype('bool')
    c = np.median(a[:,0])
    medIdx = (np.abs(a[:,0] - c)).argmin()
    # print(a[medIdx,0],a[medIdx,1])
    print(np.std(df['Diff']))
    d = -a[medIdx,1]/a[medIdx,0]
    # print(np.median(a[:,0]),a[medIdx,0],b[medIdx])
    # df1['Teff_apo_corr'] = 0 - (a[medIdx,0]*df[param[0]] + a[medIdx,1])
    df['corr'] = df[param[0]] - (a[medIdx,0]*df[param[0]] + a[medIdx,1])
    # df1['corr'] = 0.
    # for i in range(len(df)):
    #     if (df[param[0]].iloc[i] > d) & ((a[medIdx,0]*df[param[0]].iloc[i] + a[medIdx,1]) > 0.):
    #         df1['corr'].iloc[i] = df1[param[1]].iloc[i] + (a[medIdx,0]*df[param[0]].iloc[i] + a[medIdx,1])
    #     elif (df[param[0]].iloc[i] > d) & ((a[medIdx,0]*df[param[0]].iloc[i] + a[medIdx,1]) < 0.):
    #         df1['corr'].iloc[i] = df1[param[1]].iloc[i] - (a[medIdx,0]*df[param[0]].iloc[i] + a[medIdx,1])
    #     elif (df[param[0]].iloc[i] < d) & ((a[medIdx,0]*df[param[0]].iloc[i] + a[medIdx,1]) > 0.):
    #         df1['corr'].iloc[i] = df1[param[1]].iloc[i] - (a[medIdx,0]*df[param[0]].iloc[i] + a[medIdx,1])
    #     elif (df[param[0]].iloc[i] < d) & ((a[medIdx,0]*df[param[0]].iloc[i] + a[medIdx,1]) < 0.):
    #         df1['corr'].iloc[i] = df1[param[1]].iloc[i] + (a[medIdx,0]*df[param[0]].iloc[i] + a[medIdx,1])

    # plt.figure()
    # plt.scatter(df[param[0]],df1['corr'])
    # plt.scatter(df[param[0]],df['corr']-df1[param[1]])
    ''' Plot the inliers, outliers and optimal fit to the data '''
    # plt.figure()
    # lw = 2
    ax.plot(line_X,line_X,color='r', linewidth=3,alpha=0.5,linestyle='--',label='__nolegend__')
    ax.plot(line_X, a[medIdx,0]*line_X + a[medIdx,1], color='orange', linewidth=3, label=r'%.5s$*x$ + %.5s'%(a[medIdx,0],a[medIdx,1]))
    # plt.axhline(y=0., color='r', linestyle='-')
    # plt.axvline(x=d, color='r', linestyle='-')

    ax.errorbar(x, Y, xerr=df1[uncert[1]], yerr=df[uncert[0]], color='blue', marker='.', label='__nolegend__',fmt='o',alpha=0.25)

    # plt.errorbar(X[b[medIdx]], y[b[medIdx]], yerr=df['sig_Teff'][b[medIdx]], color='yellowgreen', marker='.', label='Inliers',fmt='o',alpha=0.5)
    ax.scatter(x, Y, color='blue', label='__nolegend__')
    # plt.scatter(X[b[medIdx]], y[b[medIdx]], color='yellowgreen', label='Inliers')
    # plt.xlabel(label[0])
    # plt.ylabel(label[1])
    # plt.xlim(X.min()-f, X.max()+f)
    ax.legend(frameon=False)
    # plt.show()
    # df['outlier_flag'] = 1.
    # df['outlier_flag'][b[medIdx]] = 0.
    # print(df['outlier_flag'])
    # plt.savefig('/home/bmr135/spectro_comp/C3_RAVE_GES_'+param[1]+'.png')

def odr_fit(df,df1,param,label,f1,ax,uncert):
    ''' Using orthogonal distance regression alogrithm to perform the fit '''
    def f(Par,z):
        return Par[0]*z + Par[1]
    mpar, cpar, empar, ecpar = [], [], [], []
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(df1[param[1]],df[param[0]],sx=df1[uncert[1]],sy=df[uncert[0]])
    myodr = odrpack.ODR(mydata, linear, beta0=[1., 0.],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    # print(np.sqrt(np.diag(myoutput.cov_beta * myoutput.res_var)),[empar,ecpar])
    Z = np.sqrt(np.diag(myoutput.cov_beta * myoutput.res_var))
    # myoutput.pprint()
    lw = 2
    means = np.array(myoutput.beta)
    cov = np.array(myoutput.cov_beta)
    np.diag(cov)**0.5
    print(means, cov)
    x = np.linspace(df1[param[1]].min()-f1, df1[param[1]].max()+f1, 1000)
    aa, bb = np.random.multivariate_normal(means, cov, size=len(x)).T # Take samples of a random distribution centred on the best fit,
                                                                      # using the covariance matrix to vary the values.
    tmp = [np.std(aa * i + bb) for i in x] # Take std dev of the distribution of possible models to provide confidence interval.
    a = means[0]
    b = means[1]
    y = a * x + b

    # fig = plt.figure()
    # fig, ax = plt.subplots()
    #ax.errorbar(x, y, yerr=yerr, fmt='ko', label='Data')
    ax.plot(x, a * x + b, 'r-', label=r'%.5s$*x$ + %.5s'%(a,b))
    ax.fill_between(x, a * x + b - tmp, a * x + b + tmp, alpha=0.25, label='Confidence interval',color='gray')
    ax.legend(loc=2)
    # plt.show()
    # sys.exit()

    line_X = np.arange(min(df1[param[1]])-f1, max(df1[param[1]])+f1, 0.01)#[:, np.newaxis]
    ax.plot(line_X,line_X,color='r', linewidth=3,alpha=0.5,linestyle='--',label='__nolegend__')
    ax.errorbar(df1[param[1]], df[param[0]], xerr=df1[uncert[1]], yerr=df[uncert[0]], color='blue', marker='.', label='__nolegend__',fmt='o',alpha=0.25)
    ax.scatter(df1[param[1]], df[param[0]], color='blue', label='__nolegend__')
    ax.legend(frameon=False)



if __name__ == "__main__":
    ''' Read in Data '''
    # ext = '/media/bmr135/SAMSUNG/GA/K2Poles'
    ext = '/media/ben/SAMSUNG/GA/K2Poles'
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
    a = R3[R3['Teff_RAVE'] < 4400]
    A3['FE_H_ERR'] = A3['FE_H_ERR']*20
    # print(a[['Teff_RAVE','logg']])
    R3 = R3[R3['Teff_RAVE'] > 4000]
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
    # print(len(ag),len(gr),len(ar))
    # print(A3['FE_H_ERR'])
    # print(G.columns.values)
    # print(R3.columns.values)
    # sys.exit()

    # print(rg[['Teff_RAVE','logg_RAVE','logg','radius']])
    # print(gr[['TEFF','LOGG','logg','radius']])
    # sys.exit()

    ''' Compare Teff, [Fe/H], logg of given datasets '''

    # ''' APOGEE vs EPIC - C3 '''
    # # ransac_fit(a3,ac3,['TEFF','Teff'],[r'$T_{\rm{eff}}$ - APOGEE, C3',r'$\Delta(T_{\rm{eff}})_{APO-EPIC}$'],10)#,['TEFF_ERR','sig_Teff'])
    # # ransac_fit(a3,ac3,['LOGG','logg'],[r'log$_{10}$(g) - APOGEE, C3',r'$\Delta($log$_{10}$(g)$)_{APO-EPIC}$'],0.1)#,['LOGG_ERR','sig_logg'])
    # # ransac_fit(a3,ac3,['FE_H','[Fe/H]'],[r'[Fe/H] - APOGEE, C3',r'$\Delta($[Fe/H]$)_{APO-EPIC}$'],0.1)#,['FE_H_ERR','sig_feh'])
    # # plt.show()
    #
    # ''' APOGEE vs EPIC - C6 '''
    # # ransac_fit(a6,ac6,['TEFF','Teff'],[r'$T_{\rm{eff}}$ - APOGEE, C6',r'$\Delta(T_{\rm{eff}})_{APO-EPIC}$'],10)
    # # ransac_fit(a6,ac6,['LOGG','logg'],[r'log$_{10}$(g) - APOGEE, C6',r'$\Delta($log$_{10}$(g)$)_{APO-EPIC}$'],0.1)
    # # ransac_fit(a6,ac6,['FE_H','[Fe/H]'],[r'[Fe/H] - APOGEE, C6',r'$\Delta($[Fe/H]$)_{APO-EPIC}$'],0.1)
    # # plt.show()
    #
    # ''' Gaia-ESO vs EPIC - C3 '''
    # # ransac_fit(g,gc,['TEFF','Teff'],[r'$T_{\rm{eff}}$ - Gaia-ESO, C3',r'$\Delta(T_{\rm{eff}})_{GES-EPIC}$'],10)
    # # ransac_fit(g,gc,['LOGG','logg'],[r'log$_{10}$(g) - Gaia-ESO, C3',r'$\Delta($log$_{10}$(g)$)_{GES-EPIC}$'],0.1)
    # # ransac_fit(g,gc,['FEH','[Fe/H]'],[r'[Fe/H] - Gaia-ESO, C3',r'$\Delta($[Fe/H]$)_{GES-EPIC}$'],0.1)
    # # plt.show()
    #
    # ''' RAVE vs EPIC - C3 '''
    # # ransac_fit(r3,rc3,['Teff_RAVE','Teff'],[r'$T_{\rm{eff}}$ - RAVE, C3',r'$\Delta(T_{\rm{eff}})_{RAVE-EPIC}$'],10)
    # # ransac_fit(r3,rc3,['logg_RAVE','logg'],[r'log$_{10}$(g) - RAVE, C3',r'$\Delta($log$_{10}$(g)$)_{RAVE-EPIC}$'],0.1)
    # # ransac_fit(r3,rc3,['[Fe/H]_RAVE','[Fe/H]'],[r'[Fe/H] - RAVE, C3',r'$\Delta($[Fe/H]$)_{RAVE-EPIC}$'],0.1)
    # # plt.show()
    #
    # ''' RAVE vs EPIC - C6 '''
    # # ransac_fit(r6,rc6,['Teff_RAVE','Teff'],[r'$T_{\rm{eff}}$ - RAVE, C6',r'$\Delta(T_{\rm{eff}})_{RAVE-EPIC}$'],10)
    # # ransac_fit(r6,rc6,['logg_RAVE','logg'],[r'log$_{10}$(g) - RAVE, C6',r'$\Delta($log$_{10}$(g)$)_{RAVE-EPIC}$'],0.1)
    # # ransac_fit(r6,rc6,['[Fe/H]_RAVE','[Fe/H]'],[r'[Fe/H] - RAVE, C6',r'$\Delta($[Fe/H]$)_{RAVE-EPIC}$'],0.1)
    # # plt.show()
    #
    # ''' LAMOST vs EPIC - C6 '''
    # # ransac_fit(l6,lc6,['teff_L','Teff'],[r'$T_{\rm{eff}}$ - LAMOST, C6',r'$\Delta(T_{\rm{eff}})_{LAMOST-EPIC}$'],10)
    # # ransac_fit(l6,lc6,['logg_L','logg'],[r'log$_{10}$(g) - LAMOST, C6',r'$\Delta($log$_{10}$(g)$)_{LAMOST-EPIC}$'],0.1)
    # # ransac_fit(l6,lc6,['feh_L','[Fe/H]'],[r'[Fe/H] - LAMOST, C6',r'$\Delta($[Fe/H]$)_{LAMOST-EPIC}$'],0.1)
    # # plt.show()
    #
    # ''' APOGEE vs RAVE - C3 '''
    # # if len(ar3) > 10:
    # #     print('wooo')
    # #     ransac_fit(ar3,ra3,['TEFF','Teff_RAVE'],[r'$T_{\rm{eff}}$ - APOGEE, C3',r'$\Delta(T_{\rm{eff}})_{APO-RAVE}$'],10)
    # #     ransac_fit(ar3,ra3,['LOGG','logg_RAVE'],[r'log$_{10}$(g) - APOGEE, C3',r'$\Delta($log$_{10}$(g)$)_{APO-RAVE}$'],0.1)
    # #     ransac_fit(ar3,ra3,['FE_H','[Fe/H]_RAVE'],[r'[Fe/H] - APOGEE, C3',r'$\Delta($[Fe/H]$)_{APO-RAVE}$'],0.1)
    # # plt.show()
    #
    # ''' APOGEE vs RAVE - C6 '''
    # # if len(ar) > 10:
    # #     ransac_fit(ar,ra,['TEFF','Teff_RAVE'],[r'$T_{\rm{eff}}$ - APOGEE, C6',r'$\Delta(T_{\rm{eff}})_{APO-RAVE}$'],10)
    # #     plt.show()
    # #     sys.exit()
    # #     ransac_fit(ar,ra,['LOGG','logg_RAVE'],[r'log$_{10}$(g) - APOGEE, C6',r'$\Delta($log$_{10}$(g)$)_{APO-RAVE}$'],0.1)
    # #     # plt.show()
    # #     ransac_fit(ar,ra,['FE_H','[Fe/H]_RAVE'],[r'[Fe/H] - APOGEE, C6',r'$\Delta($[Fe/H]$)_{APO-RAVE}$'],0.1)
    # #     plt.show()

    ''' APOGEE vs Gaia-ESO - C3 '''
    fig, ((ax,ax1,ax2),(ax3,ax4,ax5)) = plt.subplots(2,3,sharex='col',figsize=(15,10))
    ax.set_ylabel(r'$T_{\rm{eff}}$ - APOGEE, C3',fontsize=15)
    ax1.set_ylabel(r'log$_{10}$(g) - APOGEE, C3',fontsize=15)
    ax2.set_ylabel(r'[Fe/H] - APOGEE, C3',fontsize=15)
    ax3.set_ylabel(r'$T_{\rm{eff}}$ - RAVE, C3',fontsize=15)
    ax4.set_ylabel(r'log$_{10}$(g) - RAVE, C3',fontsize=15)
    ax5.set_ylabel(r'[Fe/H] - RAVE, C3',fontsize=15)
    ax3.set_xlabel(r'$T_{\rm{eff}}$ - GES, C3',fontsize=15)
    ax4.set_xlabel(r'log$_{10}$(g) - GES, C3',fontsize=15)
    ax5.set_xlabel(r'[Fe/H] - GES, C3',fontsize=15)

    if len(ag) > 10:
        # ransac_fit(ag,ga,['TEFF','TEFF'],[r'$T_{\rm{eff}}$ - APOGEE, C3',r'$T_{\rm{eff}}$ - GES, C3'],500,ax,['TEFF_ERR','sig_Teff'])#r'$\Delta(T_{\rm{eff}})_{APO-GES}$'],10)
        odr_fit(ag,ga,['TEFF','TEFF'],[r'$T_{\rm{eff}}$ - APOGEE, C3',r'$T_{\rm{eff}}$ - GES, C3'],500,ax,['TEFF_ERR','sig_TEFF'])#r'$\Delta(T_{\rm{eff}})_{APO-GES}$'],10)
        # plt.show()
        # ransac_fit(ag,ga,['LOGG','LOGG'],[r'log$_{10}$(g) - APOGEE, C3',r'log$_{10}$(g) - GES, C3'],0.45,ax1,['LOGG_ERR','sig_logg'])#r'$\Delta($log$_{10}$(g)$)_{APO-GES}$'],0.1)
        odr_fit(ag,ga,['LOGG','LOGG'],[r'log$_{10}$(g) - APOGEE, C3',r'log$_{10}$(g) - GES, C3'],0.45,ax1,['LOGG_ERR','sig_LOGG'])#r'$\Delta($log$_{10}$(g)$)_{APO-GES}$'],0.1)
        # plt.show()
        # ransac_fit(ag,ga,['FE_H','FEH'],[r'[Fe/H] - APOGEE, C3',r'[Fe/H] - GES, C3'],1.,ax2,['FE_H_ERR','sig_feh'])#r'$\Delta($[Fe/H]$)_{APO-GES}$'],0.1)
        odr_fit(ag,ga,['FE_H','FEH'],[r'[Fe/H] - APOGEE, C3',r'[Fe/H] - GES, C3'],1.,ax2,['FE_H_ERR','sig_FEH'])#r'$\Delta($[Fe/H]$)_{APO-GES}$'],0.1)
        # plt.show()
        # sys.exit()

    ''' RAVE vs Gaia-ESO - C3 '''
    if len(rg) > 10:
        # ransac_fit(rg,gr,['Teff_RAVE','TEFF'],[r'$T_{\rm{eff}}$ - RAVE, C3',r'$T_{\rm{eff}}$ - GES, C3'],500,ax3,['sig_Teff','sig_Teff'])#,r'$\Delta(T_{\rm{eff}})_{RAVE-GES}$'],10)
        odr_fit(rg,gr,['Teff_RAVE','TEFF'],[r'$T_{\rm{eff}}$ - RAVE, C3',r'$T_{\rm{eff}}$ - GES, C3'],500,ax3,['sig_Teff','sig_TEFF'])#,r'$\Delta(T_{\rm{eff}})_{RAVE-GES}$'],10)
        # plt.show()
        # ransac_fit(rg,gr,['logg_RAVE','LOGG'],[r'log$_{10}$(g) - RAVE, C3',r'log$_{10}$(g) - GES, C3'],0.45,ax4,['sig_logg','sig_logg'])#,r'$\Delta($log$_{10}$(g)$)_{RAVE-GES}$'],0.1)
        odr_fit(rg,gr,['logg_RAVE','LOGG'],[r'log$_{10}$(g) - RAVE, C3',r'log$_{10}$(g) - GES, C3'],0.45,ax4,['sig_logg','sig_LOGG'])#,r'$\Delta($log$_{10}$(g)$)_{RAVE-GES}$'],0.1)
        # plt.show()
        # ransac_fit(rg,gr,['[Fe/H]_RAVE','FEH'],[r'[Fe/H] - RAVE, C3',r'[Fe/H] - GES, C3'],1.,ax5,['sig_feh','sig_feh'])#,r'$\Delta($[Fe/H]$)_{RAVE-GES}$'],0.1)
        odr_fit(rg,gr,['[Fe/H]_RAVE','FEH'],[r'[Fe/H] - RAVE, C3',r'[Fe/H] - GES, C3'],1.,ax5,['sig_feh','sig_FEH'])#,r'$\Delta($[Fe/H]$)_{RAVE-GES}$'],0.1)
        # plt.show()

    ax3.set_xlim(4400, 5200)
    ax.set_ylim(4400, 5200)
    ax3.set_ylim(4400, 5200)
    ax4.set_xlim(2.1,3.4)
    ax4.set_ylim(2.1,3.4)
    ax1.set_ylim(2.1,3.4)
    ax5.set_xlim(-1.2, 0.6)
    ax5.set_ylim(-1.2, 0.6)
    ax2.set_ylim(-1.2, 0.6)
    plt.subplots_adjust(hspace=0.07,wspace=0.27)
    fig.savefig('spec_comps.pdf',bbox_inches='tight')
    plt.show()

    ''' APOGEE vs LAMOST - C6 '''
    # if len(al) > 2:
    #     ransac_fit(al,la,['TEFF','teff_L'],[r'$T_{\rm{eff}}$ - APOGEE, C6',r'$\Delta(T_{\rm{eff}})_{APO-LAMOST}$'],10)
    #     ransac_fit(al,la,['LOGG','logg_L'],[r'log$_{10}$(g) - APOGEE, C6',r'$\Delta($log$_{10}$(g)$)_{APO-LAMOST}$'],0.1)
    #     ransac_fit(al,la,['FE_H','feh_L'],[r'[Fe/H] - APOGEE, C6',r'$\Delta($[Fe/H]$)_{APO-LAMOST}$'],0.1)
