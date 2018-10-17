''' Program for determining offsets between spectroscopic and photometric values '''

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
# import emcee
# import corner
import numpy as np
import scipy.optimize as op
import scipy.odr.odrpack as odrpack
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
    print(len(a),len(b))

    return a,b

def comp_scat(A,B,teff,logg,feh,title,xsub,ysub,fit,save):
    ''' comparative scatter plots between surveys '''
    fig, ax = plt.subplots(3)
    plt.suptitle(title)
    if fit == 0:
        dt = fitting0(A,B,teff,'teff')
        dg = fitting0(A,B,logg,'logg')
        df = fitting0(A,B,feh,'feh')
        print(dt,dg,df)
    elif fit == 1:
        dt = fitting(A,B,teff,'teff')
        dg = fitting(A,B,logg,'logg')
        df = fitting(A,B,feh,'feh')
        print(dt,dg,df)

    ax[0].scatter(A[teff[0]],(A[teff[0]]-B[teff[1]]))
    ax[0].plot([min(A[teff[0]])-20,max(A[teff[0]])+20],[0,0],color='k')
    a = np.linspace(min(A[teff[0]])-20,max(A[teff[0]])+20,200)
    ax[0].plot(a,dt[0]*a+dt[1],color='m')
    ax[0].set_xlabel(r'T$_{\rm{eff},%s}$ [K]'%(xsub))
    ax[0].set_xlim(min(A[teff[0]])-20,max(A[teff[0]])+20)
    ax[0].set_ylabel(r'$\Delta$T$_{\rm{eff,%s}}$ [K]'%(ysub))

    ax[1].scatter(A[logg[0]],(A[logg[0]]-B[logg[1]]))
    ax[1].plot([min(A[logg[0]])-0.1,max(A[logg[0]])+0.1],[0,0],color='k')
    b = np.linspace(min(A[logg[0]])-0.1,max(A[logg[0]])+0.1,200)
    ax[1].plot(b,dg[0]*b+dg[1],color='m')
    ax[1].set_xlabel(r'log$_{10}$(g)$_{\rm{%s}}$'%(xsub))
    ax[1].set_xlim(min(A[logg[0]])-0.1,max(A[logg[0]])+0.1)
    ax[1].set_ylabel(r'$\Delta$log$_{10}$(g)$_{%s}$'%(ysub))

    ax[2].scatter(A[feh[0]],(A[feh[0]]-B[feh[1]]))
    ax[2].plot([min(A[feh[0]])-0.1,max(A[feh[0]])+0.1],[0,0],color='k')
    c = np.linspace(min(A[feh[0]])-0.1,max(A[feh[0]])+0.1,200)
    ax[2].plot(c,df[0]*c+df[1],color='m')
    ax[2].set_xlabel(r'[Fe/H]$_{\rm{%s}}$'%(xsub))
    ax[2].set_xlim(min(A[feh[0]])-0.1,max(A[feh[0]])+0.1)
    ax[2].set_ylabel(r'$\Delta$[Fe/H]$_{\rm{%s}}$'%(ysub))

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    plt.savefig('/media/bmr135/SAMSUNG/GA/K2Poles/Spec_Comp/'+save+'.png')

def fitting(df1,df2,param,param_str):
    ''' Least squares fitting approach to calculating offsets betrween params '''
    def f(Par,x):
        return Par[0]*x + Par[1]
    alt = 'alt_' + param_str
    df2[alt] = df1[param[0]] - df2[param[1]]
    # print(np.mean(df2[alt]),np.median(df2[alt]), len(df2), len(df1), df2[alt])
    df1['null'] = 0
    mpar, cpar, empar, ecpar = [], [], [], []
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(df1[param[0]], df2[alt])#, sy=df2['feh_err'])
    myodr = odrpack.ODR(mydata, linear, beta0=[0., 0.],maxit=100000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    # myoutput.pprint()

    return mpar,cpar,ecpar

def fitting0(df1,df2,param,param_str):
    ''' Least squares fitting approach to calculating offsets betrween params '''
    def f(Par,x):
        return Par[0]*x + Par[1]
    alt = 'alt_' + param_str
    df2[alt] = df1[param[0]] - df2[param[1]]
    df1['null'] = 0
    mpar, cpar, empar, ecpar = [], [], [], []
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(df1['null'], df2[alt])#, sy=df2['feh_err'])
    myodr = odrpack.ODR(mydata, linear, beta0=[0., 0.],maxit=100000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    # myoutput.pprint()

    return mpar,cpar,ecpar


if __name__ == "__main__":
    ''' Read in Data '''
    # ext = '/home/bmr135/GA/K2Poles'
    ext = '/media/ben/SAMSUNG/GA/K2Poles'
    APO = pd.read_csv(ext+'/matlab_in/APOGEE_full.csv')
    APO = APO.dropna(subset=['TEFF','LOGG','FE_H'])
    G = pd.read_csv(ext+'/Gaia_ESO/GES_full.csv')
    G = G.dropna(subset=['TEFF','LOGG','FEH'])
    R3 = pd.read_csv(ext+'/matlab_in/RC3.csv')
    R3 = R3.dropna(subset=['Teff_RAVE','logg_RAVE','[Fe/H]_RAVE'])
    R3 = R3[R3['logg_RAVE'] > 0]
    R6 = pd.read_csv(ext+'/matlab_in/RC6.csv')
    R6 = R6.dropna(subset=['Teff_RAVE','logg_RAVE','[Fe/H]_RAVE'])
    R6 = R6[R6['logg_RAVE'] > 0]
    # L3 = pd.read_csv(ext+'/matlab_in/LAMOST3_multidet.csv')
    # L3 = L3.dropna(subset=['teff_L','logg_L','feh_L'])
    L6 = pd.read_csv(ext+'/matlab_in/LAMOST6_multidet.csv')
    L6 = L6.dropna(subset=['teff_L','logg_L','feh_L'])
    # LR6 = pd.read_csv(ext+'/matlab_in/LAMOST_RAVE_C6.csv')
    # LR6 = LR6.dropna(subset=['teff_L','logg_L','feh_L','Teff_RAVE','logg_RAVE','[Fe/H]_RAVE'])
    # LR6 = LR6[LR6['logg_RAVE'] > 0]
    C3 = pd.read_csv(ext+'/matlab_in/C3.csv')
    C3 = C3.dropna(subset=['Teff','logg','[Fe/H]'])
    C6 = pd.read_csv(ext+'/matlab_in/C6.csv')
    C6 = C6.dropna(subset=['Teff','logg','[Fe/H]'])

    a, ac = truncate(APO,C6)
    g, gc = truncate(G,C3)
    r3, rc3 = truncate(R3,C3)
    r6, rc6 = truncate(R6,C6)
    # l3, lc3 = truncate(L3,C3)
    l6, lc6 = truncate(L6,C6)
    ar, ra = truncate(APO,R6)
    gr, rg = truncate(G,R3)


    ''' Compare Teff, [Fe/H], logg of given datasets '''
    # comp_scat(LR6,LR6,['teff_L','Teff_RAVE'],['logg_L','logg_RAVE'],['feh_L','[Fe/H]_RAVE'], \
    # r'LAMOST-RAVE: C6','L','L-R',1,'LAMOST_RAVE_C6')

    comp_scat(a,ac,['TEFF','Teff'],['LOGG','logg'],['FE_H','[Fe/H]'], \
    r'APOGEE-EPIC: C6','A','A-Ep',1,'APOGEE_EPIC')

    comp_scat(g,gc,['TEFF','Teff'],['LOGG','logg'],['FEH','[Fe/H]'], \
    r'GES-EPIC: C3','G','G-Ep',1,'GES_EPIC')

    comp_scat(r3,rc3,['Teff_RAVE','Teff'],['logg_RAVE','logg'],['[Fe/H]_RAVE','[Fe/H]'], \
    r'RAVE-EPIC: C3','R','R-Ep',1,'RAVE_EPIC_C3')

    comp_scat(r6,rc6,['Teff_RAVE','Teff'],['logg_RAVE','logg'],['[Fe/H]_RAVE','[Fe/H]'], \
    r'RAVE-EPIC: C6','R','R-Ep',1,'RAVE_EPIC_C6')

    # comp_scat(l3,lc3,['teff_L','Teff'],['logg_L','logg'],['feh_L','[Fe/H]'], \
    # r'LAMOST-EPIC: C3','L','L-Ep',1,'LAMOST_EPIC_C3')

    comp_scat(l6,lc6,['teff_L','Teff'],['logg_L','logg'],['feh_L','[Fe/H]'], \
    r'LAMOST-EPIC: C6','L','L-Ep',1,'LAMOST_EPIC_C6')

    comp_scat(ar,ra,['TEFF','Teff_RAVE'],['LOGG','logg_RAVE'],['FE_H','[Fe/H]_RAVE'], \
    r'APOGEE-RAVE: C6','A','A-R',0,'APOGEE_RAVE')

    comp_scat(gr,rg,['TEFF','Teff_RAVE'],['LOGG','logg_RAVE'],['FEH','[Fe/H]_RAVE'], \
    r'GES-RAVE: C3','G','G-R',0,'GES_RAVE')

    # plt.show()
