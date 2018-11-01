''' Understanding Selection Function Biases '''
''' sys.argv[1] == 0 => data needs reading in and saving '''
''' sys.argv[1] == 1 => data already processed, run program '''
''' sys.argv[1] == 2 => comparison of original field and seismic detections '''
''' sys.argv[1] == 3 => comparison of seismic detections and TRILEGAL '''
''' sys.argv[1] == 4 => calculation of statistical metrics to compare simulated
                        and observed distributions '''

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import sys
import numpy as np
from numbers import Number
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import K2_data as dat
import K2_properties as prop
import Detection_prob_K2 as DP
import K2_constants as const
from numbers import Number
import colormaps
import time

mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'

def seismic_params(data,numax,dnu):
    ''' Calculate the asteroseimic parameters for each star with relevant seismic constants '''

    for i in range(0,len(data),1):
        X = data[i]
        X['slogg'] = np.log10(27400 * (numax[i][i]/X['numax_const']) * np.sqrt(X['Teff']/5777))
        X['mass'] = (numax[i][i]/X['numax_const'])**3 * (dnu[i][i]/X['dnu_const'])**-4 * (X['Teff']/5777)**1.5
        X['radius'] = (numax[i][i]/X['numax_const']) * (dnu[i][i]/X['dnu_const'])**-2 * (X['Teff']/5777)**0.5
        X['L'] = X['radius']**2 * (X['Teff']/5777)**4
        X = det_prob(X,numax[i][i],X['numax_const'][0],X['dnu_const'][0])
        data[i] = X

    return data

def det_prob(X,numax,Numax,Dnu):
    ''' Run detection probability code from Mat, calculating teffred to
    start with and then running the global detection code with added
    noise. Output is the detection probability and SNR

    numax = stellar numax
    Numax = calibration numax
    Dnu = calibration numax '''

    X['Tred'] = DP.properties(X, constants_only=False)
    X['prob_s'], X['SNR'] = DP.globalDetections(X['imag'],X['KepMag'],X['L'],X['radius'],X['Teff'],numax,1.0,X['Tred'], \
                            const.teffred_solar,const.solar_Teff,Numax,Dnu,const.sys_limit,const.dilution,const.vnyq, \
                            const.cadence, vary_beta=True)

    return X

def individ(df,df1,df2):
    ''' Generate outfile containing EPIC numbers that only appear in single
        datasets '''

    Y = pd.DataFrame()
    Y['EPIC'] = df1['EPIC']
    Y['Y'] = 1.0

    B = pd.DataFrame()
    B['EPIC'] = df['EPIC']
    B['B'] = 1.0

    S = pd.DataFrame()
    S['EPIC'] = df2['EPIC']
    S['S'] = 1.0

    Y = Y.join(B[['EPIC','B']].set_index('EPIC'),on='EPIC')
    Y = Y.join(S[['EPIC','S']].set_index('EPIC'),on='EPIC')
    Y['B'].fillna(value=0.0,inplace=True)
    Y['S'].fillna(value=0.0,inplace=True)

    B = B.join(Y[['EPIC','Y']].set_index('EPIC'),on='EPIC')
    B = B.join(S[['EPIC','S']].set_index('EPIC'),on='EPIC')
    B['Y'].fillna(value=0.0,inplace=True)
    B['S'].fillna(value=0.0,inplace=True)

    S = S.join(B[['EPIC','B']].set_index('EPIC'),on='EPIC')
    S = S.join(Y[['EPIC','Y']].set_index('EPIC'),on='EPIC')
    S['B'].fillna(value=0.0,inplace=True)
    S['Y'].fillna(value=0.0,inplace=True)

    Y['Nseismo'] = Y['Y'] + Y['B'] + Y['S']
    B['Nseismo'] = B['Y'] + B['B'] + B['S']
    S['Nseismo'] = S['Y'] + S['B'] + S['S']

    df = pd.merge(df,B,how='inner',on=['EPIC'])
    df = df.reset_index(drop=True)
    df1 = pd.merge(df1,Y,how='inner',on=['EPIC'])
    df1 = df1.reset_index(drop=True)
    df2 = pd.merge(df2,S,how='inner',on=['EPIC'])
    df2 = df2.reset_index(drop=True)

    return df, df1, df2

def sigma_clip(X,a,b,c,d,a1,b1,c1,d1,sigma):
    ''' Reduce data set to those values that fall within an n sigma difference
        of each other for numax and dnu '''
    x = len(X['EPIC'])
    Y = X
    Y['comp_err_Nmx'] = np.sqrt(Y[a1]**2+Y[b1]**2)
    Y = Y[Y['comp_err_Nmx'] >= -4]
    Y = Y[Y['comp_err_Nmx'] <= 4]
    # print( max(Y['comp_err_Nmx']))
    Y['Nmx_Diff'] = ((Y[a]-Y[b])/Y['comp_err_Nmx'])
    # print( Y['Nmx_Diff'])
    Y = Y[Y['Nmx_Diff'] >= -sigma]
    Y = Y[Y['Nmx_Diff'] <= sigma]


    Y['comp_err_Dnu'] = np.sqrt(Y[c1]**2+Y[d1]**2)
    Y = Y[Y['comp_err_Dnu'] >= -0.9]
    Y = Y[Y['comp_err_Dnu'] <= 0.9]
    # print( max(Y['comp_err_Dnu']))
    Y['Dnu_Diff'] = ((Y[c]-Y[d])/Y['comp_err_Dnu'])
    # print( Y['Dnu_Diff'])
    Y = Y[Y['Dnu_Diff'] >= -sigma]
    Y = Y[Y['Dnu_Diff'] <= sigma]
    # print( len(Y), len(X))
    fract_outl = (float(len(X))-float(len(Y)))/float(x)
    # print( fract_outl)
    ''' Z stores the values that were rejected '''
    Z = pd.concat([X,Y],ignore_index=True)
    Z = Z.drop_duplicates(subset=['EPIC'],keep=False)
    Z = Z.reset_index(drop=True)
    X = Y

    return X

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

    hist0[0], d[0], c = ax0.hist(df['mass'],bins=bins[0],histtype='step',linewidth=2)
    hist[0], b[0], c = ax0.hist(df1['mass'],bins=bins[0],histtype='step',linewidth=2)
    ax0.set_xlabel(r'Mass [M$_{\odot}$]')
    if n == 6:
        ax0.set_xlim(0.5,2.5)

    hist0[1], d[1], c = ax1.hist(df['Radius'],bins=bins[1],histtype='step',linewidth=2)
    hist[1], b[1], c = ax1.hist(df1['Radius'],bins=bins[1],histtype='step',linewidth=2)
    ax1.set_xlabel(r'Radius [R$_{\odot}$]')
    if n == 6:
        ax1.set_xlim(0,30)

    hist0[2], d[2], c = ax2.hist(df['Teff'],bins=bins[2],histtype='step',linewidth=2)
    hist[2], b[2], c = ax2.hist(df1['Teff'],bins=bins[2],histtype='step',linewidth=2)
    ax2.set_xlabel(r'T$_{\rm{eff}}$ [K]')
    if n == 6:
        ax2.set_xlim(4000,5250)

    hist0[3], d[3], c = ax3.hist(df['[Fe/H]'],bins=bins[3],histtype='step',linewidth=2)
    hist[3], b[3], c = ax3.hist(df1['[Fe/H]'],bins=bins[3],histtype='step',linewidth=2)
    ax3.set_xlabel(r'[Fe/H] [dex]')
    if n == 6:
        ax3.set_xlim(-1.2,0.2)

    hist0[4], d[4], c = ax4.hist(df['logg'],bins=bins[4],histtype='step',linewidth=2)
    hist[4], b[4], c = ax4.hist(df1['logg'],bins=bins[4],histtype='step',linewidth=2)
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
    # plt.savefig('/home/bmr135/sel_func_comp/'+str(n)+'_'+time.strftime("%d%m%Y")+'_C6.png')

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
        # plt.savefig('/home/bmr135/sel_func_comp/Percentage_reductions_'+time.strftime("%d%m%Y")+'_C6.png')
        # plt.show()


    return hist, b

def percentage_decrease(h0,h1,h2,h3,h4,h5,h6,b,field):
    '''
    - Calculate percentage decrease of stars in a bin for each parameter.
    '''
    mass = pd.DataFrame()
    mass['numax'] = np.nan_to_num(100*(h0[0]-h1[0])/h0[0])
    mass['mag'] = np.nan_to_num(100*(h1[0]-h2[0])/h1[0])
    mass['JK'] = np.nan_to_num(100*(h2[0]-h3[0])/h2[0])
    mass['detP'] = np.nan_to_num(100*(h3[0]-h4[0])/h3[0])
    mass['Nseis'] = np.nan_to_num(100*(h4[0]-h5[0])/h4[0])
    mass['Sig'] = np.nan_to_num(100*(h5[0]-h6[0])/h5[0])
    mass['OVERALL'] = np.nan_to_num(100*(h0[0]-h6[0])/h0[0])
    mass = mass.rename(lambda x: b[0][x])

    Radius = pd.DataFrame()
    Radius['numax'] = np.nan_to_num(100*(h0[1]-h1[1])/h0[1])
    Radius['mag'] = np.nan_to_num(100*(h1[1]-h2[1])/h1[1])
    Radius['JK'] = np.nan_to_num(100*(h2[1]-h3[1])/h2[1])
    Radius['detP'] = np.nan_to_num(100*(h3[1]-h4[1])/h3[1])
    Radius['Nseis'] = np.nan_to_num(100*(h4[1]-h5[1])/h4[1])
    Radius['Sig'] = np.nan_to_num(100*(h5[1]-h6[1])/h5[1])
    Radius['OVERALL'] = np.nan_to_num(100*(h0[1]-h6[1])/h0[1])
    Radius = Radius.rename(lambda x: b[1][x])

    Teff = pd.DataFrame()
    Teff['numax'] = np.nan_to_num(100*(h0[2]-h1[2])/h0[2])
    Teff['mag'] = np.nan_to_num(100*(h1[2]-h2[2])/h1[2])
    Teff['JK'] = np.nan_to_num(100*(h2[2]-h3[2])/h2[2])
    Teff['detP'] = np.nan_to_num(100*(h3[2]-h4[2])/h3[2])
    Teff['Nseis'] = np.nan_to_num(100*(h4[2]-h5[2])/h4[2])
    Teff['Sig'] = np.nan_to_num(100*(h5[2]-h6[2])/h5[2])
    Teff['OVERALL'] = np.nan_to_num(100*(h0[2]-h6[2])/h0[2])
    Teff = Teff.rename(lambda x: b[2][x])

    feh = pd.DataFrame()
    feh['numax'] = np.nan_to_num(100*(h0[3]-h1[3])/h0[3])
    feh['mag'] = np.nan_to_num(100*(h1[3]-h2[3])/h1[3])
    feh['JK'] = np.nan_to_num(100*(h2[3]-h3[3])/h2[3])
    feh['detP'] = np.nan_to_num(100*(h3[3]-h4[3])/h3[3])
    feh['Nseis'] = np.nan_to_num(100*(h4[3]-h5[3])/h4[3])
    feh['Sig'] = np.nan_to_num(100*(h5[3]-h6[3])/h5[3])
    feh['OVERALL'] = np.nan_to_num(100*(h0[3]-h6[3])/h0[3])
    feh = feh.rename(lambda x: b[3][x])

    logg = pd.DataFrame()
    logg['numax'] = np.nan_to_num(100*(h0[4]-h1[4])/h0[4])
    logg['mag'] = np.nan_to_num(100*(h1[4]-h2[4])/h1[4])
    logg['JK'] = np.nan_to_num(100*(h2[4]-h3[4])/h2[4])
    logg['detP'] = np.nan_to_num(100*(h3[4]-h4[4])/h3[4])
    logg['Nseis'] = np.nan_to_num(100*(h4[4]-h5[4])/h4[4])
    logg['Sig'] = np.nan_to_num(100*(h5[4]-h6[0])/h5[0])
    logg['OVERALL'] = np.nan_to_num(100*(h0[4]-h6[4])/h0[4])
    logg = logg.rename(lambda x: b[4][x])

    if field == 3:
        mass.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3/mass_percentage_reduction'+time.strftime("%d%m%Y"))
        Radius.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3/Radius_percentage_reduction'+time.strftime("%d%m%Y"))
        Teff.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3/Teff_percentage_reduction'+time.strftime("%d%m%Y"))
        feh.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3/feh_percentage_reduction'+time.strftime("%d%m%Y"))
        logg.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3/logg_percentage_reduction'+time.strftime("%d%m%Y"))
        np.savetxt('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3/bins'+time.strftime("%d%m%Y"),b,delimiter=',')
    if field == 6:
        mass.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6/mass_percentage_reduction'+time.strftime("%d%m%Y"))
        Radius.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6/Radius_percentage_reduction'+time.strftime("%d%m%Y"))
        Teff.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6/Teff_percentage_reduction'+time.strftime("%d%m%Y"))
        feh.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6/feh_percentage_reduction'+time.strftime("%d%m%Y"))
        logg.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6/logg_percentage_reduction'+time.strftime("%d%m%Y"))
        np.savetxt('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6/bins'+time.strftime("%d%m%Y"),b,delimiter=',')

def extents(f):
  delta = f[1] - f[0]
  return [f[0] - delta/2, f[-1] + delta/2]

def likeli(X,mu,sig,param,prob):
    ''' Gaussian likelihood '''
    X['diff'] = (X[param] - mu)**2
    # X[prob] = -(len(X)/2)*np.log(2*np.pi*sig) - (1/(2*sig))*X['diff'] # varying sigma
    X[prob] = - (1/(2*sig))*X['diff'] # fixed sigma
    # X[prob] = (1/(np.sqrt(2*np.pi*sig**2))) * np.exp(-0.5*((X['diff'] - mu)/sig)**2)
    X[prob] = np.exp(X[prob])
    return X

def likeli_plt(sim,dat,param1,param2,prob):
    mean, bin_edges, binnumber = scipy.stats.binned_statistic(sim[param1],sim[prob],'mean',bins=50)
    plt.figure()
    plt.hist(sim[param1], bins=50, normed=True, histtype='stepfilled',alpha=0.2, label='histogram of data')
    plt.hist(dat[param2], bins=50, normed=True, histtype='stepfilled',alpha=0.2, label='histogram of data')
    plt.hlines(mean, bin_edges[:-1], bin_edges[1:], colors='g', lw=2, label='binned statistic of data')
    plt.legend(fontsize=10)

def single_seismo(df,keys,output):
    ''' Generation of single column of seismic values '''
    a,b,c = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
    a['EPIC'],b['EPIC'],c['EPIC'] = df['EPIC'],df['EPIC'],df['EPIC']
    a[output] = df[keys[0]]
    a = a.dropna(subset=[output])
    b[output] = df[keys[1]].dropna(axis=0,how='any')
    b = b.dropna(subset=[output])
    c[output] = df[keys[2]].dropna(axis=0,how='any')
    c = c.dropna(subset=[output])
    df1 = pd.concat([a,b,c],ignore_index=True)
    df1 = df1.drop_duplicates(subset=['EPIC'])
    df1 = df1.reset_index(drop=True)
    if len(df) == len(df1):
        df = pd.merge(df,df1,how='inner',on=['EPIC'])
        df = df.reset_index(drop=True)
    else:
        df = df.join(df1.set_index('EPIC'),on='EPIC')
    # print(df[output])
    return df

if __name__ == '__main__':

    ''' Read in necessary data files '''
    if sys.argv[1] == '0':

        GAP3, GAP6 = dat.K2_GAP()
        Yvonne_C3, Yvonne_C6, Yvonne_EC6, Yvonne_EC3 = dat.Yvonne()
        Savita_C3, Savita_C6, Savita_EC3, Savita_EC6 = dat.Savita()
        Benoit_C3, Benoit_C6, Everest_C3, Everest_C6 = dat.Benoit()
        # GAP3 = GAP3.drop(columns=['Radius'])
        # GAP3.rename(columns={'radius_val':'Radius'},inplace=True)
        # GAP6 = GAP6.drop(columns=['Radius'])
        # GAP6.rename(columns={'radius_val':'Radius'},inplace=True)

        C3 = [Benoit_C3,Yvonne_C3,Savita_C3]
        YC3, BC3, SC3 = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        C3m = [BC3,YC3,SC3] # asteroseismic values merged with GAP
        C6 = [Benoit_C6,Yvonne_C6,Savita_C6]
        YC6, BC6, SC6 = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        C6m = [BC6,YC6,SC6] # asteroseismic values merged with GAP


        C3[0],C3[1],C3[2] = individ(C3[0],C3[1],C3[2])
        C6[0],C6[1],C6[2] = individ(C6[0],C6[1],C6[2])

        for i in range(3):
            C3m[i] = pd.merge(GAP3,C3[i],how='inner',on=['EPIC'])
            C3m[i] = C3m[i].reset_index(drop=True)

        for i in range(3):
            C6m[i] = pd.merge(GAP6,C6[i],how='inner',on=['EPIC'])
            C6m[i] = C6m[i].reset_index(drop=True)

        C3m[0]['numax_const'], C6m[0]['numax_const'] = 3104, 3104
        C3m[0]['dnu_const'], C6m[0]['dnu_const'] = 138.8, 138.8
        C3m[1]['numax_const'], C6m[1]['numax_const'] = 3135, 3135
        C3m[1]['dnu_const'], C6m[1]['dnu_const'] = 135.045, 135.045
        C3m[2]['numax_const'], C6m[2]['numax_const'] = 3097.33, 3097.33
        C3m[2]['dnu_const'], C6m[2]['dnu_const'] = 135.2, 135.2

        nmxB3, nmxB6, nmxY3, nmxY6, nmxS3, nmxS6 = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        dnuB3, dnuB6, dnuY3, dnuY6, dnuS3, dnuS6 = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        nmxB3[0], nmxY3[1], nmxS3[2] = C3m[0]['Bnumax'], C3m[1]['nmx'], C3m[2]['Snumax']
        nmxB6[0], nmxY6[1], nmxS6[2] = C6m[0]['Bnumax'], C6m[1]['nmx'], C6m[2]['Snumax']
        dnuB3[0], dnuY3[1], dnuS3[2] = C3m[0]['BDnu'], C3m[1]['dnu'], C3m[2]['SDnu']
        dnuB6[0], dnuY6[1], dnuS6[2] = C6m[0]['BDnu'], C6m[1]['dnu'], C6m[2]['SDnu']

        numax_C3 = [nmxB3,nmxY3,nmxS3]
        numax_C6 = [nmxB6,nmxY6,nmxS6]
        dnu_C3 = [dnuB3,dnuY3,dnuS3]
        dnu_C6 = [dnuB6,dnuY6,dnuS6]

        C3m = seismic_params(C3m,numax_C3,dnu_C3)
        C6m = seismic_params(C6m,numax_C6,dnu_C6)

        ''' Single Numax/Dnu creation for C3 '''
        K2_3 = pd.concat([C3m[0],C3m[1],C3m[2]],ignore_index=True)
        K2_3 = K2_3.drop_duplicates(subset=['EPIC'])
        K2_3 = K2_3.reset_index(drop=True)
        # K2_3['Glogg'] = np.log10((6.6742e-8*K2_3['mass']*1.9884e33)/(K2_3['Radius']*6.957e10)**2)

        K2_3 = single_seismo(K2_3,['Bnumax','nmx','Snumax'],'NUMAX')
        K2_3 = single_seismo(K2_3,['BDnu','dnu','SDnu'],'DNU')
        K2_3 = single_seismo(K2_3,['e_Bnumax','nmx_err','e_Snumax'],'NUMAX_err')
        K2_3 = single_seismo(K2_3,['e_BDnu','dnu_err','e_SDnu'],'DNU_err')


        ''' Single Numax/Dnu creation for C6 '''
        K2_6 = pd.concat([C6m[0],C6m[1],C6m[2]],ignore_index=True)
        K2_6 = K2_6.drop_duplicates(subset=['EPIC'])
        K2_6 = K2_6.reset_index(drop=True)
        # K2_6['Glogg'] = np.log10((6.6742e-8*K2_6['mass']*1.9884e33)/(K2_6['Radius']*6.957e10)**2)

        K2_6 = single_seismo(K2_6,['Bnumax','nmx','Snumax'],'NUMAX')
        K2_6 = single_seismo(K2_6,['BDnu','dnu','SDnu'],'DNU')
        K2_6 = single_seismo(K2_6,['e_Bnumax','nmx_err','e_Snumax'],'NUMAX_err')
        K2_6 = single_seismo(K2_6,['e_BDnu','dnu_err','e_SDnu'],'DNU_err')


        ''' Save full list of stars here for all tests bar Sigma Clip '''
        K2_3.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3_full_'+time.strftime("%d%m%Y"),index=False,na_rep='NaN')
        K2_6.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6_full_'+time.strftime("%d%m%Y"),index=False,na_rep='NaN')

        Camp3 = pd.DataFrame()
        Camp3['EPIC'] = K2_3['EPIC']
        Camp3['NUMAX'] = K2_3['NUMAX']
        Camp3['DNU'] = K2_3['DNU']
        Camp3['NUMAX_err'] = K2_3['NUMAX_err']
        Camp3['DNU_err'] = K2_3['DNU_err']
        Camp6 = pd.DataFrame()
        Camp6['EPIC'] = K2_6['EPIC']
        Camp6['NUMAX'] = K2_6['NUMAX']
        Camp6['DNU'] = K2_6['DNU']
        Camp6['NUMAX_err'] = K2_6['NUMAX_err']
        Camp6['DNU_err'] = K2_6['DNU_err']

        BY3,YS3,BS3,BY6,YS6,BS6 = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        C3m2 = [BY3,YS3,BS3] # merging of asteroseismic pipelines
        C6m2 = [BY6,YS6,BS6] # merging of asteroseismic pipelines

        ''' C3 detection mergers '''
        cols_to_use = C3m[1].columns.difference(C3m[0].columns)
        cols_to_use = cols_to_use.union(['EPIC'])
        C3m2[0] = pd.merge(C3m[0],C3m[1][cols_to_use],how='inner',on=['EPIC'])
        cols_to_use = C3m[2].columns.difference(C3m[1].columns)
        cols_to_use = cols_to_use.union(['EPIC'])
        C3m2[1] = pd.merge(C3m[1],C3m[2][cols_to_use],how='inner',on=['EPIC'])
        cols_to_use = C3m[2].columns.difference(C3m[0].columns)
        cols_to_use = cols_to_use.union(['EPIC'])
        C3m2[2] = pd.merge(C3m[0],C3m[2][cols_to_use],how='inner',on=['EPIC'])

        ''' C6 detection mergers '''
        cols_to_use = C6m[1].columns.difference(C6m[0].columns)
        cols_to_use = cols_to_use.union(['EPIC'])
        C6m2[0] = pd.merge(C6m[0],C6m[1][cols_to_use],how='inner',on=['EPIC'])
        cols_to_use = C6m[2].columns.difference(C6m[1].columns)
        cols_to_use = cols_to_use.union(['EPIC'])
        C6m2[1] = pd.merge(C6m[1],C6m[2][cols_to_use],how='inner',on=['EPIC'])
        cols_to_use = C6m[2].columns.difference(C6m[0].columns)
        cols_to_use = cols_to_use.union(['EPIC'])
        C6m2[2] = pd.merge(C6m[0],C6m[2][cols_to_use],how='inner',on=['EPIC'])

        ''' Sigma Clipping '''

        C3m2[0] = sigma_clip(C3m2[0],'Bnumax','nmx','BDnu','dnu', \
                    'e_Bnumax','nmx_err','e_BDnu','dnu_err',3)
        C3m2[1] = sigma_clip(C3m2[1],'nmx','Snumax','dnu','SDnu', \
            'nmx_err','e_Snumax','dnu_err','e_SDnu',3)
        C3m2[2] = sigma_clip(C3m2[2],'Bnumax','Snumax','BDnu','SDnu', \
                    'e_Bnumax','e_Snumax','e_BDnu','e_SDnu',3)

        C6m2[0] = sigma_clip(C6m2[0],'Bnumax','nmx','BDnu','dnu', \
                    'e_Bnumax','nmx_err','e_BDnu','dnu_err',3)
        C6m2[1] = sigma_clip(C6m2[1],'nmx','Snumax','dnu','SDnu', \
                    'nmx_err','e_Snumax','dnu_err','e_SDnu',3)
        C6m2[2] = sigma_clip(C6m2[2],'Bnumax','Snumax','BDnu','SDnu', \
                    'e_Bnumax','e_Snumax','e_BDnu','e_SDnu',3)

        ''' SIGMA CLIP NUMAX/DNU MERGER '''

        K2_3m = pd.concat([C3m2[0],C3m2[1],C3m2[2]],ignore_index=True)
        K2_3m = K2_3m.drop_duplicates(subset=['EPIC'])
        K2_3m = K2_3m.reset_index(drop=True)
        K2_3m = pd.merge(K2_3m,Camp3,how='inner',on=['EPIC'])
        K2_3m = K2_3m.reset_index(drop=True)

        K2_6m = pd.concat([C6m2[0],C6m2[1],C6m2[2]],ignore_index=True)
        K2_6m = K2_6m.drop_duplicates(subset=['EPIC'])
        K2_6m = K2_6m.reset_index(drop=True)
        K2_6m = pd.merge(K2_6m,Camp6,how='inner',on=['EPIC'])
        K2_6m = K2_6m.reset_index(drop=True)

        ''' Save Sigma Clipped Stars '''
        K2_3m.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3_full_sig_clip_'+time.strftime("%d%m%Y"),index=False,na_rep='NaN')
        K2_6m.to_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6_full_sig_clip_'+time.strftime("%d%m%Y"),index=False,na_rep='NaN')

    ''' Step-by-Step Breakdown of Selection Function '''
    if sys.argv[1] == '1':

        C3orig = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3_full_01112018')
        C6orig = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6_full_01112018')
        SigC3 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3_full_sig_clip_01112018')
        SigC6 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6_full_sig_clip_01112018')

        ''' Selection Function Pre-Applied to Sigma-Clipped Stars '''
        SigC3 = SigC3[ (SigC3['NUMAX'] > 10) & (SigC3['NUMAX'] < 280) & (SigC3['Hmag'] > 7) & (SigC3['Hmag'] < 12) & (SigC3['JK'] > 0.5) & (SigC3['prob_s'] >= 0.95) ]
        SigC6 = SigC6[ (SigC6['NUMAX'] > 10) & (SigC6['NUMAX'] < 280) & (SigC6['Vcut'] > 9) & (SigC6['Vcut'] < 15) & (SigC6['JK'] > 0.5) & (SigC6['prob_s'] >= 0.95) ]

        C3orig = C3orig[C3orig['mass']>0.0]
        C6orig = C6orig[C6orig['mass']>0.0]
        C6orig = C6orig.dropna(subset=['Radius'])
        param = ['mass','radius','Teff','[Fe/H]','slogg']
        fancy = [r'Mass [M$_{\odot}$]',r'Radius [R$_{\odot}$]',r'T$_{\rm{eff}}$',r'[Fe/H]',r'log$_{10}$(g)']
        ext3 = '/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3/'
        ext6 = '/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6/'
        ext = [ext3,ext6]
        inputs = [C3orig,C6orig]
        field = [3,6]
        histo1 = [0,0,0,0,0]

        smol_rad3 = C3orig[C3orig['Radius'] < 2]
        print(smol_rad3[['EPIC','Radius','logg','NUMAX']])
        smol_rad6 = C6orig[C6orig['Radius'] < 2]
        print(smol_rad6[['EPIC','Radius','logg','NUMAX']])
        sys.exit()

        j=0
        for i in inputs:

            ''' Bins '''
            if field[j] == 3:
                bins = [np.linspace(min(C3orig['mass']),max(C3orig['mass']),50), \
                        np.linspace(min(C3orig['Radius']),max(C3orig['Radius']),50), \
                        np.linspace(min(C3orig['Teff']),max(C3orig['Teff']),50), \
                        np.linspace(min(C3orig['[Fe/H]']),max(C3orig['[Fe/H]']),50), \
                        np.linspace(min(C3orig['logg']),max(C3orig['logg']),50)]
                hist0, b0 = hist_orig(C3orig,C3orig,r'None',bins,ext[0],0)
            if field[j] == 6:
                bins = [np.linspace(min(C6orig['mass']),max(C6orig['mass']),50), \
                        np.linspace(min(C6orig['Radius']),25,50), \
                        np.linspace(min(C6orig['Teff']),max(C6orig['Teff']),50), \
                        np.linspace(min(C6orig['[Fe/H]']),max(C6orig['[Fe/H]']),50), \
                        np.linspace(min(C6orig['logg']),max(C6orig['logg']),50)]
                hist0, b0 = hist_orig(C6orig,C6orig,r'None',bins,ext[1],0)
            '''
            NuMax cuts
            cut[0] => reduced data
            cut[1] => removed data set
            '''
            nm1 = nm2 = pd.DataFrame()
            cut = [0, 0]
            cut[0] = i[i['NUMAX'] <= 280]
            cut[0] = cut[0][cut[0]['NUMAX'] >= 10]
            nm1 = i[i['NUMAX'] > 280]
            nm2 = i[i['NUMAX'] < 10]
            cut[1] = pd.concat([nm1,nm2],ignore_index=True)
            cut[1] = cut[1].drop_duplicates(subset=['EPIC'],keep=False)
            cut[1] = cut[1].reset_index(drop=True)
            if field[j] == 3:
                hist1, b = hist_orig(i,cut[0],r'$\nu_{\rm{max}}$ (\#1)',bins,ext[0],1)
            if field[j] == 6:
                hist1, b = hist_orig(i,cut[0],r'$\nu_{\rm{max}}$ (\#1)',bins,ext[1],1)
            numax = cut[1]
            cut[1] = 0.0
            m=cut[0]
            # fig, axes = plt.subplots(3,2)
            # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
            # ax = [ax0,ax1,ax2,ax3,ax4,ax5]
            # for k in range(len(param)):
            #     ax[k].hist([numax[param[k]]],bins=bins[k],stacked=True,label=[r'$\nu_{\rm{max}}$'])
            #     ax[k].set_xlabel(fancy[k])
            # box = ax[4].get_position()
            # ax[4].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
            # ax[4].legend(loc='upper center', bbox_to_anchor=(1.6, 1.0),fancybox=True, shadow=True, ncol=1)
            # ax[5].axis('off')
            # plt.tight_layout()
            # if field[j] == 3:
            #     # plt.savefig(ext3+'numax.png')
            # if field[j] == 6:
            #     # plt.savefig(ext6+'numax.png')

            ''' Magnitude Cuts '''
            if field[j] == 3:
                cut[1] = cut[0][cut[0]['Hmag'] > 12]
                cut[1] = cut[0][cut[0]['Hmag'] < 7]
                cut[0] = cut[0][cut[0]['Hmag'] <= 12]
                cut[0] = cut[0][cut[0]['Hmag'] >= 7]
                hist2, b = hist_orig(i,cut[0],r'H-band (\#2)',bins,ext[0],2)
                mag = cut[1]
                m1=cut[0]
                cut[1] = 0.0
                # fig, axes = plt.subplots(3,2)
                # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
                # ax = [ax0,ax1,ax2,ax3,ax4,ax5]
                # for k in range(len(param)):
                #     ax[k].hist([numax[param[k]],mag[param[k]]],bins=bins[k],stacked=True,label=[r'$\nu_{\rm{max}}$',r'H-Band'])
                #     ax[k].set_xlabel(fancy[k])
                # box = ax[4].get_position()
                # ax[4].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
                # ax[4].legend(loc='upper center', bbox_to_anchor=(1.6, 1.0),fancybox=True, shadow=True, ncol=1)
                # ax[5].axis('off')
                # plt.tight_layout()
                # plt.savefig(ext3+'numax_mag.png')

            if field[j] == 6:
                cut[1] = cut[0][cut[0]['Vcut'] > 15]
                cut[1] = cut[0][cut[0]['Vcut'] < 9]
                cut[0] = cut[0][cut[0]['Vcut'] < 15]
                cut[0] = cut[0][cut[0]['Vcut'] > 9]
                hist2, b = hist_orig(i,cut[0],r'H-band (\#2)',bins,ext[1],2)
                mag = cut[1]
                m1=cut[0]
                cut[1] = 0.0
                # fig, axes = plt.subplots(3,2)
                # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
                # ax = [ax0,ax1,ax2,ax3,ax4,ax5]
                # for k in range(len(param)):
                #     ax[k].hist([numax[param[k]],mag[param[k]]],bins=bins[k],stacked=True,label=[r'$\nu_{\rm{max}}$',r'V-Band'])
                #     ax[k].set_xlabel(fancy[k])
                # box = ax[4].get_position()
                # ax[4].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
                # ax[4].legend(loc='upper center', bbox_to_anchor=(1.6, 1.0),fancybox=True, shadow=True, ncol=1)
                # ax[5].axis('off')
                # plt.tight_layout()
                # plt.savefig(ext6+'numax_mag.png')

            ''' J-K '''
            cut[1] = cut[0][cut[0]['JK'] < 0.5]
            cut[0] = cut[0][cut[0]['JK'] >= 0.5]
            if field[j] == 3:
                hist3, b = hist_orig(i,cut[0],r'J-K (\#3)',bins,ext[0],3)
            if field[j] == 6:
                hist3, b = hist_orig(i,cut[0],r'J-K (\#3)',bins,ext[1],3)
            JK = cut[1]
            m2=cut[0]
            cut[1] = 0.0
            # fig, axes = plt.subplots(3,2)
            # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
            # ax = [ax0,ax1,ax2,ax3,ax4,ax5]
            # for k in range(len(param)):
            #     ax[k].hist([numax[param[k]],mag[param[k]],JK[param[k]]],bins=bins[k],stacked=True,label=[r'$\nu_{\rm{max}}$',r'H-Band',r'J-K'])
            #     ax[k].set_xlabel(fancy[k])
            # box = ax[4].get_position()
            # ax[4].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
            # ax[4].legend(loc='upper center', bbox_to_anchor=(1.6, 1.0),fancybox=True, shadow=True, ncol=1)
            # ax[5].axis('off')
            # plt.tight_layout()
            # if field[j] == 3:
            #     # plt.savefig(ext3+'numax_mag_JK.png')
            # if field[j] == 6:
            #     # plt.savefig(ext6+'numax_mag_JK.png')


            ''' Detection Probability '''
            cut[1] = cut[0][cut[0]['prob_s'] < 0.95]
            cut[0] = cut[0][cut[0]['prob_s'] >= 0.95]
            if field[j] == 3:
                hist4, b = hist_orig(i,cut[0],r'Det. Prob. (\#4)',bins,ext[0],4)
            if field[j] == 6:
                hist4, b = hist_orig(i,cut[0],r'Det. Prob. (\#4)',bins,ext[1],4)
            detPro = cut[1]
            m3=cut[0]
            cut[1] = 0.0
            # fig, axes = plt.subplots(3,2)
            # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
            # ax = [ax0,ax1,ax2,ax3,ax4,ax5]
            # for k in range(len(param)):
            #     ax[k].hist([numax[param[k]],mag[param[k]],JK[param[k]],detPro[param[k]]],bins=bins[k],stacked=True,label=[r'$\nu_{\rm{max}}$',r'H-Band',r'J-K',r'Det. Prob.'])
            #     ax[k].set_xlabel(fancy[k])
            # box = ax[4].get_position()
            # ax[4].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
            # ax[4].legend(loc='upper center', bbox_to_anchor=(1.6, 1.0),fancybox=True, shadow=True, ncol=2)
            # ax[5].axis('off')
            # plt.tight_layout()
            # if field[j] == 3:
            #     # plt.savefig(ext3+'numax_mag_JK_detP.png')
            # if field[j] == 6:
            #     # plt.savefig(ext6+'numax_mag_JK_detP.png')

            ''' Mulitple Detections '''
            cut[1] = cut[0][cut[0]['Nseismo'] < 2]
            cut[0] = cut[0][cut[0]['Nseismo'] >= 2]
            if field[j] == 3:
                hist5, b = hist_orig(i,cut[0],r'Number Seismic Dets. (\#5)',bins,ext[0],5)
            if field[j] == 6:
                hist5, b = hist_orig(i,cut[0],r'Number Seismic Dets. (\#5)',bins,ext[1],5)
            Nseis = cut[1]
            m4=cut[0]
            cut[1] = 0.0
            # fig, axes = plt.subplots(3,2)
            # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
            # ax = [ax0,ax1,ax2,ax3,ax4,ax5]
            # for k in range(len(param)):
            #     ax[k].hist([numax[param[k]],mag[param[k]],JK[param[k]],detPro[param[k]],Nseis[param[k]]],bins=bins[k],stacked=True,label=[r'$\nu_{\rm{max}}$',r'H-Band',r'J-K',r'Det. Prob.',r'N$_{seismic}$'])
            #     ax[k].set_xlabel(fancy[k])
            # box = ax[4].get_position()
            # ax[4].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
            # ax[4].legend(loc='upper center', bbox_to_anchor=(1.6, 1.0),fancybox=True, shadow=True, ncol=2)
            # ax[5].axis('off')
            # plt.tight_layout()
            # if field[j] == 3:
            #     # plt.savefig(ext3+'numax_mag_JK_detP_Nseis.png')
            # if field[j] == 6:
            #     # plt.savefig(ext6+'numax_mag_JK_detP_Nseis.png')

            ''' Sigma Clip '''
            if field[j] ==3:
                cutt = pd.merge(cut[0],SigC3[['EPIC']],how='inner',on=['EPIC'])
                cutt = cutt.reset_index(drop=True)
                cut[1] = pd.concat([cut[0],SigC3],ignore_index=True)
                cut[1] = cut[1].drop_duplicates(subset=['EPIC'],keep=False)
                cut[1] = cut[1].reset_index(drop=True)
                hist6, b = hist_orig(i,cutt,r'Sigma Clip. (\#6)',bins,ext[0],6)
                m5=cut[0]
                SC = cut[1]
                # fig, axes = plt.subplots(3,2)
                # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
                # ax = [ax0,ax1,ax2,ax3,ax4,ax5]
                # for k in range(len(param)):
                #     ax[k].hist([numax[param[k]],mag[param[k]],JK[param[k]],detPro[param[k]],Nseis[param[k]],SC[param[k]]],bins=bins[k],stacked=True,label=[r'$\nu_{\rm{max}}$',r'H-Band',r'J-K',r'Det. Prob.',r'N$_{seismic}$',r'Sig. Clip'])
                #     ax[k].set_xlabel(fancy[k])
                # box = ax[4].get_position()
                # ax[4].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
                # ax[4].legend(loc='upper center', bbox_to_anchor=(1.6, 1.0),fancybox=True, shadow=True, ncol=2)
                # ax[5].axis('off')
                # plt.tight_layout()
                # plt.savefig(ext3+'numax_mag_JK_detP_Nseis_Sig.png')

            if field[j] ==6:
                cutt = pd.merge(cut[0],SigC6[['EPIC']],how='inner',on=['EPIC'])
                cutt = cutt.reset_index(drop=True)
                cut[1] = pd.concat([cut[0],SigC6],ignore_index=True)
                cut[1] = cut[1].drop_duplicates(subset=['EPIC'],keep=False)
                cut[1] = cut[1].reset_index(drop=True)
                hist6, b = hist_orig(i,cutt,r'Sigma Clip. (\#6)',bins,ext[1],6)
                m5=cut[0]
                SC = cut[1]
                # fig, axes = plt.subplots(3,2)
                # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
                # ax = [ax0,ax1,ax2,ax3,ax4,ax5]
                # for k in range(len(param)):
                #     ax[k].hist([numax[param[k]],mag[param[k]],JK[param[k]],detPro[param[k]],Nseis[param[k]],SC[param[k]]],bins=50,stacked=True,label=[r'$\nu_{\rm{max}}$',r'H-Band',r'J-K',r'Det. Prob.',r'N$_{seismic}$',r'Sig. Clip'])
                #     ax[k].set_xlabel(fancy[k])
                # box = ax[4].get_position()
                # ax[4].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
                # ax[4].legend(loc='upper center', bbox_to_anchor=(1.6, 1.0),fancybox=True, shadow=True, ncol=2)
                # ax[5].axis('off')
                # plt.tight_layout()
                # plt.savefig(ext6+'numax_mag_JK_detP_Nseis_Sig.png')

            percentage_decrease(hist0,hist1,hist2,hist3,hist4,hist5,hist6,b0,field[j])
            j+=1

            fig, ax = plt.subplots()
            ax.scatter(i['Teff'],i['slogg'],label=r'Full')
            ax.scatter(m['Teff'],m['slogg'],label=r'$\nu_{\rm{max}}$')
            ax.scatter(m1['Teff'],m1['slogg'],label=r'Mag.')
            ax.scatter(m2['Teff'],m2['slogg'],label=r'Colour')
            ax.scatter(m3['Teff'],m3['slogg'],label=r'Det. Prob.')
            ax.scatter(m4['Teff'],m4['slogg'],label=r'N$_{\rm{det}}$')
            ax.scatter(m5['Teff'],m5['slogg'],label=r'$\sigma_{\rm{clip}}$')
            ax.set_xlabel(r'$T_{\rm{eff}}$')
            ax.set_ylabel(r'log$_{10}$(g)')
            ax.legend()
            plt.gca().invert_xaxis()
            plt.gca().invert_yaxis()
            # plt.savefig('/home/bmr135/sel_func_comp/Reduction_HRD_'+time.strftime("%d%m%Y")+'_C6.png')
            # plt.show()
            # sys.exit()




        # plt.show()

    ''' Direct Comparison of Original Fields and Seismic Data '''
    if sys.argv[1] == '2':
        C3orig = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3_full')
        C6orig = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6_full')
        GAP3, GAP6 = dat.K2_GAP()

        C3orig = C3orig[ (C3orig['Hmag'] > 7) & (C3orig['Hmag'] < 12) & (C3orig['JK'] > 0.5) & (C3orig['mass'] > 0)]
        C6orig = C6orig[ (C6orig['Vcut'] > 9) & (C6orig['Vcut'] < 15) & (C6orig['JK'] > 0.5) & (C6orig['mass'] > 0)]
        GAP3 = GAP3[ (GAP3['Hmag'] > 7) & (GAP3['Hmag'] < 12) & (GAP3['JK'] > 0.5) & (GAP3['mass'] > 0)]
        GAP6 = GAP6[ (GAP6['Vcut'] > 9) & (GAP6['Vcut'] < 15) & (GAP6['JK'] > 0.5) & (GAP6['mass'] > 0)]

        ext_fig = '/home/bmr135/Dropbox/K2Poles/pop_trends/231018/GAP_vs_Dets/'

        ''' Mag vs J-K '''
        fig, axes = plt.subplots(2,1,sharex=True)
        ax0,ax1 = axes.flatten()
        ax0.scatter(GAP3['JK'],GAP3['Hmag'],label=r'Full Field')
        ax0.set_title(r'C3')
        ax0.set_ylabel(r'H')
        ax0.scatter(C3orig['JK'],C3orig['Hmag'],label=r'With Detections')
        ax0.legend()
        ax1.scatter(GAP6['JK'],GAP6['Vcut'])
        ax1.set_title(r'C6')
        ax1.set_ylabel(r'V')
        ax1.set_xlabel(r'J-K')
        ax1.scatter(C6orig['JK'],C6orig['Vcut'])
        plt.tight_layout()
        plt.savefig(ext_fig+'K2_Obs_Dets.png')

        plt.figure()
        hist1, xb1, yb1, im1 = plt.hist2d(C3orig['JK'],C3orig['Hmag'],bins=49,cmap=colormaps.parula)#,normed=True)
        hist2, xb2, yb2, im2 = plt.hist2d(GAP3['JK'],GAP3['Hmag'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
        hist = hist2-hist1
        plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
        cbar = plt.colorbar()
        cbar.set_label(r'$N_{\rm{obs}} - N_{\rm{seis}}$', rotation=270, fontsize=20, labelpad=25)
        cbar.ax.tick_params(labelsize=20)
        plt.ylabel(r'H',fontsize=20, labelpad=20)
        plt.xlabel(r'J-K',fontsize=20, labelpad=10)
        plt.ylim(max(yb1),min(yb1))
        plt.title(r'C3 Number Comparison',fontsize=20)
        plt.tick_params(labelsize=15)
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_2d_hist.png')

        plt.figure()
        hist1, xb1, yb1, im1 = plt.hist2d(C6orig['JK'],C6orig['Vcut'],bins=49,cmap=colormaps.parula)#,normed=True)
        hist2, xb2, yb2, im2 = plt.hist2d(GAP6['JK'],GAP6['Vcut'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
        hist = hist2-hist1
        plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
        cbar = plt.colorbar()
        cbar.set_label(r'$N_{\rm{obs}} - N_{\rm{seis}}$', rotation=270, fontsize=20, labelpad=25)
        cbar.ax.tick_params(labelsize=20)
        plt.ylabel(r'V',fontsize=20, labelpad=20)
        plt.xlabel(r'J-K',fontsize=20, labelpad=10)
        plt.ylim(max(yb1),min(yb1))
        plt.title(r'C6 Number Comparison',fontsize=20)
        plt.tick_params(labelsize=15)
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_2d_hist.png')


        ''' Mass '''
        plt.figure()
        hist0, bins, patches = plt.hist(GAP3['mass'],bins=50,label=r'GAP sample')
        plt.hist(C3orig['mass'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'Mass [M$_{\odot}$]')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_mass.png')

        plt.figure()
        hist0, bins, patches = plt.hist(GAP6['mass'],bins=50,label=r'GAP sample')
        plt.hist(C6orig['mass'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'Mass [M$_{\odot}$]')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_mass.png')

        ''' [Fe/H] '''
        plt.figure()
        hist0, bins, patches = plt.hist(GAP3['[Fe/H]'],bins=50,label=r'GAP sample')
        plt.hist(C3orig['[Fe/H]'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'[Fe/H]')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_feh.png')

        plt.figure()
        hist0, bins, patches = plt.hist(GAP6['[Fe/H]'],bins=50,label=r'GAP sample')
        plt.hist(C6orig['[Fe/H]'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'[Fe/H]')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_feh.png')

        ''' Radius '''
        plt.figure()
        hist0, bins, patches = plt.hist(GAP3['Radius'],bins=50,label=r'GAP sample')
        plt.hist(C3orig['Radius'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'Radius [R$_{\odot}$]')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_rad.png')

        plt.figure()
        hist0, bins, patches = plt.hist(GAP6['Radius'],bins=50,label=r'GAP sample')
        plt.hist(C6orig['Radius'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'Radius [R$_{\odot}$]')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_rad.png')

        ''' Teff '''
        plt.figure()
        hist0, bins, patches = plt.hist(GAP3['Teff'],bins=50,label=r'GAP sample')
        plt.hist(C3orig['Teff'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'T$_{\rm{eff}}$ [K]')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_teff.png')

        plt.figure()
        hist0, bins, patches = plt.hist(GAP6['Teff'],bins=50,label=r'GAP sample')
        plt.hist(C6orig['Teff'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'T$_{\rm{eff}}$ [K]')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_teff.png')

        ''' logg '''
        plt.figure()
        hist0, bins, patches = plt.hist(GAP3['logg'],bins=50,label=r'GAP sample')
        plt.hist(C3orig['logg'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'logg')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_logg.png')

        plt.figure()
        hist0, bins, patches = plt.hist(GAP6['logg'],bins=50,label=r'GAP sample')
        plt.hist(C6orig['logg'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'logg')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_logg.png')

        ''' NuMax Comparisons '''
        C3n, C6n = pd.DataFrame(), pd.DataFrame()
        C3n['EPIC'] = GAP3['EPIC']
        C3n['RAD'] = GAP3['Radius']
        C6n['EPIC'] = GAP6['EPIC']
        C6n['RAD'] = GAP6['Radius']
        C3n['Numax'] = GAP3['mass'] * (GAP3['Radius']**-2) * ((GAP3['Teff']/5777)**-0.5) * 3090
        C6n['Numax'] = GAP6['mass'] * (GAP6['Radius']**-2) * ((GAP6['Teff']/5777)**-0.5) * 3090

        C3m = pd.merge(C3n,C3orig,how='inner',on=['EPIC'])
        C3m = C3m.reset_index(drop=True)
        C3m = C3m.fillna(value='NaN',method=None)
        C6m = pd.merge(C6n,C6orig,how='inner',on=['EPIC'])
        C6m = C6m.reset_index(drop=True)
        C6m = C6m.fillna(value='NaN',method=None)

        C3m['diff'] = C3m['NUMAX'] - C3m['Numax']
        C6m['diff'] = C6m['NUMAX'] - C6m['Numax']

        plt.figure()
        plt.scatter(C3m['NUMAX'],C3m['diff'])
        plt.plot([min(C3m['NUMAX']),max(C3m['NUMAX'])],[0,0],color='orange')
        plt.xlabel(r'$\nu_{\rm{max}}$ [$\mu$Hz]')
        plt.ylabel(r'$\nu_{\rm{max}}$ - $\nu_{\rm{max},SR}$')
        plt.xlim(min(C3m['NUMAX']),max(C3m['NUMAX']))
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_numax_comp.png')

        plt.figure()
        plt.scatter(C6m['NUMAX'],C6m['diff'])
        plt.plot([min(C6m['NUMAX']),max(C6m['NUMAX'])],[0,0],color='orange')
        plt.xlabel(r'$\nu_{\rm{max}}$ [$\mu$Hz]')
        plt.ylabel(r'$\nu_{\rm{max}}$ - $\nu_{\rm{max},SR}$')
        plt.xlim(min(C6m['NUMAX']),max(C6m['NUMAX']))
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_numax_comp.png')

        # plt.show()

    ''' Comparisons with TRILEGAL '''
    if sys.argv[1] == '3':

        C3orig = pd.read_csv('/home/bmr135/GA/K2Poles/Selection_Function_Biases/C3_full')
        C6orig = pd.read_csv('/home/bmr135/GA/K2Poles/Selection_Function_Biases/C6_full')
        # C3orig = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3_full')
        # C6orig = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6_full')
        GAP3, GAP6 = dat.K2_GAP()
        T3, T6 = dat.TRILEGAL()

        C3orig = C3orig[ (C3orig['Hmag'] > 7) & (C3orig['Hmag'] < 12) & (C3orig['JK'] > 0.5) & (C3orig['mass'] > 0)]
        C3orig = C3orig.reset_index(drop=True)
        C6orig = C6orig[ (C6orig['Vcut'] > 9) & (C6orig['Vcut'] < 15) & (C6orig['JK'] > 0.5) & (C6orig['mass'] > 0)]
        C6orig = C6orig.reset_index(drop=True)
        GAP3 = GAP3[ (GAP3['Hmag'] > 7) & (GAP3['Hmag'] < 12) & (GAP3['JK'] > 0.5) & (GAP3['mass'] > 0)]
        GAP3 = GAP3.reset_index(drop=True)
        GAP6 = GAP6[ (GAP6['Vcut'] > 9) & (GAP6['Vcut'] < 15) & (GAP6['JK'] > 0.5) & (GAP6['mass'] > 0)]
        GAP6 = GAP6.reset_index(drop=True)
        T3 = T3[ (T3['Hmag'] > 7) & (T3['Hmag'] < 12) & (T3['JK'] > 0.5) & (T3['Mass'] > 0)]
        T3 = T3.reset_index(drop=True)
        T6 = T6[ (T6['Vcut'] > 9) & (T6['Vcut'] < 15) & (T6['JK'] > 0.5) & (T6['Mass'] > 0)]
        T6 = T6.reset_index(drop=True)

        ext_fig = '/home/bmr135/Dropbox/K2Poles/pop_trends/261017/Sim_vs_Obs/'
        # ext_fig = '/home/bmr135/Dropbox/K2Poles/pop_trends/261017/Sim_vs_Obs/'

        ''' Comparison to Whole Samples '''
        fig, axes = plt.subplots(2,1,sharex=True)
        ax0,ax1 = axes.flatten()
        ax0.scatter(T3['JK'],T3['Hmag'],label=r'TRILEGAL')
        ax0.scatter(GAP3['JK'],GAP3['Hmag'],label=r'GAP Field')
        ax0.set_title(r'C3')
        ax0.set_ylabel(r'H')
        ax0.legend()
        ax1.scatter(T6['JK'],T6['Vcut'])
        ax1.scatter(GAP6['JK'],GAP6['Vcut'])
        ax1.set_title(r'C6')
        ax1.set_ylabel(r'V')
        ax1.set_xlabel(r'J-K')
        plt.tight_layout()
        plt.savefig(ext_fig+'K2_comp.png')


        ''' Comparison to Sample with Detections '''
        fig, axes = plt.subplots(2,1,sharex=True)
        ax0,ax1 = axes.flatten()
        ax0.scatter(T3['JK'],T3['Hmag'],label=r'TRILEGAL')
        ax0.scatter(C3orig['JK'],C3orig['Hmag'],label=r'With Dets.')
        ax0.set_title(r'C3')
        ax0.set_ylabel(r'H')
        ax0.legend()
        ax1.scatter(T6['JK'],T6['Vcut'])
        ax1.scatter(C6orig['JK'],C6orig['Vcut'])
        ax1.set_title(r'C6')
        ax1.set_ylabel(r'V')
        ax1.set_xlabel(r'J-K')
        plt.tight_layout()
        plt.savefig(ext_fig+'K2_sim_dets.png')

        ''' 2D Histograms '''
        plt.figure()
        hist1, xb1, yb1, im1 = plt.hist2d(GAP3['JK'],GAP3['Hmag'],bins=49,cmap=colormaps.parula)#,normed=True)
        hist2, xb2, yb2, im2 = plt.hist2d(T3['JK'],T3['Hmag'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
        hist = hist2-hist1
        plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
        cbar = plt.colorbar()
        cbar.set_label(r'$N_{\rm{sim}} - N_{\rm{GAP}}$', rotation=270, fontsize=20, labelpad=25)
        cbar.ax.tick_params(labelsize=20)
        plt.ylabel(r'H',fontsize=20, labelpad=20)
        plt.xlabel(r'J-K',fontsize=20, labelpad=10)
        plt.ylim(max(yb1),min(yb1))
        plt.title(r'C3 Number Comparison',fontsize=20)
        plt.tick_params(labelsize=15)
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_2d_hist_GAP_sim.png')

        plt.figure()
        hist1, xb1, yb1, im1 = plt.hist2d(GAP6['JK'],GAP6['Vcut'],bins=49,cmap=colormaps.parula)#,normed=True)
        hist2, xb2, yb2, im2 = plt.hist2d(T6['JK'],T6['Vcut'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
        hist = hist2-hist1
        plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
        cbar = plt.colorbar()
        cbar.set_label(r'$N_{\rm{sim}} - N_{\rm{GAP}}$', rotation=270, fontsize=20, labelpad=25)
        cbar.ax.tick_params(labelsize=20)
        plt.ylabel(r'Vcut',fontsize=20, labelpad=20)
        plt.xlabel(r'J-K',fontsize=20, labelpad=10)
        plt.ylim(max(yb1),min(yb1))
        plt.title(r'C6 Number Comparison',fontsize=20)
        plt.tick_params(labelsize=15)
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_2d_hist_GAP_sim.png')

        plt.figure()
        hist1, xb1, yb1, im1 = plt.hist2d(C3orig['JK'],C3orig['Hmag'],bins=49,cmap=colormaps.parula)#,normed=True)
        hist2, xb2, yb2, im2 = plt.hist2d(T3['JK'],T3['Hmag'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
        hist = hist2-hist1
        plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
        cbar = plt.colorbar()
        cbar.set_label(r'$N_{\rm{sim}} - N_{\rm{seis}}$', rotation=270, fontsize=20, labelpad=25)
        cbar.ax.tick_params(labelsize=20)
        plt.ylabel(r'H',fontsize=20, labelpad=20)
        plt.xlabel(r'J-K',fontsize=20, labelpad=10)
        plt.ylim(max(yb1),min(yb1))
        plt.title(r'C3 Number Comparison',fontsize=20)
        plt.tick_params(labelsize=15)
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_2d_hist_dets_sim.png')

        plt.figure()
        hist1, xb1, yb1, im1 = plt.hist2d(C6orig['JK'],C6orig['Vcut'],bins=49,cmap=colormaps.parula)#,normed=True)
        hist2, xb2, yb2, im2 = plt.hist2d(T6['JK'],T6['Vcut'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
        hist = hist2-hist1
        plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
        cbar = plt.colorbar()
        cbar.set_label(r'$N_{\rm{sim}} - N_{\rm{seis}}$', rotation=270, fontsize=20, labelpad=25)
        cbar.ax.tick_params(labelsize=20)
        plt.ylabel(r'Vcut',fontsize=20, labelpad=20)
        plt.xlabel(r'J-K',fontsize=20, labelpad=10)
        plt.ylim(max(yb1),min(yb1))
        plt.title(r'C6 Number Comparison',fontsize=20)
        plt.tick_params(labelsize=15)
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_2d_hist_dets_sim.png')

        ''' TRILEGAL Prior to Det. Prob. Cut '''
        plt.figure()
        hist0, bins, patches = plt.hist(T3['Teff'],bins=50,label=r'TRILEGAL')
        plt.hist(GAP3['Teff'],bins=bins,label=r'GAP')
        plt.xlabel(r'T$_{\rm{eff}}$ [K]')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_teff_no_cut.png')

        plt.figure()
        hist0, bins, patches = plt.hist(T6['Teff'],bins=50,label=r'TRILEGAL')
        plt.hist(GAP6['Teff'],bins=bins,label=r'GAP')
        plt.xlabel(r'T$_{\rm{eff}}$ [K]')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_teff_no_cut.png')

        plt.figure()
        hist0, bins, patches = plt.hist(T3['logg'],bins=50,label=r'TRILEGAL')
        plt.hist(GAP3['logg'],bins=bins,label=r'GAP')
        plt.xlabel(r'logg')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_logg_no_cut.png')

        plt.figure()
        hist0, bins, patches = plt.hist(T6['logg'],bins=50,label=r'TRILEGAL')
        plt.hist(GAP6['logg'],bins=bins,label=r'GAP')
        plt.xlabel(r'logg')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_logg_no_cut.png')

        ''' Detection Probability Cut '''
        T3 = prop.det_prob(T3,'numax',3090.0,135.1)
        T6 = prop.det_prob(T6,'numax',3090.0,135.1)

        C3orig = C3orig[C3orig['prob_s'] > 0.95]
        T3 = T3[T3['prob_s'] > 0.95]
        C6orig = C6orig[C6orig['prob_s'] > 0.95]
        T6 = T6[T6['prob_s'] > 0.95]

        plt.figure()
        hist1, xb1, yb1, im1 = plt.hist2d(C3orig['JK'],C3orig['Hmag'],bins=49,cmap=colormaps.parula)#,normed=True)
        hist2, xb2, yb2, im2 = plt.hist2d(T3['JK'],T3['Hmag'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
        hist = hist2-hist1
        plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
        cbar = plt.colorbar()
        cbar.set_label(r'$N_{\rm{sim}} - N_{\rm{seis}}$', rotation=270, fontsize=20, labelpad=25)
        cbar.ax.tick_params(labelsize=20)
        plt.ylabel(r'H',fontsize=20, labelpad=20)
        plt.xlabel(r'J-K',fontsize=20, labelpad=10)
        plt.ylim(max(yb1),min(yb1))
        plt.title(r'C3 Number Comparison (det. prob.)',fontsize=20)
        plt.tick_params(labelsize=15)
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_2d_hist_det_prob.png')

        plt.figure()
        hist1, xb1, yb1, im1 = plt.hist2d(C6orig['JK'],C6orig['Vcut'],bins=49,cmap=colormaps.parula)#,normed=True)
        hist2, xb2, yb2, im2 = plt.hist2d(T6['JK'],T6['Vcut'],bins=[xb1,yb1],cmap=colormaps.parula)#,normed=True)
        hist = hist2-hist1
        plt.imshow(hist.T,interpolation='none',cmap=colormaps.parula,extent=[min(xb1),max(xb1),max(yb1),min(yb1)],aspect='auto')
        cbar = plt.colorbar()
        cbar.set_label(r'$N_{\rm{sim}} - N_{\rm{seis}}$', rotation=270, fontsize=20, labelpad=25)
        cbar.ax.tick_params(labelsize=20)
        plt.ylabel(r'Vcut',fontsize=20, labelpad=20)
        plt.xlabel(r'J-K',fontsize=20, labelpad=10)
        plt.ylim(max(yb1),min(yb1))
        plt.title(r'C6 Number Comparison (det. prob.)',fontsize=20)
        plt.tick_params(labelsize=15)
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_2d_hist_det_prob.png')

        ''' Individual Parameters '''
        ''' Mass '''
        plt.figure()
        hist0, bins, patches = plt.hist(T3['Mass'],bins=50,label=r'TRILEGAL')
        plt.hist(C3orig['mass'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'Mass [M$_{\odot}$]')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_mass.png')

        plt.figure()
        hist0, bins, patches = plt.hist(T6['Mass'],bins=50,label=r'TRILEGAL')
        plt.hist(C6orig['mass'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'Mass [M$_{\odot}$]')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_mass.png')

        ''' [Fe/H] '''
        plt.figure()
        hist0, bins, patches = plt.hist(T3['M_H'],bins=50,label=r'TRILEGAL')
        plt.hist(C3orig['[Fe/H]'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'[Fe/H]')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_feh.png')

        plt.figure()
        hist0, bins, patches = plt.hist(T6['M_H'],bins=50,label=r'TRILEGAL')
        plt.hist(C6orig['[Fe/H]'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'[Fe/H]')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_feh.png')

        ''' Radius '''
        plt.figure()
        hist0, bins, patches = plt.hist(T3['radius'],bins=50,label=r'TRILEGAL')
        plt.hist(C3orig['Radius'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'Radius [R$_{\odot}$]')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_rad.png')

        plt.figure()
        hist0, bins, patches = plt.hist(T6['radius'],bins=50,label=r'TRILEGAL')
        plt.hist(C6orig['Radius'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'Radius [R$_{\odot}$]')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_rad.png')

        ''' Teff '''
        plt.figure()
        hist0, bins, patches = plt.hist(T3['Teff'],bins=50,label=r'TRILEGAL')
        plt.hist(C3orig['Teff'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'T$_{\rm{eff}}$ [K]')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_teff.png')

        plt.figure()
        hist0, bins, patches = plt.hist(T6['Teff'],bins=50,label=r'TRILEGAL')
        plt.hist(C6orig['Teff'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'T$_{\rm{eff}}$ [K]')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_teff.png')

        ''' logg '''
        plt.figure()
        hist0, bins, patches = plt.hist(T3['logg'],bins=50,label=r'TRILEGAL')
        plt.hist(C3orig['logg'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'logg')
        plt.legend()
        plt.title(r'C3')
        plt.tight_layout()
        plt.savefig(ext_fig+'C3_logg.png')

        plt.figure()
        hist0, bins, patches = plt.hist(T6['logg'],bins=50,label=r'TRILEGAL')
        plt.hist(C6orig['logg'],bins=bins,label=r'With Dets.')
        plt.xlabel(r'logg')
        plt.legend()
        plt.title(r'C6')
        plt.tight_layout()
        plt.savefig(ext_fig+'C6_logg.png')

        ''' Combined Plot '''
        fig, axes = plt.subplots(2,1)
        T3 = T3[T3['Mass'] < 4]
        T3 = T3[T3['radius'] < 20]
        C3orig = C3orig[C3orig['mass'] < 4]
        C3orig = C3orig[C3orig['Radius'] < 20]
        ax0,ax1 = axes.flatten()
        hist0, bins, patches = ax0.hist(T3['radius'],histtype='stepped',bins=50,label=r'TRILEGAL',normed=True)
        ax0.hist(C3orig['Radius'],bins=bins,histtype='stepped',label=r'With Dets.',normed=True)
        ax0.set_xlabel(r'Radius [R$_{\odot}$]')
        ax0.legend()
        hist1, bins1, patches1 = ax1.hist(T3['Mass'],bins=50,histtype='step',label=r'TRILEGAL',normed=True)
        ax1.hist(C3orig['mass'],bins=bins1,histtype='step',label=r'With Dets.',normed=True)
        ax1.set_xlabel(r'Mass [M$_{\odot}$]')
        # fig.subplots_adjust(top=3.0)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])


        plt.show()

    ''' Statistics - Sim and Obs Comparisons '''
    if sys.argv[1] == '4':
        ''' Data and Probability Cut '''
        C3orig = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C3_full')
        C6orig = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/Selection_Function_Biases/C6_full')
        GAP3, GAP6 = dat.K2_GAP()
        T3, T6 = dat.TRILEGAL()

        C3orig = C3orig[ (C3orig['Hmag'] > 7) & (C3orig['Hmag'] < 12) & (C3orig['JK'] > 0.5) & (C3orig['mass'] > 0)]
        C3orig = C3orig.reset_index(drop=True)
        C6orig = C6orig[ (C6orig['Vcut'] > 9) & (C6orig['Vcut'] < 15) & (C6orig['JK'] > 0.5) & (C6orig['mass'] > 0)]
        C6orig = C6orig.reset_index(drop=True)
        GAP3 = GAP3[ (GAP3['Hmag'] > 7) & (GAP3['Hmag'] < 12) & (GAP3['JK'] > 0.5) & (GAP3['mass'] > 0)]
        GAP3 = GAP3.reset_index(drop=True)
        GAP6 = GAP6[ (GAP6['Vcut'] > 9) & (GAP6['Vcut'] < 15) & (GAP6['JK'] > 0.5) & (GAP6['mass'] > 0)]
        GAP6 = GAP6.reset_index(drop=True)
        T3 = T3[ (T3['Hmag'] > 7) & (T3['Hmag'] < 12) & (T3['JK'] > 0.5) & (T3['Mass'] > 0)]
        T3 = T3.reset_index(drop=True)
        T6 = T6[ (T6['Vcut'] > 9) & (T6['Vcut'] < 15) & (T6['JK'] > 0.5) & (T6['Mass'] > 0)]
        T6 = T6.reset_index(drop=True)

        x0,x1,x2,x3,x4 = pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
        x0 = T3[T3['label'] == 0]
        x1 = T3[T3['label'] == 1]
        x2 = T3[T3['label'] == 2]
        x3 = T3[T3['label'] == 3]
        x4 = T3[T3['label'] == 4]

        y0,y1,y2,y6,y4 = pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
        y0 = T6[T6['label'] == 0]
        y1 = T6[T6['label'] == 1]
        y2 = T6[T6['label'] == 2]
        y3 = T6[T6['label'] == 3]
        y4 = T6[T6['label'] == 4]

        ''' Comparison to Sample with Detections '''
        fig, axes = plt.subplots(2,1,sharex=True)
        ax0,ax1 = axes.flatten()
        ax0.scatter(x4['JK'],x4['Hmag'],label=r'TRILEGAL 4')
        ax0.scatter(x3['JK'],x3['Hmag'],label=r'TRILEGAL 3')
        ax0.scatter(x2['JK'],x2['Hmag'],label=r'TRILEGAL 2')
        ax0.scatter(x1['JK'],x1['Hmag'],label=r'TRILEGAL 1')
        ax0.scatter(x0['JK'],x0['Hmag'],label=r'TRILEGAL 0')
        # ax0.scatter(C3orig['JK'],C3orig['Hmag'],label=r'With Dets.')
        ax0.set_title(r'C3')
        ax0.set_ylabel(r'H')
        ax0.legend()
        ax1.scatter(y4['JK'],y4['Vcut'],label=r'TRILEGAL 4')
        ax1.scatter(y3['JK'],y3['Vcut'],label=r'TRILEGAL 3')
        ax1.scatter(y2['JK'],y2['Vcut'],label=r'TRILEGAL 2')
        ax1.scatter(y1['JK'],y1['Vcut'],label=r'TRILEGAL 1')
        ax1.scatter(y0['JK'],y0['Vcut'],label=r'TRILEGAL 0')
        # ax1.scatter(C6orig['JK'],C6orig['Vcut'])
        ax1.set_title(r'C6')
        ax1.set_ylabel(r'V')
        ax1.set_xlabel(r'J-K')
        plt.tight_layout()
        # plt.savefig(ext_fig+'K2_sim_dets.png')

        T3 = prop.det_prob(T3,'numax',3090.0,135.1)
        T6 = prop.det_prob(T6,'numax',3090.0,135.1)

        C3orig = C3orig[C3orig['prob_s'] > 0.95]
        C3orig = C3orig.reset_index(drop=True)
        T3 = T3[T3['prob_s'] > 0.95]
        T3 = T3.reset_index(drop=True)
        C6orig = C6orig[C6orig['prob_s'] > 0.95]
        C6orig = C6orig.reset_index(drop=True)
        T6 = T6[T6['prob_s'] > 0.95]
        T6 = T6.reset_index(drop=True)

        x0,x1,x2,x3,x4 = pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
        x0 = T3[T3['label'] == 0]
        x1 = T3[T3['label'] == 1]
        x2 = T3[T3['label'] == 2]
        x3 = T3[T3['label'] == 3]
        x4 = T3[T3['label'] == 4]

        y0,y1,y2,y6,y4 = pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
        y0 = T6[T6['label'] == 0]
        y1 = T6[T6['label'] == 1]
        y2 = T6[T6['label'] == 2]
        y3 = T6[T6['label'] == 3]
        y4 = T6[T6['label'] == 4]

        ''' Comparison to Sample with Detections '''
        fig, axes = plt.subplots(2,1,sharex=True)
        ax0,ax1 = axes.flatten()
        ax0.scatter(x4['JK'],x4['Hmag'],label=r'TRILEGAL 4')
        ax0.scatter(x3['JK'],x3['Hmag'],label=r'TRILEGAL 3')
        ax0.scatter(x2['JK'],x2['Hmag'],label=r'TRILEGAL 2')
        ax0.scatter(x1['JK'],x1['Hmag'],label=r'TRILEGAL 1')
        ax0.scatter(x0['JK'],x0['Hmag'],label=r'TRILEGAL 0')
        ax0.scatter(C3orig['JK'],C3orig['Hmag'],label=r'With Dets.')
        ax0.set_title(r'C3')
        ax0.set_ylabel(r'H')
        ax0.legend()
        ax1.scatter(y4['JK'],y4['Vcut'],label=r'TRILEGAL 4')
        ax1.scatter(y3['JK'],y3['Vcut'],label=r'TRILEGAL 3')
        ax1.scatter(y2['JK'],y2['Vcut'],label=r'TRILEGAL 2')
        ax1.scatter(y1['JK'],y1['Vcut'],label=r'TRILEGAL 1')
        ax1.scatter(y0['JK'],y0['Vcut'],label=r'TRILEGAL 0')
        ax1.scatter(C6orig['JK'],C6orig['Vcut'],label=r'With Dets.')
        ax1.set_title(r'C6')
        ax1.set_ylabel(r'V')
        ax1.set_xlabel(r'J-K')
        plt.tight_layout()
        # plt.savefig(ext_fig+'K2_sim_dets.png'


        ''' Calculate Mean, Median, Mode and Standard Deviation '''
        ''' Mass '''
        frames = [C3orig,C6orig,T3,T6]
        mass = pd.DataFrame()
        mass['Data'] = ['C3','C6','T3','T6']
        mass['Mean'], mass['Median'], mass['Mode'], mass['Std_Dev'] = 0.0, 0.0, 0.0, 0.0
        a = ['mass','mass','Mass','Mass']

        for i in range(4):
            X = frames[i]
            mass['Mean'][i] = X[a[i]].mean()
            mass['Median'][i] = X[a[i]].median()
            # mass['Mode'][i] = X['mass'].mode()
            mass['Std_Dev'][i] = X[a[i]].std()
            frames[i] = X
        # print(mass)

        ''' [Fe/H] '''
        feh = pd.DataFrame()
        feh['Data'] = ['C3','C6','T3','T6']
        feh['Mean'], feh['Median'], feh['Mode'], feh['Std_Dev'] = 0.0, 0.0, 0.0, 0.0
        a = ['[Fe/H]','[Fe/H]','M_H','M_H']

        for i in range(4):
            X = frames[i]
            feh['Mean'][i] = X[a[i]].mean()
            feh['Median'][i] = X[a[i]].median()
            # [Fe/H]['Mode'][i] = X['[Fe/H]'].mode()
            feh['Std_Dev'][i] = X[a[i]].std()
            frames[i] = X
        # print(feh)

        ''' Teff '''
        Teff = pd.DataFrame()
        Teff['Data'] = ['C3','C6','T3','T6']
        Teff['Mean'], Teff['Median'], Teff['Mode'], Teff['Std_Dev'] = 0.0, 0.0, 0.0, 0.0
        a = ['Teff','Teff','Teff','Teff']

        for i in range(4):
            X = frames[i]
            Teff['Mean'][i] = X[a[i]].mean()
            Teff['Median'][i] = X[a[i]].median()
            # Teff['Mode'][i] = X['Teff'].mode()
            Teff['Std_Dev'][i] = X[a[i]].std()
            frames[i] = X
        # print(Teff)

        ''' Radius '''
        radius = pd.DataFrame()
        radius['Data'] = ['C3','C6','T3','T6']
        radius['Mean'], radius['Median'], radius['Mode'], radius['Std_Dev'] = 0.0, 0.0, 0.0, 0.0
        a = ['Radius','Radius','radius','radius']

        for i in range(4):
            X = frames[i]
            radius['Mean'][i] = X[a[i]].mean()
            radius['Median'][i] = X[a[i]].median()
            # radius['Mode'][i] = X['radius'].mode()
            radius['Std_Dev'][i] = X[a[i]].std()
            frames[i] = X
        # print(radius)

        ''' logg '''
        logg = pd.DataFrame()
        logg['Data'] = ['C3','C6','T3','T6']
        logg['Mean'], logg['Median'], logg['Mode'], logg['Std_Dev'] = 0.0, 0.0, 0.0, 0.0
        a = ['logg','logg','logg','logg']

        for i in range(4):
            X = frames[i]
            logg['Mean'][i] = X[a[i]].mean()
            logg['Median'][i] = X[a[i]].median()
            # logg['Mode'][i] = X['logg'].mode()
            logg['Std_Dev'][i] = X[a[i]].std()
            frames[i] = X
        # print(logg)

        ''' Likelihood Calculations '''
        # C3orig['lnP_m'],C3orig['lnP_f'],C3orig['lnP_t'],C3orig['lnP_r'],C3orig['lnP_g'] = 0,0,0,0,0
        # C6orig['lnP_m'],C6orig['lnP_f'],C6orig['lnP_t'],C6orig['lnP_r'],C6orig['lnP_g'] = 0,0,0,0,0
        T3['lnP_m'],T3['lnP_f'],T3['lnP_t'],T3['lnP_r'],T3['lnP_g'] = 0,0,0,0,0
        T6['lnP_m'],T6['lnP_f'],T6['lnP_t'],T6['lnP_r'],T6['lnP_g'] = 0,0,0,0,0

        T3 = likeli(T3,mass['Mean'][0],mass['Std_Dev'][0],'Mass','lnP_m')
        T6 = likeli(T6,mass['Mean'][1],mass['Std_Dev'][1],'Mass','lnP_m')

        T3 = likeli(T3,radius['Mean'][0],radius['Std_Dev'][0],'radius','lnP_r')
        T6 = likeli(T6,radius['Mean'][1],radius['Std_Dev'][1],'radius','lnP_r')

        T3 = likeli(T3,feh['Mean'][0],feh['Std_Dev'][0],'M_H','lnP_f')
        T6 = likeli(T6,feh['Mean'][1],feh['Std_Dev'][1],'M_H','lnP_f')

        T3 = likeli(T3,Teff['Mean'][0],Teff['Std_Dev'][0],'Teff','lnP_t')
        T6 = likeli(T6,Teff['Mean'][1],Teff['Std_Dev'][1],'Teff','lnP_t')

        T3 = likeli(T3,logg['Mean'][0],logg['Std_Dev'][0],'logg','lnP_g')
        T6 = likeli(T6,logg['Mean'][1],logg['Std_Dev'][1],'logg','lnP_g')

        # likeli_plt(T3,C3orig,'Mass','mass','lnP_m')
        # likeli_plt(T3,C3orig,'Teff','Teff','lnP_t')
        # likeli_plt(T3,C3orig,'radius','Radius','lnP_r')
        # likeli_plt(T3,C3orig,'M_H','[Fe/H]','lnP_f')
        # likeli_plt(T3,C3orig,'logg','logg','lnP_g')

        plt.show()
