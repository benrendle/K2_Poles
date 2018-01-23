''' Calculation of different properties required in
    K2_seismo_comp.py '''

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import sys
from pandas import DataFrame, read_csv
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import scipy.odr.odrpack as odrpack
import numpy as np
from scipy.stats import norm
import Detection_prob_K2 as DP
import Detection_prob_kep as DP2
import K2_constants as const
''' Credit to M. Schofield for creation of the detection recipe code, 01/2017 '''

def selection_function(data,numax):
    ''' Selection function used for K2 GAP and unbiasing the
        asteroseismic dectections. Based on proposals by Dennis Stello for
        the C3 and C6 campaigns '''
    for i in range(0,len(data),1):
        df = data[i]
        if i < 8: # C3
            df = df[ (df[numax[i]] > 10) & (df[numax[i]] < 280) & (df['Hmag'] > 7) & (df['Hmag'] < 12) & (df['JK'] > 0.5) & (df['prob_s'] >= 0.95) ]
            # df = df[ (df['Hmag'] > 7) & (df['Hmag'] < 12) & (df['JK'] > 0.5) ]
        elif i >= 8: # C6
            df = df[ (df[numax[i]] > 10) & (df[numax[i]] < 280) & (df['Vcut'] > 9) & (df['Vcut'] < 15) & (df['JK'] > 0.5) & (df['prob_s'] >= 0.95) ]
            # df = df[ (df['Vcut'] > 9) & (df['Vcut'] < 15) & (df['JK'] > 0.5) ]

        df = df.reset_index(drop=True)
        data[i] = df

    return data

def lmrl_comps(data,numax,dnu,Numax,Dnu,GAP):#,Numax,Dnu,Teff,g):

    ''' Function to calculate logg, mass, radius and luminosity
    Numax and Dnu are relative solar values for each pipeline '''

    for i in range(0,len(data),1):
        X = data[i]
        X['slogg'] = np.log10(const.solar_g * (X[numax[i]]/Numax[i]) * np.sqrt(X['Teff']/const.solar_Teff))
        X['mass'] = (X[numax[i]]/Numax[i])**3 * (X[dnu[i]]/Dnu[i])**-4 * (X['Teff']/const.solar_Teff)**1.5
        X['radius'] = (X[numax[i]]/Numax[i]) * (X[dnu[i]]/Dnu[i])**-2 * (X['Teff']/const.solar_Teff)**0.5
        X['L'] = X['radius']**2 * (X['Teff']/const.solar_Teff)**4
        if GAP == 1:
            X['e_Teff'] = np.sqrt(X['em_teff']**2 + X['ep_teff']**2 + 59**2) # 59K from Torres et al. 2012
            X['e_FeH'] = np.sqrt(X['em_[Fe/H]']**2 + X['ep_[Fe/H]']**2 + 0.06**2) # 0.06K from Torres et al. 2012
            X['e_logg'] = np.sqrt(X['em_logg']**2 + X['ep_logg']**2) # 59K from Torres et al. 2012
        # X['Tred'] = DP.properties(X, constants_only=False)
        X = det_prob(X,numax[i],Numax[i],Dnu[i])
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
    X['prob_s'], X['SNR'] = DP.globalDetections(X['imag'],X['KepMag'],X['L'],X['radius'],X['Teff'],X[numax],1.0,X['Tred'], \
                    const.teffred_solar,const.solar_Teff,Numax,Dnu,const.sys_limit,const.dilution,const.vnyq, \
                    const.cadence, vary_beta=True)
    return X

def det_prob_GAP(X,numax,Numax,Dnu):
    ''' Run detection probability code from Mat, calculating teffred to
    start with and then running the global detection code with added
    noise. Output is the detection probability and SNR

    numax = stellar numax
    Numax = calibration numax
    Dnu = calibration numax '''

    X['Tred'] = DP.properties_GAP(X, constants_only=False)
    X['prob_s'], X['SNR'] = DP.globalDetections(X['imag'],X['KepMag'],X['Lumo'],X['Radius'],X['Teff'],X[numax],1.0,X['Tred'], \
                    const.teffred_solar,const.solar_Teff,Numax,Dnu,const.sys_limit,const.dilution,const.vnyq, \
                    const.cadence,vary_beta=True)
    return X

def det_prob_Kepler(X,numax,Numax,Dnu):
    ''' Run detection probability code from Mat, calculating teffred to
    start with and then running the global detection code with added
    noise. Output is the detection probability and SNR

    numax = stellar numax
    Numax = calibration numax
    Dnu = calibration numax '''

    X['Tred'] = DP2.properties(X, constants_only=False)
    X['prob_s'], X['SNR'] = DP2.globalDetections(X['imag'],X['KepMag'],X['L'],X['radius'],X['Teff'],X[numax],1.0,X['Tred'], \
                    const.teffred_solar,const.solar_Teff,Numax,Dnu,const.sys_limit,const.dilution,const.vnyq, \
                    const.cadence, vary_beta=True)
    return X

def det_fract(data):
    ''' Dectection fraction for a given bin of magnitude '''
    bins=[0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9]
    data['binned'] = np.digitize(data.JK.values, bins=bins)

    ''' loop where you iterate over the number of entries in in bins and
    calculate the fraction of detections from there '''
    z=[]
    for i in range(0,len(bins)-1,1):
        x = data[data['prob_s'] > 0.95]
        # print( len(x))
        x = x[x['binned'] == i]
        y = data[data['binned'] == i]
        # print( len(y))
        if len(y) == 0:
            fract = 0.0
        elif len(y) > 0:
            fract = np.float(len(x))/np.float(len(y))
        z.append(fract)
    # print( z)
    return z

def galactic_coords(data):
    ''' Function to convert RA and DEC to galactic coordinates l and b.
    This function requires an input list of data frames to work '''
    for i in data:
        c_icrs = SkyCoord(ra=i['RA'].as_matrix()*u.degree, dec=i['Dec'].as_matrix()*u.degree, frame='icrs')
        c = c_icrs.galactic
        d = coord.ICRS(ra=i['RA'].as_matrix()*u.degree, dec=i['Dec'].as_matrix()*u.degree, distance=i['Distance'].as_matrix()*u.pc)
        e = d.transform_to(coord.Galactocentric)
        i['Glon'] = c.l*u.degree # l is galactic longitude
        i['Glat'] = c.b*u.degree # b is galactic latitiude
        i['Grad'] = np.sqrt((e.x*u.pc)**2 + (e.y*u.pc)**2)# + (e.z*u.pc)**2) <- Add for distance
        i['GZ'] = e.z*u.pc
    return data

def galactic_coords2(data):
    ''' Function to convert RA and DEC to galactic coordinates l and b.
    The is method will work for a single data frame '''

    c_icrs = SkyCoord(ra=data['RA'].as_matrix()*u.degree, dec=data['Dec'].as_matrix()*u.degree, frame='icrs')
    c = c_icrs.galactic

    data['Glon'] = c.l*u.degree # l is galactic longitude
    data['Glat'] = c.b*u.degree # b is galactic latitiude

    return data

def least_squares_2err(Matched):
    ''' Least squares fitting routine '''
    def f(Par,x):
        return Par[0]*x + Par[1]
    # runs=200
    # samp_size=len(Matched['mass_x'])/5#int(np.sqrt(len(x)))
    # Matched['xerr'] = (Matched['e_mass1_x'] + Matched['e_mass2_x'])/2.0
    # Matched['yerr'] = (Matched['e_mass1_y'] + Matched['e_mass2_y'])/2.0
    mpar,cpar=[],[] # y = mx + c
    empar,ecpar=[],[]
    # for i in range(0,runs,1):
    #     # idx=np.random.choice(len(x),samp_size,replace=False)
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(Matched['BDnu'], Matched['EDnu'], sx=Matched['e_BDnu'], sy=Matched['e_EDnu'])#, sy=Matched['yerr'])
    myodr = odrpack.ODR(mydata, linear, beta0=[1., 0.],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    myoutput.pprint()
    # print( mydata)
    return mpar, cpar, empar, ecpar

def least_squares_1err(Matched):
    ''' Least squares fitting routine '''
    def f(Par,x):
        return Par[0]*x + Par[1]
    # runs=200
    # samp_size=len(Matched['mass_x'])/5#int(np.sqrt(len(x)))
    Matched['yerr'] = (Matched['e_mass1'] + Matched['e_mass2'])/2.0
    mpar,cpar=[],[] # y = mx + c
    empar,ecpar=[],[]
    # for i in range(0,runs,1):
    #     # idx=np.random.choice(len(x),samp_size,replace=False)
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(Matched['mass_x'], Matched['mass_y'], sy=Matched['yerr'])#, sx=Matched['xerr']
    myodr = odrpack.ODR(mydata, linear, beta0=[1., 0.],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    myoutput.pprint()
    # print( mydata)
    return mpar, cpar, empar, ecpar

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
    Z = pd.concat([X,Y],ignore_index=True)
    Z = Z.drop_duplicates(subset=['EPIC'],keep=False)
    Z = Z.reset_index(drop=True)
    X = Y

    return X, Z

def temp_JK(I,met,JK):
    ''' Calculate temperature from J-K, Alonso 1999 + 2001 erratum
        theta = Teff/5040'''
    a0 = 0.5816
    a1 = 0.9134
    a2 = -0.1443
    a3 = 0.0
    a4 = 0.0
    a5 = 0.0
    X = JK
    F = met

    theta = a0 + a1*X + a2*X**2 + a3*F + a4*F + a5*F**2
    I['Teff_JK'] = 5040/theta

    return I

def individ(df,df1,df2,pipeline):
    ''' Generate outfile containing EPIC numbers that only appear in single
        datasets '''
    Y = pd.DataFrame([])
    Y['EPIC'] = df['EPIC']
    Y['Y'] = 1.0
    Y['B'] = 0.0
    Y['S'] = 0.0
    # Y['YDnu'] = df['dnu']
    # Y['e_YDnu'] = df['dnu_err']
    # Y['Ynumax'] = df['nmx']
    # Y['e_Ynumax'] = df['nmx_err']
    # Y['BDnu'] = -99.9
    # Y['e_BDnu'] = -99.9
    # Y['Bnumax'] = -99.9
    # Y['e_Bnumax'] = -99.9
    # Y['SDnu'] = -99.9
    # Y['e_SDnu'] = -99.9
    # Y['Snumax'] = -99.9
    # Y['e_Snumax'] = -99.9

    B = pd.DataFrame([])
    B['EPIC'] = df1['EPIC']
    B['Y'] = 0.0
    B['B'] = 1.0
    B['S'] = 0.0
    # B['BDnu'] = df1['EDnu']
    # B['e_BDnu'] = df1['e_EDnu']
    # B['Bnumax'] = df1['Enumax']
    # B['e_Bnumax'] = df1['e_Enumax']
    # B['YDnu'] = -99.9
    # B['e_YDnu'] = -99.9
    # B['Ynumax'] = -99.9
    # B['e_Ynumax'] = -99.9
    # B['SDnu'] = -99.9
    # B['e_SDnu'] = -99.9
    # B['Snumax'] = -99.9
    # B['e_Snumax'] = -99.9

    S = pd.DataFrame([])
    S['EPIC'] = df2['EPIC']
    S['Y'] = 0.0
    S['B'] = 0.0
    S['S'] = 1.0
    # S['SDnu'] = df2['SDnu']
    # S['e_SDnu'] = df2['e_SDnu']
    # S['Snumax'] = df2['Snumax']
    # S['e_Snumax'] = df2['e_Snumax']
    # S['BDnu'] = -99.9
    # S['e_BDnu'] = -99.9
    # S['Bnumax'] = -99.9
    # S['e_Bnumax'] = -99.9
    # S['YDnu'] = -99.9
    # S['e_YDnu'] = -99.9
    # S['Ynumax'] = -99.9
    # S['e_Ynumax'] = -99.9

    for i in range(0,len(Y['EPIC']),1):
        for j in range(0,len(B['EPIC']),1):
            if Y['EPIC'][i] == B['EPIC'][j]:
                Y['B'][i] = 1
                # Y['BDnu'][i] = df1['EDnu'][j]
                # Y['e_BDnu'][i] = df1['e_EDnu'][j]
                # Y['Bnumax'][i] = df1['Enumax'][j]
                # Y['e_Bnumax'][i] = df1['e_Enumax'][j]
        for k in range(0,len(S['EPIC']),1):
            if Y['EPIC'][i] == S['EPIC'][k]:
                Y['S'][i] = 1
                # Y['SDnu'][i] = df2['SDnu'][k]
                # Y['e_SDnu'][i] = df2['e_SDnu'][k]
                # Y['Snumax'][i] = df2['Snumax'][k]
                # Y['e_Snumax'][i] = df2['e_Snumax'][k]

    for i in range(0,len(B['EPIC']),1):
        for j in range(0,len(Y['EPIC']),1):
            if B['EPIC'][i] == Y['EPIC'][j]:
                B['Y'][i] = 1
                # B['YDnu'][i] = df['dnu'][k]
                # B['e_YDnu'][i] = df['dnu_err'][k]
                # B['Ynumax'][i] = df['nmx'][k]
                # B['e_Ynumax'][i] = df['nmx_err'][k]
        for k in range(0,len(S['EPIC']),1):
            if B['EPIC'][i] == S['EPIC'][k]:
                B['S'][i] = 1
                # B['SDnu'][i] = df2['SDnu'][k]
                # B['e_SDnu'][i] = df2['e_SDnu'][k]
                # B['Snumax'][i] = df2['Snumax'][k]
                # B['e_Snumax'][i] = df2['e_Snumax'][k]

    for i in range(0,len(S['EPIC']),1):
        for j in range(0,len(B['EPIC']),1):
            if S['EPIC'][i] == B['EPIC'][j]:
                S['B'][i] = 1
                # S['BDnu'][i] = df1['EDnu'][j]
                # S['e_BDnu'][i] = df1['e_EDnu'][j]
                # S['Bnumax'][i] = df1['Enumax'][j]
                # S['e_Bnumax'][i] = df1['e_Enumax'][j]
        for k in range(0,len(Y['EPIC']),1):
            if S['EPIC'][i] == Y['EPIC'][k]:
                S['Y'][i] = 1
                # S['YDnu'][i] = df['dnu'][k]
                # S['e_YDnu'][i] = df['dnu_err'][k]
                # S['Ynumax'][i] = df['nmx'][k]
                # S['e_Ynumax'][i] = df['nmx_err'][k]

    Y['Nseismo'] = Y['Y'] + Y['B'] + Y['S']
    B['Nseismo'] = B['Y'] + B['B'] + B['S']
    S['Nseismo'] = S['Y'] + S['B'] + S['S']


    df = pd.merge(df,Y,how='inner',on=['EPIC'])
    df = df.reset_index(drop=True)
    df1 = pd.merge(df1,B,how='inner',on=['EPIC'])
    df1 = df1.reset_index(drop=True)
    df2 = pd.merge(df2,S,how='inner',on=['EPIC'])
    df2 = df2.reset_index(drop=True)

    # All = pd.concat([Y,B,S],ignore_index=True)
    # All = All.drop_duplicates(subset=['EPIC'])
    # All = All.reset_index(drop=True)

    # Yo = All[All['Y'] == 1.0]
    # Yo = Yo[Yo['B'] == 0.0]
    # Yo = Yo[Yo['S'] == 0.0]
    # Yno = All[All['Y'] == 0.0]
    # Bo = All[All['B'] == 0.0]
    # So = All[All['S'] == 0.0]

    # Yo.to_csv('/home/bmr135/Dropbox/K2Poles/Data0405/PARAM_out_MDs/Poles/YvonneOnly_seismic_'+pipeline,index=False)
    # Yno.to_csv('/home/bmr135/Dropbox/K2Poles/Data0405/PARAM_out_MDs/Poles/YvonneNoDet_seismic_'+pipeline,index=False)
    # B.to_csv('/home/bmr135/Dropbox/K2Poles/Data0405/PARAM_out_MDs/Poles/Benoit_nseismo_'+pipeline,index=False)
    # S.to_csv('/home/bmr135/Dropbox/K2Poles/Data0405/PARAM_out_MDs/Poles/Savita_nseismo_'+pipeline,index=False)
    # All.to_csv('/home/bmr135/Dropbox/K2Poles/Data0405/PARAM_out_MDs/Poles/All_nseismo_'+pipeline,index=False)
    print( pipeline)

    return df, df1, df2

def nyq_refl(x,dnu,t,name):
    ''' Check if power is potentially reflected in the Nyquist limit or not
        x = input list of dataframes
        dnu = labels for relevant dnu values
        t = cadence in minutes
        name = name of pipeline/campaign '''

    nyq = 1e6/(2*t*60)
    dnu_nyq = 0.263*nyq**0.77 # Stello 2009 dnu/numax relation
    print( dnu_nyq)
    j=0
    for i in x:
        df = i
        # print( df.columns.values)
        df['Nyq'] = 0.0
        for k in range(0,len(df['EPIC']),1):
            if df[dnu[j]][k] > dnu_nyq:
                df['Nyq'][k] = 1.0
                print( df['EPIC'][k], " is Super-Nyquist, ", name[j], ", Nseismo = ", df['Nseismo'][k])
        i = df
        j+=1

    print( "Nyquist Checked")
    return x

def single_seismo(df,keys,output):
    ''' Generation of single column of seismic values '''
    a,b,c = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
    a['EPIC'],b['EPIC'],c['EPIC'] = df['EPIC'],df['EPIC'],df['EPIC']
    a[output] = df[keys[0]]
    a = a[a[output] != 'NaN']
    b[output] = df[keys[1]].dropna(axis=0,how='any')
    b = b[b[output] != 'NaN']
    c[output] = df[keys[2]].dropna(axis=0,how='any')
    c = c[c[output] != 'NaN']
    df1 = pd.concat([a,b,c],ignore_index=True)
    df1 = df1.drop_duplicates(subset=['EPIC'])
    df1 = df1.reset_index(drop=True)
    if len(df) == len(df1):
        df = pd.merge(df,df1,how='inner',on=['EPIC'])
        df = df.reset_index(drop=True)
    return df

def met_filter(df):
    ''' Filter out [Fe/H] values from EPIC that shouldn't be trusted '''
    for i in range(len(df)):
        if (df['stpropflag'][i] == 'rpm') or (df['stpropflag'][i] == 'col') or (df['stpropflag'][i] == 'plx'):
            df['[Fe/H]'][i] = -99.9
            df['em_[Fe/H]'][i] = -99.9
            df['ep_[Fe/H]'][i] = -99.9
    # print(df['[Fe/H]'])
    # df['[Fe/H]'] = df.where(df['stpropflag'] == 'rpm', -99.9)
    # print(df['[Fe/H]'])

    return df
