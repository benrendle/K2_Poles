''' Functions specifically for use in param_out.py '''

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import sys
from pandas import DataFrame, read_csv
import scipy.odr.odrpack as odrpack
import K2_constants as const

def param_in(lists,name):
    ''' Reading in of PARAM processed files '''
    j=0
    for i in range(0,len(lists),1):
        lists[i] = pd.read_csv('/home/bmr135/GA/K2Poles/param_outputs/'+name[j]+'.in.me',delimiter=r'\s+')
        # lists[i] = lists[i][lists[i].dist > -10] # Read in a cut BC3 param output file
        try: lists[i].rename(columns={'#Id':'EPIC'},inplace=True)
        except: pass
        # lists[i].reset_index(inplace=True)
        j+=1
    return lists

def param_inputs(lists,name):
    ''' Reading in of files prior to PARAM run '''
    j=0
    for i in range(0,len(lists),1):
        if i < 7:
            lists[i] = pd.read_csv('/home/bmr135/GA/K2Poles/matlab_in/'+name[j]+'.csv',delimiter=r'\s+')
        if i >= 7:
            lists[i] = pd.read_csv('/home/bmr135/GA/K2Poles/matlab_in/'+name[j]+'.csv')
        # lists[i] = lists[i][lists[i].dist > -10] # Read in a cut BC3 param output file
        # lists[i].reset_index(inplace=True)
        j+=1
    return lists

def param_inputs_csv(lists,name):
    ''' Reading in of files prior to PARAM run '''
    j=0
    for i in range(0,len(lists),1):
        lists[i] = pd.read_csv('/home/bmr135/GA/K2Poles/matlab_in/'+name[j]+'.csv')
        j+=1
    return lists

def selection_function(data,numax):
    ''' Selection function used for K2 GAP and unbiasing the
        asteroseismic dectections. Based on proposals by Dennis Stello for
        the C3 and C6 campaigns '''
    for i in range(0,len(data),1):
        df = data[i]
        if i < 1: # C3
            df = df[ (df[numax[i]] > 10) & (df[numax[i]] < 280) & (df['Hmag'] > 7) & (df['Hmag'] < 12) & (df['JK'] > 0.5) & (df['prob_s'] >= 0.95) ]
            # df = df[ (df['Hmag'] > 7) & (df['Hmag'] < 12) & (df['JK'] > 0.5) ]
        elif i >= 1: # C6
            df = df[ (df[numax[i]] > 10) & (df[numax[i]] < 280) & (df['Vcut'] > 9) & (df['Vcut'] < 15) & (df['JK'] > 0.5) & (df['prob_s'] >= 0.95) ]
            # df = df[ (df['Vcut'] > 9) & (df['Vcut'] < 15) & (df['JK'] > 0.5) ]

        data[i] = df

    return data

def uncerts(param,uncert,pd):
    ''' Calculate the uncertainties on each PARAM output
        using the difference between the value and 68th
        percentiles '''
    # print pd
    for k in range(0,len(pd),1):
        X = pd[k]
        j=0
        for i in param:
            # print i
            # print X[i]
            X[uncert[j]] = X[i] - X[i+'_68L']
            X[uncert[j+1]] = X[i+'_68U'] - X[i]
            j+=2
    return pd

def least_squares_2err(Matched):
    ''' Least squares fitting routine '''
    def f(Par,x):
        return Par[0]*x + Par[1]
    # runs=200
    # samp_size=len(Matched['mass_x'])/5#int(np.sqrt(len(x)))
    Matched['xerr'] = (Matched['e_mass1_x'] + Matched['e_mass2_x'])/2.0
    Matched['yerr'] = (Matched['e_mass1_y'] + Matched['e_mass2_y'])/2.0
    mpar,cpar=[],[] # y = mx + c
    empar,ecpar=[],[]
    # for i in range(0,runs,1):
    #     # idx=np.random.choice(len(x),samp_size,replace=False)
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(Matched['mass_x'], Matched['mass_y'], sx=Matched['xerr'], sy=Matched['yerr'])#, sy=Matched['yerr'])
    myodr = odrpack.ODR(mydata, linear, beta0=[1., 0.],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    myoutput.pprint()
    # print mydata
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
    # print mydata
    return mpar, cpar, empar, ecpar
