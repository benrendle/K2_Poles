''' Functions file for mass_distr.py '''

import pandas as pd
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import scipy.odr.odrpack as odrpack
from scipy import integrate
from scipy.stats import beta
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import seaborn as sns
import colormaps
# from path import path
from pyqt_fit import kde
import sim_pert as sp
from decimal import Decimal
import pystan
from pystan import stan

def p_in(name,df1,field):
    ''' Extract input data from input PARAM files as necessary.

        epic = df containing EPIC numbers
        name = name of file
        df1 = input df containing data
    '''
    EPIC = pd.DataFrame()
    EPIC['#Id'] = df1['#Id']
    df = pd.read_csv('/home/bmr135/GA/K2Poles/param_inputs/Poles/'+name+'.in',delimiter=r'\s+')
    # df = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/param_inputs/Poles/'+name+'.in',delimiter=r'\s+')
    df.rename(columns={'ID':'#Id'},inplace=True)
    df = pd.merge(df,EPIC,how='inner',on=['#Id'])
    df.reset_index(drop=True)
    df1['nmx'] = df['numax']
    df1['nmx_err'] = df['enumax']
    df1['dnu'] = df['Dnu']
    df1['dnu_err'] = df['eDnu']
    df1['Teff'] = df['teff']
    df1['Teff_err'] = df['eteff']
    df1['Glat'] = df['GLAT']
    df1['Glon'] = df['GLON']
    df1['feh'] = df['feh']
    df1['feh_err'] = df['efeh']
    df1['alpha'] = df['alpha']
    df1['m_seismo'] = (df1['nmx']/3090)**3 * (df1['dnu']/135.1)**-4 * (df1['Teff']/5777)**1.5
    df1 = df1[df1['dist'] > -90.0]
    df1 = df1[df1.mass > 0.0]
    df1 = df1[df1.feh > -90.0]
    df1 = df1[df1.age > 0.0]
    df1['logAge'] = np.log10(df1['age']*10**9)
    # if name[:2] != 'TC':
    #     df1 = Nseismo_link(field,df1)
    # df1 = df1[df1.Nseismo == 1.0] # Filter for a specific number of seismic pipelines with detections

    return df1

def Nseismo_link(field,df):
    ''' Match number of seismic detections to EPIC numbers '''
    a = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/PARAM_out_MDs/Poles/All_nseismo_K2P2_'+field)
    # a = pd.read_csv('/home/ben/Dropbox/K2Poles/Data0405/PARAM_out_MDs/Poles/All_nseismo_K2P2_'+field)
    a.rename(columns={'EPIC':'#Id'},inplace=True)
    df = pd.merge(df,a,how='inner',on=['#Id'])
    df.reset_index(drop=True)

    return df

def MZ_DistDiff(spec,MAST):
    ''' Mass vs vertical height plot using SR masses (const. mass) '''
    spec['Z'] = spec['dist']*np.sin(spec['Glat']*np.pi/180)*1e-3
    spec['Z_err1'] = (spec['dist']-spec['dist_68L'])*(spec['Z']/spec['dist'])
    spec['Z_err2'] = (spec['dist_68U']-spec['dist'])*(spec['Z']/spec['dist'])
    spec['mass_err1'] = spec['mass']-spec['mass_68L']
    spec['mass_err2'] = spec['mass_68U']-spec['mass']
    MAST['Z'] = MAST['dist']*np.sin(MAST['Glat']*np.pi/180)*1e-3
    MAST['Z_err1'] = (MAST['dist']-MAST['dist_68L'])*(MAST['Z']/MAST['dist'])
    MAST['Z_err2'] = (MAST['dist_68U']-MAST['dist'])*(MAST['Z']/MAST['dist'])
    MAST['mass_err1'] = MAST['mass']-MAST['mass_68L']
    MAST['mass_err2'] = MAST['mass_68U']-MAST['mass']
    plt.figure()
    plt.errorbar(spec['m_seismo'],spec['Z'],yerr=spec['Z_err1'],alpha=0.7,label=r'M$_{ws}$',fmt='o')
    # plt.scatter(spec['mass'],spec['Z'],alpha=0.7,label=r'M$_{ws}$',color='r')
    plt.errorbar(MAST['m_seismo'],MAST['Z'],yerr=MAST['Z_err1'],alpha=0.7,label=r'M$_{wos}$',color='r',fmt='o')
    plt.xlabel(r'Mass, [$M_{\odot}$]')
    plt.ylabel(r'Z [kpc]')
    plt.legend()
    plt.tight_layout()

def MZ(spec1,spec2,MAST,APK):
    ''' Mass vs vertical height plot '''
    spec1['Z'] = spec1['dist']*np.sin(spec1['Glat']*np.pi/180)*1e-3
    spec1['Z_err1'] = (spec1['dist']-spec1['dist_68L'])*(spec1['Z']/spec1['dist'])
    spec1['Z_err2'] = (spec1['dist_68U']-spec1['dist'])*(spec1['Z']/spec1['dist'])
    spec1['mass_err1'] = spec1['mass']-spec1['mass_68L']
    spec1['mass_err2'] = spec1['mass_68U']-spec1['mass']

    spec2['Z'] = spec2['dist']*np.sin(spec2['Glat']*np.pi/180)*1e-3
    spec2['Z_err1'] = (spec2['dist']-spec2['dist_68L'])*(spec2['Z']/spec2['dist'])
    spec2['Z_err2'] = (spec2['dist_68U']-spec2['dist'])*(spec2['Z']/spec2['dist'])
    spec2['mass_err1'] = spec2['mass']-spec2['mass_68L']
    spec2['mass_err2'] = spec2['mass_68U']-spec2['mass']

    MAST['Z'] = MAST['dist']*np.sin(MAST['Glat']*np.pi/180)*1e-3
    MAST['Z_err1'] = (MAST['dist']-MAST['dist_68L'])*(MAST['Z']/MAST['dist'])
    MAST['Z_err2'] = (MAST['dist_68U']-MAST['dist'])*(MAST['Z']/MAST['dist'])
    MAST['mass_err1'] = MAST['mass']-MAST['mass_68L']
    MAST['mass_err2'] = MAST['mass_68U']-MAST['mass']

    # spec1 = spec1[spec1.Z > 0] # Toggle for C3/C6
    # spec2 = spec2[spec2.Z < 0]
    # spec1 = spec1[spec1.alpha > -5.0]
    # spec2 = spec2[spec2.alpha > -5.0]
    APK = APK[APK['mass'] > 0]

    plt.figure()
    plt.scatter(APK['mass'],APK['Z']*1e-3,alpha=0.2,label=r'Kepler',color='grey')
    plt.scatter(spec1['mass'],spec1['Z'],alpha=0.4,label=r'K2 C6',color='b')
    plt.scatter(spec2['mass'],spec2['Z'],alpha=0.4,label=r'K2 C3',color='orange')
    ''' Errorbar plots '''
    # plt.errorbar(spec1['mass'],spec1['Z'],xerr=[spec1['mass_err1'],spec1['mass_err2']], \
    #             yerr=[spec1['Z_err1'],spec1['Z_err2']],alpha=0.2,label=r'RC6$',fmt='.')
    # plt.errorbar(spec2['mass'],spec2['Z'],xerr=[spec2['mass_err1'],spec2['mass_err2']], \
    #             yerr=[spec2['Z_err1'],spec2['Z_err2']],alpha=0.7,label=r'APO$_{wos}$',color='r',fmt='o')
    # plt.scatter(spec1['mass'],spec1['Z'],c=spec1['age'],cmap=colormaps.parula,s=30,label='_nolegend_')
    # cbar = plt.colorbar()
    # cbar.set_label(r'Age [Gyr]', rotation=270, fontsize=20, labelpad=25)
    # cbar.ax.tick_params(labelsize=20)
    plt.xlabel(r'Mass, [$M_{\odot}$]',fontsize=20)
    plt.ylabel(r'Z [kpc]',fontsize=20)
    plt.tick_params(labelsize=15)
    plt.ylim(-10,10)
    plt.legend()
    plt.tight_layout()

def spectro(spec,K2):
    ''' Deals with spectroscopic data, comparing EPICS with and without
        spectroscopic data included in the PARAM input file '''
    EPIC = pd.DataFrame()
    EPIC['#Id'] = spec['#Id']
    K2 = pd.merge(K2,EPIC,how='inner',on=['#Id'])
    K2.reset_index(drop=True)
    ID = pd.DataFrame()
    ID['#Id'] = K2['#Id']
    spec = pd.merge(spec,ID,how='inner',on=['#Id'])
    # spec.rename(columns={'#Id':'EPIC'},inplace=True)
    spec.reset_index(drop=True)
    K2['Glat'] = spec['Glat']
    # spec['m_seismo'] = (spec['nmx']/3090)**3 * (spec['dnu']/135.5)**-4 * (spec['Teff']/5777)**1.5
    # spec.to_csv('/home/bmr135/GA/K2Poles/Gaia_ESO/GES_C3_ws.csv',index=False)
    # K2.to_csv('/home/bmr135/GA/K2Poles/Gaia_ESO/C3_ws.csv',index=False)

    return spec, K2

def pop_dist(df1,colour,label,scat):
    ''' C3/C6 distributions on a single plot '''
    plt.figure()
    j=0
    for i in df1:
        df = i
        if label[j] == 'APOKASC':
            plt.scatter(df['mass'],df['Z']*1e-3,alpha=0.2,label=label[j],color=colour[j])
            j+=1
        else:
            df['Z'] = df['dist']*np.sin(df['Glat']*np.pi/180)*1e-3
            df['Z_err1'] = (df['dist']-df['dist_68L'])*(df['Z']/df['dist'])
            df['Z_err2'] = (df['dist_68U']-df['dist'])*(df['Z']/df['dist'])
            df['mass_err1'] = df['mass']-df['mass_68L']
            df['mass_err2'] = df['mass_68U']-df['mass']
            plt.scatter(df['mass'],df['Z'],alpha=0.7,label=label[j],color=colour[j])
            j+=1
    plt.xlabel(r'Mass, [$M_{\odot}$]',fontsize=20)
    plt.ylabel(r'Z [kpc]',fontsize=20)
    plt.tick_params(labelsize=15)
    plt.legend()
    plt.tight_layout()

def c3_6_param(df,param,label):
    ''' Scatter plots with colourbars for combined C3/C6 data '''
    # df = pd.concat([df3,df6],ignore_index=True)
    # df.reset_index(drop=True)
    df['Z'] = df['dist']*np.sin(df['Glat']*np.pi/180)*1e-3
    df = df[df.alpha > -5]
    plt.figure()
    plt.scatter(df['mass'],df['Z'],c=df[param],s=20,cmap=colormaps.parula)
    cbar = plt.colorbar()
    cbar.set_label(label, rotation=270, fontsize=20, labelpad=25)
    cbar.ax.tick_params(labelsize=20)
    plt.xlabel(r'Mass, [$M_{\odot}$]',fontsize=20)
    plt.ylabel(r'Z [kpc]',fontsize=20)
    plt.tick_params(labelsize=15)
    # plt.legend()
    plt.tight_layout()

def met_comp(X,name,p1,p2,colour,label):
    ''' Alpha/Fe Fe/H plot or comparison of spectro comp and age '''
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
    j=0
    for i in X:
        df = i
        df = df[df[p2] > -5]
        axScatter.scatter(df[p1],df[p2],alpha=0.3,label=name[j],color=colour[j])
        # plt.scatter(df[p1],df[p2],alpha=0.3,label=name[j],color=colour[j])

        axScatter.set_xlabel(label[0],fontsize=20)
        axScatter.set_ylabel(label[1],fontsize=20)
        axScatter.tick_params(labelsize=15)

        # plt.xlabel(label[0])
        # plt.ylabel(label[1])

        # bin_list = np.arange(-3, 1.5, 0.1)
        # axFEH.hist(df[p1][np.isfinite(df[p1])], bins=bin_list, normed=1, facecolor=colour[j], label=name[j], histtype='step')
        # axALP.hist(df[p2][np.isfinite(df[p2])], bins=bin_list, normed=1, facecolor=colour[j], label=name[j],orientation='horizontal', histtype='step')

        densityFe = kde.KDE1D(df[p1])
        densityAl = kde.KDE1D(df[p2])
        x1 = np.r_[min(df[p1]):max(df[p1]):1024j]
        x2 = np.r_[min(df[p2]):max(df[p2]):1024j]
        axFEH.plot(x1,densityFe(x1),color=colour[j],linewidth=2)
        axALP.plot(densityAl(x2),x2,color=colour[j],linewidth=2)

        # axFEH.set_xlim(axScatter.get_xlim())
        # axALP.set_ylim(axScatter.get_ylim())
        axScatter.legend(loc=3,prop={'size':15})
        # plt.legend(loc=3,prop={'size':15})

        j+=1

def histo(df,name,bins,label,log,tag):
    ''' Distribution histogram '''
    # plt.figure()
    n, bins, patches = plt.hist(df[name][np.isfinite(df[name])], bins=bins, label=tag, linewidth=2, histtype='step', normed=True)
    x=np.zeros(len(n))
    for i in np.arange(0,len(bins),1):
        try: x[i] = (bins[i] + bins[i+1])/2
        except: pass
    err = 1.0/np.sqrt(n)
    # plt.errorbar(x,n,yerr=[err,err],fmt='.')
    plt.ticklabel_format(useOffset=False)	# Should turn off the scientific label formatting (+n*10^x notations)
    plt.xlabel(label,fontsize=20)
    # plt.xlim(0.7,20)
    plt.tick_params(labelsize=15)
    if log == 1:
        plt.xscale('log')
    plt.tight_layout()

def least_squares(df1,df2):
    ''' Least squares fitting approach to calculating offsets betrween Teffs '''
    def f(Par,x):
        return Par[0]*x + Par[1]
    df2['alt_feh'] = df2['feh'] - df1['feh']
    df1['null'] = 0
    mpar, cpar, empar, ecpar = [], [], [], []
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(df1['null'], df2['alt_feh'])#, sy=df2['feh_err'])
    myodr = odrpack.ODR(mydata, linear, beta0=[0., 0.],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    # myoutput.pprint()

    return mpar, cpar, empar, ecpar

def least_squares2(df1,df2):
    ''' Least squares fitting approach to calculating offsets between Teffs '''
    def f(Par,x):
        return Par[0]*x + Par[1]

    df2['alt_Teff'] = df2['Teff'] - df1['Teff']
    df1['null'] = 0
    mpar, cpar, empar, ecpar = [], [], [], []
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(df1['null'], df2['alt_Teff'])#, sx=df1['Teff_err'], sy=df2['Teff_err'])
    myodr = odrpack.ODR(mydata, linear, beta0=[0., 0.],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    # myoutput.pprint()

    return mpar, cpar, empar, ecpar

def uncerts(df,param,x):
    ''' Provide single percentage uncertainty from PARAM values
        x = uncertainty placeholder '''
    j=0
    for i in param:
        if i == 'logAge':
            d1 = df['age']-df['age_68L']
            d2 = df['age_68U'] - df['age']
            df['age_err'] = np.sqrt(d1**2 + d2**2)/2
            df[i+'_err'] = 100 * (df['age_err']/(df['age']*np.log(10))) / df[i]
        else:
            d1 = df[i]-df[i+'_68L']
            d2 = df[i+'_68U'] - df[i]
            df[i+'_err'] = 100*(np.sqrt(d1**2 + d2**2)/2)/df[i]
        x[j][1] = np.mean(df[i+'_err'])
        j+=1
    return x

def sim_pert(obs,sim):
    ''' Generate uncertainties for perturbing the simulations '''
    obs = obs[obs['feh'] > -5]
    param = ['mass','rad','age','logAge']
    x = [['Mass',0.0],['radius',0.0],['age',0.0],['logAge',0.0]]
    err = uncerts(obs,param,x)
    efeh = 100*(0.117/np.average(obs['feh']))
    eTeff = 100*(117/np.average(obs['Teff']))
    err.append(['M_H',abs(efeh)])
    err.append(['Teff',eTeff])
    pert = sp.PERTURBATION(sim,err)
    alt_sim = pert()

    return alt_sim

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def vert_dist(i):
    ''' Calculate true vertical distance from coordinates and PARAM output '''
    d = coord.Galactic(l=i['Glon'].as_matrix()*u.degree, b=i['Glat'].as_matrix()*u.degree, distance=i['dist'].as_matrix()*u.pc)
    f = d.transform_to(coord.Galactocentric)
    i['Z'] = (f.z*u.pc)*1e-3
    return i

def spec(i):
    ''' Specification of code to run depending upon input '''
    itrs=0
    samp=0
    clump=0
    # With Clump, no sampling
    if i == '0': itrs = 1; samp = 0; clump = 0; folder_loc = 'Clump_No_Samp/'
    # Without Clump, no sampling
    if i == '1': itrs = 1; samp = 0; clump = 1; folder_loc = 'No_Clump_No_Samp/'
    # With Clump, sampled
    if i == '2': itrs = 1; samp = 1; clump = 0; folder_loc = 'Clump_Samp/'
    # Without Clump, sampled
    if i == '3': itrs = 1; samp = 1; clump = 1; folder_loc = 'No_Clump_Samp/'

    return itrs, samp, clump, folder_loc

def Z_subplots(df1,df2,param,ran,xtag,z):
    ''' Plots of distributions of mass and age at different heights above the plane '''
    # z = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8]
    n = len(z) - 1
    m = int(n/2)
    fig, axes = plt.subplots(m,2,sharex=True,sharey=True)
    if len(z) > 9:
        ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13 = axes.flatten()
        a = [ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13]
    if len(z) == 9:
        ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7 = axes.flatten()
        a = [ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7]
    for i in range(n):
        zf1 = df1[(df1['Z'] > z[i]) & (df1['Z'] < z[i+1])]
        a[i].hist(zf1[param],bins=ran)
        # a[i].annotate(r'%s $<$ Z $<$ %s' %(z[i],z[i+1]), xy=(0.55, 0.8), xycoords='axes fraction')
        a[i].set_title(r'%s $<$ Z $<$ %s' %(z[i],z[i+1]))
        # print(min(zf1['Z']),max(zf1['Z']))
    a[len(z)-3].set_xlabel(xtag,fontsize=15)
    a[len(z)-2].set_xlabel(xtag,fontsize=15)
    # plt.suptitle(r'K2 '+param+r' Distribution as Function of Z [kpc]',fontsize=15)
    plt.tight_layout()

def sim_coords(df,field):
    ''' Give simulations values for GLAT and GLON from the input coords '''
    if field == 3:
        coords = pd.read_csv('/home/bmr135/GA/K2Poles/avk_k2_complete_c3.txt',delimiter=r'\s+')
        # coords = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/avk_k2_complete_c3.txt',delimiter=r'\s+')
    if field == 6:
        coords = pd.read_csv('/home/bmr135/GA/K2Poles/avk_k2_complete_c6.txt',delimiter=r'\s+')
        # coords = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/avk_k2_complete_c6.txt',delimiter=r'\s+')
    if field == 'K':
        coords = pd.read_csv('/home/bmr135/GA/K2Poles/avk_kep_complete.txt',delimiter=r'\s+')
        # coords = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/avk_kep_complete.txt',delimiter=r'\s+')

    df['Glon'] = 0
    df['Glat'] = 0
    df2 = pd.DataFrame([])
    for i in range(int(len(df)/len(coords))):
        df2 = df2.append(coords.sample(n=len(coords)),ignore_index=True)
    x = len(df)%len(coords)
    df2 = df2.append(coords.sample(n=x))
    df2 = df2.reset_index(drop=True)
    df2['y'] = np.random.uniform(-0.5,0.5,len(df))
    df['Glon'] = df2['l'] + df2['y']
    df['Glat'] = df2['b'] + df2['y']

def sim_dist(df):
    ''' Distances to stars using the methods laid out in Pijpiers 2014 and
        Miglio 2013 '''
    a = -0.190537e5
    b = 0.155145e5
    c = -0.421279e4
    d = 0.381476e3
    Mbol_s = 4.72 # Torres 2010, recommended solar bolometric mag. value
    # df['BCv'] = a + b*df['logTe'] + c*(df['logTe'])**2 + d*(df['logTe'])**3
    # ''' Pijpers 2014 '''
    # df['pi'] = 10**(0.5*(4 + 0.4*Mbol_s - 0.4*(df['Vcut'] - df['Av'] + df['BCv']) - df['logL']))
    # df['dist'] = 1/df['pi']
    ''' Miglio 2013 '''
    df['dist'] = 10**(1 + 2.5*(df['logTe'] - np.log10(5777)) + np.log10(df['numax']/3090) \
                    - 2*np.log10(df['deltanu']/135.1) + 0.2*(df['mbolmag'] - Mbol_s))
    # print(min(df['dist']),max(df['dist']))

def space_density(obs):
    ''' Calculate the space density for a given band above the Galactic plane '''
    # obs['Zabs'] = np.abs(obs['Z'])
    a = min(obs['Z'])
    b = max(obs['Z'])
    theta = 0.5*np.sqrt(116)
    t_theta = np.tan(theta*(np.pi/180))**2
    bins = np.linspace(np.log10(100),np.log10(5000),19)
    Z = np.log10(np.abs(obs['Z'])*1000)
    hist, bin_edges = np.histogram(Z, bins = bins)
    # print(bin_edges)
    # print(hist)
    volume = []
    # volume of square based pyramid section
    for i in range(len(bin_edges)-1):
        V1 = (4*t_theta*(10**bin_edges[i])**3)/3
        V2 = (4*t_theta*(10**bin_edges[i+1])**3)/3
        Vol = V2 - V1
        volume.append(Vol)
    # volume of cube
    # for i in range(len(bin_edges)-1):
    #     V1 = (2*np.sqrt(t_theta)*10**max(bin_edges))**2 * (10**(bin_edges[i+1]) - 10**(bin_edges[i]))
    #     volume.append(V1)
    rho = []
    rho10 = np.zeros(len(volume))
    for i in range(len(volume)):
        density = hist[i]/volume[i]
        rho10[i] = density
        rho.append(np.log10(density)) # In units number of stars/pc^3

    print(rho)
    bin_10 = 10**bin_edges[1:]
    print( len(bin_10), len(rho10))
    ''' Least squares fitting '''
    def f(Par,z):
        return Par[0]*np.exp(-z/Par[1])
    mpar, cpar, empar, ecpar = [], [], [], []
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(bin_10[5:12], rho10[5:12])#, sy=df2['feh_err'])
    myodr = odrpack.ODR(mydata, linear, beta0=[1e-3, 500],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    # myoutput.pprint()

    b = bin_10[14:]
    b = np.delete(b,[1])
    r = rho10[14:]
    r = np.delete(r,[1])
    def f1(Par,z):
        return Par[0]*np.exp(-z/Par[1])
    mpar1, cpar1, empar1, ecpar1 = [], [], [], []
    linear1 = odrpack.Model(f1)
    mydata1 = odrpack.RealData(b, r)#, sy=df2['feh_err'])
    myodr1 = odrpack.ODR(mydata1, linear1, beta0=[5e-4, 1000],maxit=20000)
    myoutput1 = myodr1.run()
    mpar1.append(myoutput1.beta[0])
    cpar1.append(myoutput1.beta[1])
    empar1.append(myoutput1.sd_beta[0])
    ecpar1.append(myoutput1.sd_beta[1])
    # myoutput1.pprint()

    # def f2(Par,z):
    #     return Par[0]*np.exp(-z/Par[1])
    # mpar2, cpar2, empar2, ecpar2 = [], [], [], []
    # linear2 = odrpack.Model(f2)
    # mydata2 = odrpack.RealData(bin_10[4:], rho10[4:])#, sy=df2['feh_err'])
    # myodr2 = odrpack.ODR(mydata2, linear2, beta0=[1e-4, 600],maxit=20000)
    # myoutput2 = myodr2.run()
    # mpar2.append(myoutput2.beta[0])
    # cpar2.append(myoutput2.beta[1])
    # empar2.append(myoutput2.sd_beta[0])
    # ecpar2.append(myoutput2.sd_beta[1])
    # myoutput2.pprint()

    plt.figure()
    plt.scatter(bin_10,rho)
    plt.scatter(bin_10[5:12],rho[5:12],color='r')
    plt.scatter(bin_10[13:],rho[13:],color='y')
    # plt.scatter(b,np.log10(r))
    plt.plot(bin_10,np.log10(myoutput.beta[0]*np.exp(-bin_10/myoutput.beta[1])),color='orange', \
            label=r'$\rho_{0} =$ %.4g pc$^{-3}$, H$_{\rm{z}} = $ %.4g pc'%(mpar[0],cpar[0]))
    plt.plot(bin_10,np.log10(myoutput1.beta[0]*np.exp(-bin_10/myoutput1.beta[1])),color='orange',linestyle='--', \
            label=r'$\rho_{0} =$ %.4g pc$^{-3}$, H$_{\rm{z}} = $ %.4g pc'%(mpar1[0],cpar1[0]))
    plt.ylim(-7.75,-3.75)
    plt.xlim(100, 5100)
    plt.tick_params(labelsize=15)
    plt.ylabel(r'Log Space Density',fontsize=20)
    plt.xlabel(r'Z [pc]',fontsize=20)
    plt.legend(prop={'size':10})
    plt.tight_layout()
    # plt.savefig('/home/bmr135/Dropbox/K2Poles/pop_trends/211117/C3_Gilmore_Reid.png')
    # plt.show()

def space_density2(obs):
    ''' Improved calculation of the volume for determining the space density '''
    theta = np.deg2rad(np.sqrt(116))
    t1 = theta/2.0
    alpha_c3 = np.deg2rad(61.4) # degrees
    alpha_c6 = np.deg2rad(50.4) # degrees

    bins = np.linspace(np.log10(100),np.log10(5000),19)
    Z = np.log10(np.abs(obs['Z'])*1000)
    hist, bin_edges = np.histogram(Z, bins = bins)
    # print(bin_edges)
    print(hist)
    volume = []

    for i in range(len(bin_edges)-1):
        x2 = lambda x: ((10**bin_edges[i])**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - x)**-3)
        x3 = lambda xa: ((10**bin_edges[i+1])**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - xa)**-3)

        vol1 = integrate.quad(x2,0,theta)
        vol2 = integrate.quad(x3,0,theta)
        Vol = vol2[0] - vol1[0]
        volume.append(Vol)

    # print(volume)
    rho = []
    rho10 = np.zeros(len(volume))
    sigma = np.zeros(len(volume))
    for i in range(len(volume)):
        density = hist[i]/volume[i]
        rho10[i] = density
        sigma[i] = np.sqrt(hist[i])
        rho.append(np.log10(density)) # In units number of stars/pc^3

    print(sigma)
    # print(rho10)

    bin_10 = 10**bin_edges[1:]
    # x = pd.DataFrame()
    # x['Z'] = bin_10
    # x['log_rho'] = np.log10(rho10)
    # x['log_sig_rho'] = np.log10(np.sqrt(hist)/volume)
    # x['rho'] = 10**np.log10(rho10)
    # x['sig_rho'] = 10**np.log10(np.sqrt(hist)/volume)
    # print(x)
    # x.to_csv('/home/bmr135/Dropbox/K2Poles/pop_trends/Ioana_Scale_Heights/C6',index=False)
    # x.to_csv('/home/ben/Dropbox/K2Poles/pop_trends/Ioana_Scale_Heights/C6',index=False)
    sig_rho = np.sqrt(hist)/volume
    print(sig_rho)

    ''' Least squares fitting '''
    def f(Par,z):
        return Par[0]*np.exp(-z/Par[1])
    mpar, cpar, empar, ecpar = [], [], [], []
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(bin_10[6:11], rho10[6:11], sy=sig_rho[6:11])
    myodr = odrpack.ODR(mydata, linear, beta0=[1e-3, 500],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    # myoutput.pprint()

    b = bin_10[12:]
    # b = np.delete(b,[1])
    r = rho10[12:]
    # r = np.delete(r,[1])
    sig = sig_rho[12:]
    # sig = np.delete(sig,[1])
    def f1(Par,z):
        return Par[0]*np.exp(-z/Par[1])
    mpar1, cpar1, empar1, ecpar1 = [], [], [], []
    linear1 = odrpack.Model(f1)
    mydata1 = odrpack.RealData(b, r, sy=sig)
    myodr1 = odrpack.ODR(mydata1, linear1, beta0=[5e-4, 1000],maxit=20000)
    myoutput1 = myodr1.run()
    mpar1.append(myoutput1.beta[0])
    cpar1.append(myoutput1.beta[1])
    empar1.append(myoutput1.sd_beta[0])
    ecpar1.append(myoutput1.sd_beta[1])
    # myoutput1.pprint()

    def f2(Par,z):
        return (beta.pdf(abs(z),Par[4],Par[5])) + (Par[0])*np.exp(-abs(z)/Par[1]) + (Par[2])*np.exp(-abs(z)/Par[3])
        # return Par[0]*((Par[1]*z + Par[2])+(Par[3]*z + Par[4]))
    par2, epar2 = [], []
    linear2 = odrpack.Model(f2)
    mydata2 = odrpack.RealData(bin_10[3:], rho10[3:], sy=sig_rho[3:])
    myodr2 = odrpack.ODR(mydata2, linear2, beta0=[4e-4, 500, 3e-5, 1000, 1., 10.],maxit=20000)
    myoutput2 = myodr2.run()
    par2.append(myoutput2.beta[0])
    par2.append(myoutput2.beta[1])
    par2.append(myoutput2.beta[2])
    par2.append(myoutput2.beta[3])
    par2.append(myoutput2.beta[4])
    par2.append(myoutput2.beta[5])
    epar2.append(myoutput2.sd_beta[0])
    epar2.append(myoutput2.sd_beta[1])
    epar2.append(myoutput2.sd_beta[2])
    epar2.append(myoutput2.sd_beta[3])
    epar2.append(myoutput2.sd_beta[4])
    epar2.append(myoutput2.sd_beta[5])
    print(par2)
    print(epar2)
    # myoutput2.pprint()

    sig_rho = (np.sqrt(hist)/volume) / (rho10 * np.log(10))
    print(sig_rho)
    plt.figure()
    plt.errorbar(bin_10,rho,yerr=sig_rho, fmt='o', color='k')
    # plt.scatter(bin_10[6:11],rho[6:11],color='r')
    # plt.scatter(bin_10[12:16],rho[12:16],color='y')
    # plt.scatter(b,np.log10(r))
    plt.plot(bin_10, \
    np.log10((beta.pdf(bin_10,myoutput2.beta[4],myoutput2.beta[5])) + (myoutput2.beta[0]*np.exp(-bin_10/myoutput2.beta[1]) + myoutput2.beta[2]*np.exp(-bin_10/myoutput2.beta[3]))), \
    label=r'\zeta(\rm{Z})')
    plt.plot(bin_10,np.log10((myoutput2.beta[0])*np.exp(-bin_10/myoutput2.beta[1])),color='orange', \
            label=r'H$_{\rm{z}} = $ %.4g pc'%(myoutput2.beta[1]))
    plt.plot(bin_10,np.log10((myoutput2.beta[2])*np.exp(-bin_10/myoutput2.beta[3])),color='orange',linestyle='--', \
            label=r'H$_{\rm{z}} = $ %.4g pc'%(myoutput2.beta[3]))
            # label=r'$\rho_{0} =$ %.4g pc$^{-3}$, H$_{\rm{z}} = $ %.4g pc'%(myoutput2.beta[3]))

    hz = (rho10)*np.exp(-abs(bin_10)/500) + ((rho10))*np.exp(-abs(bin_10)/1000)
    # plt.plot(bin_10,np.log10(hz))


    plt.ylim(-7.75,-3.0)
    plt.xlim(100, 5100)
    plt.tick_params(labelsize=15)
    plt.ylabel(r'Log Space Density',fontsize=20)
    plt.xlabel(r'Z [pc]',fontsize=20)
    plt.legend(prop={'size':10})
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.show()

def sp3(obs):
    '''
    Using stan to perform fitting
    '''
    theta = np.deg2rad(np.sqrt(116))
    t1 = theta/2.0
    alpha_c3 = np.deg2rad(61.4) # degrees
    alpha_c6 = np.deg2rad(50.4) # degrees

    bins = np.linspace(np.log10(100),np.log10(5000),19)
    Z = np.log10(np.abs(obs['Z'])*1000)
    hist, bin_edges = np.histogram(Z, bins = bins)
    theta = np.deg2rad(np.sqrt(116))
    t1 = theta/2.0
    alpha_c3 = np.deg2rad(61.4) # degrees
    alpha_c6 = np.deg2rad(50.4) # degrees

    bins = np.linspace(np.log10(100),np.log10(5000),19)
    Z = np.log10(np.abs(obs['Z'])*1000)
    hist, bin_edges = np.histogram(Z, bins = bins)
    volume = []
    for i in range(len(bin_edges)-1):
        x2 = lambda x: ((10**bin_edges[i])**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - x)**-3)
        x3 = lambda xa: ((10**bin_edges[i+1])**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - xa)**-3)
        vol1 = integrate.quad(x2,0,theta)
        vol2 = integrate.quad(x3,0,theta)
        Vol = vol2[0] - vol1[0]
        volume.append(Vol)

    rho = []
    rho10 = np.zeros(len(volume))
    for i in range(len(volume)):
        density = hist[i]/volume[i]
        rho10[i] = density
        rho.append(np.log10(density)) # In units number of stars/pc^3
    bin_10 = 10**bin_edges[1:]
    dens_model = '''
    data {
        int<lower=0> N;
        real M[N];
        real Z[N];
    }
    parameters {
        //real<lower =0.001> sigma;
        //real<lower = 0> alpha;
        real M_true_std[N];
        real<lower = 0> sh1;
        real<lower = 0> sh2;
        real<lower = 0> c1;
        real<lower = 0> c2;
    }
    //transformed parameters {
    //    real M_true[N]; // Transform from N(0,1) back to M
    //    for (i in 1:N)
    //        M_true[i] = mu + sigma * M_true_std[i];
    //}
    model {
        M_true_std ~ c1*np.exp(-abs(Z)/sh1) + c2*np.exp(-abs(Z)/sh2)
        sh1 ~ normal( 500, 200);
        sh2 ~ normal( 800, 200);
        c1 ~ normal( 0, 1);
        c2 ~ normal( 0, 1);
        //sigma ~ normal(0, 10);
    }
    '''
    sm = pystan.StanModel(model_code=dens_model, model_name='DensityModel')
    dat = {'N': len(rho10),
           'M': rho10,
           'Z': bin_10
          }
    fit = sm.sampling(data=dat, iter=5000, chains=1, pars=['sh1', 'sh2', 'c1', 'c2'])
    print(fit)
    fit.plot()

    plt.figure()
    plt.scatter(bin_10,rho)
    # plt.plot(bin_10,np.log10(myoutput2.beta[0]*np.exp(-bin_10/myoutput2.beta[1])),color='orange', \
    #         label=r'$\rho_{0} =$ %.4g pc$^{-3}$, H$_{\rm{z}} = $ %.4g pc'%(mpar2[0],cpar2[0]))

    plt.plot(bin_10, \
            np.log10(fit['c1'].mean()*np.exp(-bin_10/fit['sh1'].mean()) + fit['c2'].mean()*np.exp(-bin_10/fit['sh2'].mean())))
            #  - (*(np.tan(myoutput2.beta[5]*z))^-1 + )), \

    # plt.ylim(-7.75,-3.5)
    plt.xlim(100, 5100)
    plt.tick_params(labelsize=15)
    plt.ylabel(r'Log Space Density',fontsize=20)
    plt.xlabel(r'Z [pc]',fontsize=20)
    plt.legend(prop={'size':10})
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
