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
import scipy.odr.odrpack as odrpack
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
import TAR as tt

mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'

if __name__ == '__main__':

    ''' Set up of star/sampling type to be analysed from user input '''
    itrs, samp, clump, folder_loc = mdf.spec(sys.argv[1])
    print('Iterations = ',itrs)
    if clump == 0:
        print('Running with clump included')
    else: print('Running with clump cut (R < 9.0)')
    if samp == 0:
        print('Running without sampling')
    else: print('Running with sampling')

    # ext = '/home/bmr135/GA/K2Poles/param_outputs/Poles/'
    # ext_fig = '/home/bmr135/Dropbox/K2Poles/pop_trends/051017/'
    ext = '/media/ben/SAMSUNG/GA/K2Poles/param_outputs/Poles/'
    ext_fig = '/home/ben/Dropbox/K2Poles/pop_trends/051017/'

    print('Figures saved to: '+ ext_fig + folder_loc)

    save_out = 0 # 0 - no saving; 1 - save all figures

    ''' Read in files '''
    # APK = pd.read_csv('/home/bmr135/GA/K2Poles/APOKASC4BEN.txt')
    APK = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/APOKASC4BEN.txt')

    # a = [['APOKASC.tar.gz', 'APOKASC/', 'APOKASC.in.me'], \
    #      ['GES.tar.gz', 'GES/', 'GES.in.me'], \
    #      ['GES_v2.tar.gz', 'GES_v2/', 'GES.in.me'], \
    #      ['GES_ns.tar.gz', 'GES_ns/', 'GES.in.me'], \
    #      ['GES_A40.tar.gz', 'GES_A40/', 'GES.in.me'], \
    #      ['APO6.tar.gz', 'APO6/','APOGEE.in.me'], \
    #      ['APO6_ns.tar.gz', 'APO6_ns/', 'APOGEE.in.me'], \
    #      ['RC3.tar.gz', 'RC3/', 'RC3.in.me'], \
    #      ['RC3_ns.tar.gz', 'RC3_ns/', 'RC3.in.me'], \
    #      ['RC3_A40.tar.gz', 'RC3_A40/', 'RC3.in.me'], \
    #      ['RC6.tar.gz', 'RC6/', 'RC6.in.me'], \
    #      ['RC6_ns.tar.gz', 'RC6_ns/', 'RC6.in.me'], \
    #      ['RC6_A40.tar.gz', 'RC6_A40/', 'RC6.in.me'], \
    #      ['C3.tar.gz', 'C3/', 'C3.in.me'], \
    #      ['C3_ns.tar.gz', 'C3_ns/', 'C3.in.me'], \
    #      ['C3_A40.tar.gz', 'C3_A40/', 'C3.in.me'], \
    #      ['C3_40_06FeH.tar.gz', 'C3_40_06FeH/', 'C3_40_06FeH.in.me'], \
    #      ['C6.tar.gz', 'C6/', 'C6.in.me'], \
    #      ['C6_ns.tar.gz', 'C6_ns/', 'C6.in.me'], \
    #      ['C6_A40.tar.gz', 'C6_A40/', 'C6.in.me'], \
    #      ['C6_40_06FeH.tar.gz', 'C6_40_06FeH/', 'C6_40_06FeH.in.me'], \
    #      ['C3_TRUE.tar.gz', 'C3_TRUE/', 'C3_True.in.me'], \
    #      ['C6_TRUE.tar.gz', 'C6_TRUE/', 'C6_TRUE.in.me'], \
    #      ['C6_TRUE_2.tar.gz', 'C6_TRUE_2/', 'C6_TRUE_2.in.me'], \
    #      ]
    #
    # z=[pd.DataFrame()]*len(a)
    # for i in range(len(a)):
    #     tar = tt.TAR(ext,a[i][0],a[i][1],a[i][2],r'\s+')
    #     z[i] = tar()
    #
    # APK2, GES, GES_v2, GES_ns, GES_A40, APO, APO_ns, RC3, RC3_ns, RC3_A40, \
    # RC6, RC6_ns, RC6_A40, C3, C3_ns, C3_A40, C3_40_06FeH, C6, C6_ns, C6_A40, \
    # C6_40_06FeH, C3_True, C6_True, C6_True_2 = z
    # print(C6_True_2)
    # sys.exit()

    APK2 = pd.read_csv(ext+'APOKASC/APOKASC.in.me',delimiter=r'\s+')
    GES = pd.read_csv(ext+'GES/GES.in.me',delimiter=r'\s+')
    GES_v2 = pd.read_csv(ext+'GES_v2/GES.in.me',delimiter=r'\s+')
    GES_ns = pd.read_csv(ext+'GES_ns/GES.in.me',delimiter=r'\s+')
    GES_A40 = pd.read_csv(ext+'GES_A40/GES.in.me',delimiter=r'\s+')
    APO = pd.read_csv(ext+'APO6/APOGEE.in.me',delimiter=r'\s+')
    APO_ns = pd.read_csv(ext+'APO6_ns/APOGEE.in.me',delimiter=r'\s+')
    RC3 = pd.read_csv(ext+'RC3/RC3.in.me',delimiter=r'\s+')
    RC3_ns = pd.read_csv(ext+'RC3_ns/RC3.in.me',delimiter=r'\s+')
    RC3_A40 = pd.read_csv(ext+'RC3_A40/RC3.in.me',delimiter=r'\s+')
    RC6_A40 = pd.read_csv(ext+'RC6_A40/RC6.in.me',delimiter=r'\s+')
    RC6 = pd.read_csv(ext+'RC6/RC6.in.me',delimiter=r'\s+')
    RC6_ns = pd.read_csv(ext+'RC6_ns/RC6.in.me',delimiter=r'\s+')
    C3 = pd.read_csv(ext+'C3/C3.in.me',delimiter=r'\s+')
    C3_ns = pd.read_csv(ext+'C3_ns/C3.in.me',delimiter=r'\s+')
    C3_A40 = pd.read_csv(ext+'C3_A40/C3.in.me',delimiter=r'\s+')
    C3_40_06FeH = pd.read_csv(ext+'C3_40_06FeH/C3_40_06FeH.in.me',delimiter=r'\s+')
    C6 = pd.read_csv(ext+'C6/C6.in.me',delimiter=r'\s+')
    C6_A40 = pd.read_csv(ext+'C6_A40/C6.in.me',delimiter=r'\s+')
    C6_40_06FeH = pd.read_csv(ext+'C6_40_06FeH/C6_40_06FeH.in.me',delimiter=r'\s+')
    C6_ns = pd.read_csv(ext+'C6_ns/C6.in.me',delimiter=r'\s+')

    # TRI3 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C3_self')
    TRI3 = pd.read_csv('/home/ben/Dropbox/K2Poles/Data0405/TRILEGAL_C3_self')
    # TRI3 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C3_old_stars.csv')
    # TRI6 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C6_self')
    TRI6 = pd.read_csv('/home/ben/Dropbox/K2Poles/Data0405/TRILEGAL_C6_self')
    # Kep_Sim = pd.read_csv('/home/bmr135/GA/K2Poles/Standard_kepler_field/k1.6_K15.all.out.txt',delimiter=r'\s+')
    Kep_Sim = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/Standard_kepler_field/k1.6_K15.all.out.txt',delimiter=r'\s+')

    ''' [Fe/H] correct PARAM runs '''
    C3_True = pd.read_csv(ext+'C3_TRUE/C3_True.in.me',delimiter=r'\s+')
    C6_True = pd.read_csv(ext+'C6_TRUE/C6_TRUE.in.me',delimiter=r'\s+')
    C6_True_2 = pd.read_csv(ext+'C6_TRUE_2/C6_TRUE_2.in.me',delimiter=r'\s+')
    C6_True = pd.concat([C6_True,C6_True_2],ignore_index=True)
    C6_True = C6_True.reset_index(drop=True)


    ''' Application of K2 Selection Function to Kepler Simulation '''
    Kep_Sim['Teff'] = 10**(Kep_Sim['logTe'])
    Kep_Sim['g'] = 10**(Kep_Sim['logg'])
    Kep_Sim['L'] = 10**(Kep_Sim['logL'])
    Kep_Sim['radius'] = np.sqrt(Kep_Sim['Mass'] / (Kep_Sim['g']/4.441))
    Kep_Sim['JK'] = Kep_Sim['Jmag'] - Kep_Sim['Kmag']
    Kep_Sim['Vcut'] = Kep_Sim['Kmag'] + 2*(Kep_Sim['JK']+0.14) + 0.382*np.exp(2*(Kep_Sim['JK']-0.2))
    # Kep_Sim = prop.det_prob_Kepler(Kep_Sim,'numax',3090.0,135.1)
    # Kep_Sim = Kep_Sim[ (Kep_Sim['numax'] > 10) & (Kep_Sim['numax'] < 280) & \
    # (Kep_Sim['Vcut'] > 9) & (Kep_Sim['Vcut'] < 15) & (Kep_Sim['JK'] > 0.5)]

    ''' Additional PARAM inputs that aren't output values '''
    C3 = mdf.p_in('C3',C3,'C3')
    C3_True = mdf.p_in('C3',C3_True,'C3')
    C3_ns = mdf.p_in('C3',C3_ns,'C3')
    C3_A40 = mdf.p_in('C3',C3_A40,'C3')
    C3_40_06FeH = mdf.p_in('C3_40_06FeH',C3_40_06FeH,'C3_40_06FeH')
    GES = mdf.p_in('GES',GES,'C3')
    GES_v2 = mdf.p_in('GES',GES_v2,'C3')
    GES_ns = mdf.p_in('GES',GES_ns,'C3')
    GES_A40 = mdf.p_in('GES',GES_A40,'C3')
    RC3 = mdf.p_in('RC3',RC3,'C3')
    RC3_ns = mdf.p_in('RC3',RC3_ns,'C3')
    RC3_A40 = mdf.p_in('RC3',RC3_A40,'C3')
    C6 = mdf.p_in('C6',C6,'C6')
    C6_True = mdf.p_in('C6',C6_True,'C6')
    C6_A40 = mdf.p_in('C6',C6_A40,'C6')
    C6_40_06FeH = mdf.p_in('C6_40_06FeH',C6_40_06FeH,'C6_40_06FeH')
    APO = mdf.p_in('APOGEE',APO,'C6')
    APO_ns = mdf.p_in('APOGEE',APO_ns,'C6')
    RC6 = mdf.p_in('RC6',RC6,'C6')
    RC6_A40 = mdf.p_in('RC6',RC6_A40,'C6')

    # C3 = mdf.p_in('Pre_02_2018/C3',C3,'C3')
    # C3_True = mdf.p_in('Pre_02_2018/C3',C3_True,'C3')
    # C3_ns = mdf.p_in('Pre_02_2018/C3',C3_ns,'C3')
    # C3_A40 = mdf.p_in('Pre_02_2018/C3',C3_A40,'C3')
    # C3_40_06FeH = mdf.p_in('Pre_02_2018/C3_40_06FeH',C3_40_06FeH,'C3_40_06FeH')
    # GES = mdf.p_in('Pre_02_2018/GES',GES,'C3')
    # GES_v2 = mdf.p_in('Pre_02_2018/GES',GES_v2,'C3')
    # GES_ns = mdf.p_in('Pre_02_2018/GES',GES_ns,'C3')
    # GES_A40 = mdf.p_in('Pre_02_2018/GES',GES_A40,'C3')
    # RC3 = mdf.p_in('Pre_02_2018/RC3',RC3,'C3')
    # RC3_ns = mdf.p_in('Pre_02_2018/RC3',RC3_ns,'C3')
    # RC3_A40 = mdf.p_in('Pre_02_2018/RC3',RC3_A40,'C3')
    # C6 = mdf.p_in('Pre_02_2018/C6',C6,'C6')
    # C6_True = mdf.p_in('Pre_02_2018/C6',C6_True,'C6')
    # C6_A40 = mdf.p_in('Pre_02_2018/C6',C6_A40,'C6')
    # C6_40_06FeH = mdf.p_in('Pre_02_2018/C6_40_06FeH',C6_40_06FeH,'C6_40_06FeH')
    # APO = mdf.p_in('Pre_02_2018/APOGEE',APO,'C6')
    # APO_ns = mdf.p_in('Pre_02_2018/APOGEE',APO_ns,'C6')
    # RC6 = mdf.p_in('Pre_02_2018/RC6',RC6,'C6')
    # RC6_A40 = mdf.p_in('Pre_02_2018/RC6',RC6_A40,'C6')

    # RC6_ns = p_in('RC6',RC6_ns,'C6')
    TRI3['age'] = (10**TRI3['logAge'])/1e9
    TRI3['Vcut'] = TRI3['Kmag'] + 2*(TRI3['JK']+0.14) + 0.382*np.exp(2*(TRI3['JK']-0.2))
    # TRI3['logAge'] = np.log10(TRI3['age'])
    TRI6['age'] = (10**TRI6['logAge'])/1e9
    # TRI6['logAge'] = np.log10(TRI6['age'])

    cols_to_use = APK.columns.difference(APK2.columns)
    cols_to_use = cols_to_use.union(['KIC'])
    APK2 = pd.merge(APK2,APK[cols_to_use],how='inner',on=['KIC'])
    APK2['dist'] = APK['dist']
    APK2['Z'] = APK['Z']
    APK2 = APK2[APK2['dist'] > -99.9]

    ''' Adding coordinates to K2 simulations '''
    mdf.sim_coords(TRI3,3)
    mdf.sim_coords(TRI6,6)
    mdf.sim_coords(Kep_Sim,'K')
    mdf.sim_dist(TRI3)
    mdf.sim_dist(TRI6)
    mdf.sim_dist(Kep_Sim)
    Kep_Sim = Kep_Sim[Kep_Sim['dist'] < 2e4]
    mdf.vert_dist(TRI3)
    mdf.vert_dist(TRI6)
    mdf.vert_dist(Kep_Sim)
    mdf.vert_dist(C3_A40)
    mdf.vert_dist(C3_True)
    mdf.vert_dist(C6_A40)
    mdf.vert_dist(C6_True)
    mdf.vert_dist(C3_40_06FeH)
    mdf.vert_dist(C6_40_06FeH)
    mdf.vert_dist(C3)
    mdf.vert_dist(C6)
    mdf.vert_dist(APK2)
    mdf.vert_dist(GES_A40)
    mdf.vert_dist(APO)
    mdf.vert_dist(RC3_A40)
    mdf.vert_dist(RC6_A40)

    ''' Spectro tests '''
    ges,c3 = mdf.spectro(GES,C3)
    # ges_40, c3_40 = mdf.spectro(GES_A40,C3_A40)
    ges_40, c3_40 = mdf.spectro(GES_A40,C3_40_06FeH)
    ges_v2,c3_v2 = mdf.spectro(GES,C3)
    ges_ns,c3_ns = mdf.spectro(GES_ns,C3_ns)
    apo,c6 = mdf.spectro(APO,C6)
    apo_ns,c6_ns = mdf.spectro(APO_ns,C6)
    rc3,c3_R = mdf.spectro(RC3,C3)
    rc3_ns,c3_Rns = mdf.spectro(RC3_ns,C3_ns)
    # rc3_40, c3R_40 = mdf.spectro(RC3_A40,C3_A40)
    rc3_40, c3R_40 = mdf.spectro(RC3_A40,C3_40_06FeH)
    rc6_40, c6R_40 = mdf.spectro(RC6_A40,C6_40_06FeH)
    rc6,c6_R = mdf.spectro(RC6,C6)
    # rc6_ns,c6_Rns = spectro(RC6_ns,C6_ns)
    GR3,RG3 = mdf.spectro(GES,RC3)
    AR6,RA6 = mdf.spectro(APO,RC6)

    ''' Gaia-ESO Uncertainty Comps '''
    param = ['mass','rad','age']
    x = [['mass',0.0],['rad',0.0],['age',0.0]]
    err1 = mdf.uncerts(c3_ns,param,x)
    err2 = mdf.uncerts(ges_ns,param,x)
    err3 = mdf.uncerts(ges,param,x)
    err4 = mdf.uncerts(c3,param,x)

    ''' Combine all spectroscopic data for each campaign '''
    c_three = pd.concat([ges_40,rc3_40],ignore_index=True)
    c_three.reset_index(drop=True)
    c_six = pd.concat([apo,rc6_40],ignore_index=True)
    c_six.reset_index(drop=True)
    # c3 = pd.concat([c3_40,c3R_40],ignore_index=True)
    # c3.reset_index(drop=True)
    AS = pd.concat([c_three,c_six],ignore_index=True)
    AS.reset_index(drop=True)
    # K2 = C3_A40[C3_A40['rad'] < 8]
    # C6_A40 = C6_A40[C6_A40['rad'] < 8]
    K2_True = pd.concat([C3_True,C6_True],ignore_index=True)
    K2 = pd.concat([C3_A40,C6_A40],ignore_index=True)
    # K2 = pd.concat([C3_40_06FeH,C6_40_06FeH],ignore_index=True)
    K2.reset_index(drop=True)
    K2_True.reset_index(drop=True)
    # k2 = pd.concat([c3_40,c6_40],ignore_index=True)
    # k2.reset_index(drop=True)
    TRI = pd.concat([TRI3,TRI6],ignore_index=True)
    TRI.reset_index(drop=True)

    K2['lgs'] = np.log10(27400 * (K2['nmx']/3090) * np.sqrt(K2['Teff']/5777))
    AS['lgs'] = np.log10(27400 * (AS['nmx']/3090) * np.sqrt(AS['Teff']/5777))
    mdf.vert_dist(AS)
    # AS['feh'] = AS['feh'] + 3

    ''' Photometric data for spectroscopic K2 stars '''
    cols_to_use = AS.columns.difference(K2.columns)
    cols_to_use = cols_to_use.union(['#Id'])
    K2_AS = pd.merge(K2,AS[cols_to_use],how='inner',on=['#Id'])
    K2_AS.reset_index(drop=True)

    cols_to_use = c_three.columns.difference(C3_True.columns)
    cols_to_use = cols_to_use.union(['#Id'])
    C3_AS = pd.merge(C3_True,c_three[cols_to_use],how='inner',on=['#Id'])
    C3_AS.reset_index(drop=True)

    cols_to_use = c_six.columns.difference(C6_True.columns)
    cols_to_use = cols_to_use.union(['#Id'])
    C6_AS = pd.merge(C6_True,c_six[cols_to_use],how='inner',on=['#Id'])
    C6_AS.reset_index(drop=True)

    K2_True_AS = pd.concat([C3_AS,C6_AS],ignore_index=True)
    K2_True_AS.reset_index(drop=True)

    K2_dnu = pd.DataFrame()
    K2_dnu = K2[K2['dnu'] > 0]
    K2_dnu.reset_index(drop=True)


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
    #     mdf.c3_6_param(AS,param[i],label[i])
    # plt.show()

    #     plt.savefig('/home/bmr135/Dropbox/ALL_spectro_'+param[i])


    ''' Clump removal function '''
    if clump == 1:
        TRI = TRI[TRI['radius'] < 9.0]
        K2 = K2[K2['rad'] < 9.0]
        APK2 = APK2[APK2['rad'] < 9.0]
        Kep_Sim = Kep_Sim[Kep_Sim['radius'] < 9.0]

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
    #     # print(alt_sim.columns.values)
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

    thin = TRI[TRI['#Gc'] == 1]
    thick = TRI[TRI['#Gc'] == 2]

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
    # mdf.histo(K2,'mass',bins,r'Mass [M$_{\odot}$]',0,r'K2 PARAM')
    # mdf.histo(alt_sim,'Mass',bins,r'Mass [M$_{\odot}$]',0,r'TRILEGAL')
    # mdf.histo(K2,'m_seismo',bins,r'Mass [M$_{\odot}$]',0,r'K2 S-R')
    # plt.legend(prop={'size':15})
    # plt.show()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'Mass_comp_param_sim_SR.png')

    ''' Parameter distribution plots '''
    K2_RGB = pd.DataFrame()
    APK2 = APK2[APK2['radius'] > 0.0]
    K2_RGB = K2_AS[K2_AS['rad'] < 9.0]
    # fig, axes = plt.subplots(3,2)
    # ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
    # # plt.suptitle(r'R $< 9$', fontsize=15)
    # plt.suptitle(r'$\forall$ R', fontsize=15)
    # ax0.hist(K2_True['mass'],bins=bins,histtype='step',label=r'Corrected',normed=True)
    # ax0.hist(K2['mass'],bins=bins,histtype='step',label=r'Original',normed=True)
    # ax0.legend(prop={'size':10})
    # ax0.set_xlabel(r'Mass [M$_{\odot}$]')
    # ax1.hist(K2_True['logAge'],bins=ranges[0],histtype='step',label=r'R $< 9$',normed=True)
    # ax1.hist(K2['logAge'],bins=ranges[0],histtype='step',label=r'$\forall$ R',normed=True)
    # ax1.set_xlabel(r'log$_{10}$(Age)')
    # ax2.hist(K2_True['rad'],bins=np.linspace(3,20,50),histtype='step',label=r'Spectroscopic',normed=1.0)
    # ax2.hist(K2['rad'],bins=np.linspace(3,20,50),histtype='step',label=r'K2',normed=1.0)
    # ax2.set_xlabel(r'Radius [R$_{\odot}$]')
    # ax3.hist(np.abs(K2_True['Z']),bins=ranges_kep[1],histtype='step',label=r'Spectroscopic',normed=True)
    # ax3.hist(np.abs(K2['Z']),bins=ranges_kep[1],histtype='step',label=r'K2',normed=True)
    # ax3.set_xlabel(r'Z [kpc]')
    # ax4.hist(K2_True['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'Spectroscopic',normed=1)
    # ax4.hist(K2['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'K2',normed=1)
    # ax4.set_xlabel(r'[Fe/H]')
    # ax5.axis('off')
    # # fig.subplots_adjust(top=3.0)
    # plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.show()
    # # if save_out == 1:
    # #     plt.savefig(ext_fig+folder_loc+'Kep_K2_age_distr.png')

    ''' Percentage Difference Between PARAM and Scaling Relations '''
    # df = pd.DataFrame()
    # df['percentage'] = 100*((AS['m_seismo']/AS['mass'])-1)
    # plt.figure()
    # plt.scatter(AS['mass'],df['percentage'])
    # plt.plot([0.5,max(AS['mass'])+0.1],[0,0],linewidth=2,color='k')
    # plt.xlabel(r'Mass$_{\rm{PARAM}}$')
    # plt.xlim([0.5,max(AS['mass'])+0.1])
    # plt.ylabel(r'$\%$ Difference between PARAM and Scaling-Relations')
    # plt.savefig(ext_fig+'Mass_percent_diff_spectro')
    # plt.show()

    ''' Teff/Met comparisons '''
    ds = mdf.least_squares(c3,ges)
    dsT = mdf.least_squares2(c3,ges)
    ds1 = mdf.least_squares(c3_R,rc3)
    ds1T = mdf.least_squares2(c3_R,rc3)
    ds2 = mdf.least_squares(c6_R,rc6)
    ds2T = mdf.least_squares2(c6_R,rc6)
    ds3 = mdf.least_squares(c6,apo)
    ds3T = mdf.least_squares2(c6,apo)
    df = mdf.least_squares2(GR3,RG3)
    df1 = mdf.least_squares2(AR6,RA6)
    df2 = mdf.least_squares(GR3,RG3)
    df3 = mdf.least_squares(AR6,RA6)

    x = np.linspace(4100,5500,100)
    plt.figure()
    plt.scatter(c3['Teff'],ges['Teff'],label=r'GES C3')
    plt.scatter(c3_R['Teff'],rc3['Teff'],label=r'RAVE C3',color='g')
    plt.scatter(c6_R['Teff'],rc6['Teff'],label=r'RAVE C6',color='r')
    plt.scatter(c6['Teff'],apo['Teff'],label=r'APOGEE C6',color='m')
    plt.xlabel(r'MAST T$_{\rm{eff}}$')
    plt.ylabel(r'Spectroscopic T$_{\rm{eff}}$')
    plt.plot([4100,5500],[4100,5500],c='k')
    plt.plot(x,(x*ds[0])+ds[1],linewidth=2,c='b')
    plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2,c='g')
    plt.plot(x,(x*ds2[0])+ds2[1],linewidth=2,c='r')
    plt.plot(x,(x*ds3[0])+ds3[1],linewidth=2,c='m')
    plt.xlim(4100,5500)
    plt.ylim(4100,5500)
    plt.legend(loc=4)
    plt.show()

    plt.figure()
    plt.subplot(2,2,1)
    plt.scatter(c3['feh'],ges['feh'],label=r'GES C3')
    plt.scatter(c3_R['feh'],rc3['feh'],label=r'RAVE C3',color='g')
    plt.scatter(c6_R['feh'],rc6['feh'],label=r'RAVE C6',color='r')
    plt.scatter(c6['feh'],apo['feh'],label=r'APOGEE C6',color='m')
    plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    plt.plot(x,(x*ds[0])+ds[1],linewidth=2,c='b')
    plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2,c='g')
    plt.plot(x,(x*ds2[0])+ds2[1],linewidth=2,c='r')
    plt.plot(x,(x*ds3[0])+ds3[1],linewidth=2,c='m')
    plt.xlim(-3.0,1.0)
    plt.ylim(-3.0,1.0)
    plt.tick_params(labelsize=15)
    plt.title(r'All Pipelines',fontsize=20)

    plt.subplot(2,2,2)
    plt.scatter(c3['feh'],ges['feh'],label=r'C3')
    plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    plt.plot(x,(x*ds[0])+ds[1],linewidth=2,c='b')
    plt.xlim(-3.0,1.0)
    plt.ylim(-3.0,1.0)
    plt.text(-.1, -2.75,r'Fit: %.6sx $+$ %.6s' %(ds[0][0],ds[1][0]), ha='center', va='center',fontsize=15)
    plt.tick_params(labelsize=15)
    plt.title(r'Gaia-ESO',fontsize=20)

    plt.subplot(2,2,3)
    plt.scatter(c3_R['feh'],rc3['feh'],label=r'C3',color='g')
    plt.scatter(c6_R['feh'],rc6['feh'],label=r'C6',color='r')
    plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2,c='g')
    plt.plot(x,(x*ds2[0])+ds2[1],linewidth=2,c='r')
    plt.xlim(-3.0,1.0)
    plt.ylim(-3.0,1.0)
    plt.text(-.1, -2.5,r'Fit RC3: %.6sx $+$ %.6s' %(ds1[0][0],ds1[1][0]), ha='center', va='center',fontsize=15)
    plt.text(-.1, -2.75,r'Fit RC6: %.6sx $+$ %.6s' %(ds2[0][0],ds2[1][0]), ha='center', va='center',fontsize=15)
    plt.tick_params(labelsize=15)
    plt.title(r'RAVE',fontsize=20)
    plt.legend()

    plt.subplot(2,2,4)
    plt.scatter(c6['feh'],apo['feh'],label=r'C6',color='m')
    plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    plt.plot(x,(x*ds3[0])+ds3[1],linewidth=2,c='m')
    plt.xlim(-3.0,1.0)
    plt.ylim(-3.0,1.0)
    plt.tick_params(labelsize=15)
    plt.title(r'APOGEE',fontsize=20)
    plt.text(-.1, -2.75,r'Fit: %.6sx $+$ %.6s' %(ds3[0][0],ds3[1][0]), ha='center', va='center',fontsize=15)

    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.subplot(2,1,1)
    plt.scatter(GR3['Teff'],RG3['Teff'])
    plt.xlabel(r'Gaia-ESO T$_{\rm{eff}}$')
    plt.ylabel(r'RAVE T$_{\rm{eff}}$')
    plt.plot([4100,5500],[4100,5500],c='k')
    # plt.plot(x,(x*df2[0])+df2[1],linewidth=2)
    plt.xlim(4100,5500)
    plt.title(r'C3')
    plt.subplot(2,1,2)
    plt.scatter(AR6['Teff'],RA6['Teff'])
    plt.xlabel(r'APOGEE T$_{\rm{eff}}$')
    plt.ylabel(r'RAVE T$_{\rm{eff}}$')
    plt.plot([4100,5500],[4100,5500],c='k')
    # plt.plot(x,(x*df3[0])+df3[1],linewidth=2)
    plt.xlim(4100,5500)
    plt.title(r'C6')

    plt.tight_layout()
    plt.show()

    ''' Spectroscopic and Photometric Parameter Comparison Plots - FeH/Teff '''
    # plt.figure()
    # plt.subplot(2,1,1)
    # plt.scatter(ges['Teff'],ges['Teff']-c3['Teff'])
    # plt.xlabel(r'Gaia-ESO T$_{\rm{eff}}$')
    # plt.ylabel(r'$\Delta$T$_{\rm{eff}}$ (Gaia-ESO - Photom)')
    # plt.plot([4100,5500],[0,0],c='k')
    # plt.plot(x,(x*dsT[0])+dsT[1],linewidth=2)
    # plt.xlim(4100,5500)
    # plt.title(r'C3')
    # plt.text(5300, 400,r'Offset: %.4sK' %(dsT[1][0]), ha='center', va='center')
    #
    # plt.subplot(2,1,2)
    # plt.scatter(rc3['Teff'],rc3['Teff']-c3_R['Teff'])
    # plt.xlabel(r'RAVE T$_{\rm{eff}}$')
    # plt.ylabel(r'$\Delta$T$_{\rm{eff}}$ (RAVE - Photom)')
    # plt.plot([4100,5500],[0,0],c='k')
    # plt.plot(x,(x*ds1T[0])+ds1T[1],linewidth=2)
    # plt.xlim(4100,5500)
    # plt.title(r'C3')
    # plt.text(5300, 175,r'Offset: %.4sK' %(ds1T[1][0]), ha='center', va='center')
    #
    # plt.tight_layout()
    # plt.show()
    #
    # x = np.linspace(-1.5,.5,100)
    # plt.figure()
    # plt.subplot(2,1,1)
    # plt.scatter(ges['feh'],ges['feh']-c3['feh'])
    # plt.xlabel(r'Gaia-ESO [Fe/H]')
    # plt.ylabel(r'$\Delta$[Fe/H] (Gaia-ESO - Photom)')
    # plt.plot([-1.5,0.5],[0,0],c='k')
    # plt.plot(x,(x*ds[0])+ds[1],linewidth=2)
    # plt.xlim(-1.5,0.5)
    # plt.title(r'C3')
    # plt.text(-1.25, 1.25,r'Offset: %.4sdex' %(ds[1][0]), ha='center', va='center')
    #
    # plt.subplot(2,1,2)
    # plt.scatter(rc3['feh'],rc3['feh']-c3_R['feh'])
    # plt.xlabel(r'RAVE [Fe/H]')
    # plt.ylabel(r'$\Delta$[Fe/H] (RAVE - Photom)')
    # plt.plot([-1.5,0.5],[0,0],c='k')
    # plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2)
    # plt.xlim(-1.5,0.5)
    # plt.title(r'C3')
    # plt.text(-1.25, 0.8,r'Offset: %.5sdex' %(ds1[1][0]), ha='center', va='center')
    #
    # plt.tight_layout()
    # plt.show()

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


    ''' Radius Distributions '''
    # plt.figure()
    # mdf.histo(K2_AS,'rad',np.linspace(0,20,100),r'Radius [R$_{\odot}$]',0,r'K2')
    # mdf.histo(AS,'rad',np.linspace(0,20,100),r'Radius [R$_{\odot}$]',0,r'TRILEGAL')
    # plt.legend(prop={'size':15})
    # plt.show()
    # if save_out == 1:
    # plt.savefig(ext_fig+'K2_sim_radius_spectro.png')

    ''' KEPLER vs K2 simulations '''
    # plt.figure()
    # for i in range(100):
    #     Kep_Sim2 = alt_sim_kep.sample(n=len(alt_sim))
    #     Kep_Sim2 = Kep_Sim2.reset_index(drop=True)
    #     plt.hist(Kep_Sim2['logAge'][np.isfinite(Kep_Sim2['logAge'])],bins=ranges[0],histtype='step',color='b',label=r'Kepler Sim',normed=True)
    #     plt.hist(alt_sim['logAge'][np.isfinite(alt_sim['logAge'])],bins=ranges[0],histtype='step',color='orange',label=r'K2 Sim',normed=True)
    # plt.legend([r'Kepler Sim',r'K2 Sim'],prop={'size':15},loc=2)
    # plt.xlabel(r'log$_{10}$(Age)')
    # # plt.title(r'$\forall$ R')
    # plt.title(r'R $< 9$')
    # plt.show()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'Kep_multi_samp_K2_age_distr.png')
    #
    # thin_kep = Kep_Sim2[Kep_Sim2['#Gc'] == 1]
    # thick_kep = Kep_Sim2[Kep_Sim2['#Gc'] == 2]
    #
    # # plt.figure()
    # fig, axes = plt.subplots(2,1,sharex=True,sharey=True)
    # ax0,ax1 = axes.flatten()
    # ax0.hist([thick['logAge'],thin['logAge']],bins=ranges[0],stacked=True,label=[r'K2 Sim Thick',r'K2 Sim Thin'])
    # # plt.xlabel(r'log$_{10}$(Age)')
    # ax0.legend(prop={'size':15},loc=2)
    # ax1.hist([thick_kep['logAge'],thin_kep['logAge']],bins=ranges[0],stacked=True,label=[r'Kepler Sim Thick',r'Kepler Sim Thin'])
    # plt.xlabel(r'log$_{10}$(Age)')
    # ax1.legend(prop={'size':15},loc=2)
    # plt.show()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'Kep_K2_age_distr.png')

    ''' Moving average for age-metallicity/alpha trends '''
    # plt.figure()
    # c_three = c_three.sort_values(['age'])
    # c_three = c_three[c_three['feh'] > -5]
    # c_three = c_three[c_three['alpha'] > -5]
    # c_six = c_six.sort_values(['age'])
    # c_six = c_six[c_six['feh'] > -5]
    # c_six = c_six[c_six['alpha'] > -5]
    #
    # plt.scatter(c_three['age'],c_three['feh'])
    # y_av = mdf.movingaverage(c_three['feh'],50)
    # plt.plot(c_three['age'], y_av,'r')
    # plt.xlabel(r'Age [Gyr]')
    # plt.ylabel(r'[Fe/H]')
    # plt.show()

    ''' Kiel Diagram '''
    # plt.figure()
    # plt.scatter(K2['Teff'],K2['lgs'],c=K2['age'],cmap=colormaps.parula)
    # # plt.scatter(AS['Teff'],AS['lgs'],c=AS['feh'],cmap=colormaps.parula)
    # # plt.scatter(AS['Teff'],AS['lgs'],color='k')
    # cbar = plt.colorbar()
    # cbar.set_label(r'Age [Gyr]', rotation=270, fontsize=15, labelpad=25)
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # plt.xlabel(r'T$_{\rm{eff}}$ [K]', fontsize=15)
    # plt.ylabel(r'log$_{10}$(g)', fontsize=15)
    # plt.tight_layout()
    # plt.show()
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
    C3_True['Rs'] = (C3_True['nmx']/3090) * (C3_True['dnu']/135.1)**-2 * (C3_True['Teff']/5777)**0.5
    C6_True['Rs'] = (C6_True['nmx']/3090) * (C6_True['dnu']/135.1)**-2 * (C6_True['Teff']/5777)**0.5
    K2['Rs'] = (K2['nmx']/3090) * (K2['dnu']/135.1)**-2 * (K2['Teff']/5777)**0.5
    AS['Rs'] = (AS['nmx']/3090) * (AS['dnu']/135.1)**-2 * (AS['Teff']/5777)**0.5
    K2_AS['Rs'] = (K2_AS['nmx']/3090) * (K2_AS['dnu']/135.1)**-2 * (K2_AS['Teff']/5777)**0.5
    plt.figure()
    plt.scatter(K2['rad'],K2['rad']-K2['Rs'])
    plt.plot([0,max(K2['rad'])+0.1],[0,0])
    plt.xlim(0,max(K2['rad'])+0.1)
    plt.xlabel(r'Radius [R$_{\odot}$]', fontsize=15)
    plt.ylabel(r'R - R$_{sr}$', fontsize=15)
    plt.tight_layout()
    plt.show()

    plt.figure()
    hist, bins, patches = plt.hist(C3_True['rad'],bins=50,histtype='step',label=r'PARAM R',normed=True,linewidth=2)
    plt.hist(C3_True['Rs'],bins=bins,histtype='step',label=r'Scaling R',normed=True,linewidth=2)
    # plt.hist(alt_sim['radius'],bins=bins,histtype='step',label=r'TRI R',normed=True,linewidth=2,alpha=0.5)
    plt.xlabel(r'Radius [R$_{\odot}$]', fontsize=15)
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure()
    hist, bins, patches = plt.hist(C6_True['rad'],bins=50,histtype='step',label=r'PARAM R',normed=True,linewidth=2)
    plt.hist(C6_True['Rs'],bins=bins,histtype='step',label=r'Scaling R',normed=True,linewidth=2)
    # plt.hist(alt_sim['radius'],bins=bins,histtype='step',label=r'TRI R',normed=True,linewidth=2,alpha=0.5)
    plt.xlabel(r'Radius [R$_{\odot}$]', fontsize=15)
    plt.tight_layout()
    plt.legend()
    plt.show()

    plt.figure()
    hist, bins, patches = plt.hist(K2_AS['Rs'],bins=50,histtype='step',label=r'Photometry',normed=True,linewidth=2)
    plt.hist(AS['Rs'],bins=bins,histtype='step',label=r'Spectroscopy',normed=True,linewidth=2)
    # plt.xlabel(r'Scaling Relation Radius [R$_{\odot}$]', fontsize=15)
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure()
    hist, bins, patches = plt.hist(K2_AS['Teff'],bins=50,histtype='step',label=r'Photometry',normed=True,linewidth=2)
    plt.hist(AS['Teff'],bins=bins,histtype='step',label=r'Spectroscopy',normed=True,linewidth=2)
    plt.xlabel(r'Teff [K]', fontsize=15)
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.scatter(AS['logAge'],AS['rad'],c=AS['Teff'],cmap=colormaps.parula)
    cbar = plt.colorbar()
    cbar.set_label(r'T$_{\rm{eff}}$ [K]', rotation=270, fontsize=15, labelpad=25)
    plt.xlabel(r'log Age', fontsize=15)
    plt.ylabel(r'Radius [R$_{\odot}$]', fontsize=15)
    plt.tight_layout()
    plt.show()


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
    # mdf.space_density2(C6_True)
    # mdf.space_density2(AS)
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
    # # print(m_err,m_per,z_per,z_err)
    # K2['age'] = 10**K2['logAge'] * 1e-9
    # # K2 = K2[K2['alpha'] > -4]
    # bins = np.linspace(0.75,2.25,10)
    # mass_Z, edges, number = scipy.stats.binned_statistic(C3_True['mass'],C3_True['feh'],statistic='median',bins=bins)
    # mass_Z6, edges6, number6 = scipy.stats.binned_statistic(C6_True['mass'],C6_True['feh'],statistic='median',bins=bins)
    # print(mass_Z)
    # print(mass_Z6)
    #
    # plt.figure()
    # plt.scatter(K2_True['mass'],K2_True['Z'],c=K2_True['feh'],cmap=colormaps.parula)
    # plt.errorbar(2.25, -3, xerr=m_err, yerr=z_err*1e-3,color='k')
    # cbar = plt.colorbar()
    # cbar.set_label(r'Age [Gyr]', rotation=270, fontsize=15, labelpad=25)
    # plt.plot(edges[:-1],mass_Z,color='k',linewidth=2)
    # plt.plot(edges6[:-1],abs(mass_Z6),color='k',linewidth=2)
    # plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    # plt.ylabel(r'Z [kpc]', fontsize=15)
    # plt.tight_layout()
    # plt.show()

    ''' [Fe/H] vs [Alpha/Fe] '''
    # AS = AS[AS['alpha'] > -4]
    # plt.figure()
    # plt.scatter(AS['feh'],AS['alpha'],c=AS['age'],cmap=colormaps.parula)
    # cbar = plt.colorbar()
    # cbar.set_label(r'Age [Gyr]', rotation=270, fontsize=15, labelpad=25)
    # plt.xlabel(r'[Fe/H]', fontsize=15)
    # plt.ylabel(r'[$\alpha$/Fe]', fontsize=15)
    # plt.tight_layout()
    # plt.show()

    ''' Mag vs Z with trend lines '''
    # bins = np.linspace(-8,0,16)
    # bins6 = np.linspace(0,8,16)
    # Kp_Z, edges, number = scipy.stats.binned_statistic(C3_40_06FeH['Z'],C3_40_06FeH['Kepler'],statistic='median',bins=bins)
    # Kp_Z6, edges6, number6 = scipy.stats.binned_statistic(C6_40_06FeH['Z'],C6_40_06FeH['Kepler'],statistic='median',bins=bins6)
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
