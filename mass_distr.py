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

    ext = '/home/bmr135/GA/K2Poles/param_outputs/Poles/'
    ext_fig = '/home/bmr135/Dropbox/K2Poles/pop_trends/051017/'
    # ext = '/media/ben/SAMSUNG/GA/K2Poles/param_outputs/Poles/'
    # ext_fig = '/home/ben/Dropbox/K2Poles/pop_trends/051017/'

    print('Figures saved to: '+ ext_fig + folder_loc)

    save_out = 0 # 0 - no saving; 1 - save all figures

    ''' Read in files '''
    APK = pd.read_csv('/home/bmr135/GA/K2Poles/APOKASC4BEN.txt')
    # APK = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/APOKASC4BEN.txt')

    a = [\
        #  ['ben_k2.tgz', 'APO_080218.in_STD_20Gyr/', 'APOKASC.in.me'], \
         ['ben_k2.tgz', 'GES_080218.in_STD_20Gyr/', 'GES_080218.in.me'], \
         ['ben_k2.tgz', 'APO_080218.in_STD_20Gyr/','APO_080218.in.me'], \
         ['ben_k2.tgz', 'RC3_080218.in_STD_20Gyr/', 'RC3_080218.in.me'], \
         ['ben_k2.tgz', 'RC6_080218.in_STD_20Gyr/', 'RC6_080218.in.me'], \
         ['ben_k2.tgz', 'C3_080218.in_STD_20Gyr/', 'C3_080218.in.me'], \
         ['ben_k2.tgz', 'C6_080218_1.in_STD_20Gyr/', 'C6_080218_1.in.me'], \
         ['ben_k2.tgz', 'C6_080218_2.in_STD_20Gyr/', 'C6_080218_2.in.me'], \
         ['C3_C6.tar.gz', 'C3_new_FeH/', 'C3_220218.in.me'], \
         ['C3_C6.tar.gz', 'C6_new_FeH/', 'C6_220218.in.me'], \
         ]

    z=[pd.DataFrame()]*len(a)
    for i in range(len(a)):
        tar = tt.TAR(ext,a[i][0],a[i][1],a[i][2],r'\s+')
        z[i] = tar()

    GES_T, APO_T, RC3_T, \
    RC6_T, C3, C6_1, C6_2, \
    C3_New, C6_New = z
    # sys.exit()

    APK2 = pd.read_csv(ext+'APOKASC/APOKASC.in.me',delimiter=r'\s+')

    TRI3 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C3_self')
    # TRI3 = pd.read_csv('/home/ben/Dropbox/K2Poles/Data0405/TRILEGAL_C3_self')
    # TRI3 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C3_old_stars.csv')
    TRI6 = pd.read_csv('/home/bmr135/Dropbox/K2Poles/Data0405/TRILEGAL_C6_self')
    # TRI6 = pd.read_csv('/home/ben/Dropbox/K2Poles/Data0405/TRILEGAL_C6_self')
    Kep_Sim = pd.read_csv('/home/bmr135/GA/K2Poles/Standard_kepler_field/k1.6_K15.all.out.txt',delimiter=r'\s+')
    # Kep_Sim = pd.read_csv('/media/ben/SAMSUNG/GA/K2Poles/Standard_kepler_field/k1.6_K15.all.out.txt',delimiter=r'\s+')

    ''' [Fe/H] correct PARAM runs '''
    # C3 = pd.read_csv(ext+'C3_080218.in.me',delimiter=r'\s+')
    # C6 = pd.read_csv(ext+'C6_080218.in.me',delimiter=r'\s+')
    # GES_T = pd.read_csv(ext+'GES_080218.in.me',delimiter=r'\s+')
    # RC3_T = pd.read_csv(ext+'RC3_080218.in.me',delimiter=r'\s+')
    # RC6_T = pd.read_csv(ext+'RC6_080218.in.me',delimiter=r'\s+')
    L3_T = pd.read_csv(ext+'L3_080218.in.me',delimiter=r'\s+')
    L6_T = pd.read_csv(ext+'L6_080218.in.me',delimiter=r'\s+')
    # APO_T = pd.read_csv(ext+'APO_080218.in.me',delimiter=r'\s+')

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
    C3 = mdf.p_in('C3_080218',C3,'C3')
    C3_New = mdf.p_in('C3_220218',C3_New,'C3')
    GES = mdf.p_in('GES_080218',GES_T,'C3')
    RC3 = mdf.p_in('RC3_080218',RC3_T,'C3')
    L3 = mdf.p_in('L3_080218',L3_T,'C3')
    C6_1 = mdf.p_in('C6_080218_1',C6_1,'C6')
    C6_2 = mdf.p_in('C6_080218_2',C6_2,'C6')
    C6_New = mdf.p_in('C6_220218',C6_New,'C6')
    APO = mdf.p_in('APO_080218',APO_T,'C6')
    RC6 = mdf.p_in('RC6_080218',RC6_T,'C6')
    L6 = mdf.p_in('L6_080218',L6_T,'C6')

    C6 = pd.concat([C6_1,C6_2],ignore_index=True)

    TRI3['age'] = (10**TRI3['logAge'])/1e9
    TRI3['Vcut'] = TRI3['Kmag'] + 2*(TRI3['JK']+0.14) + 0.382*np.exp(2*(TRI3['JK']-0.2))
    TRI6['age'] = (10**TRI6['logAge'])/1e9

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
    mdf.vert_dist(C3)
    mdf.vert_dist(C3_New)
    mdf.vert_dist(C6)
    mdf.vert_dist(C6_New)
    mdf.vert_dist(APK2)
    mdf.vert_dist(GES)
    mdf.vert_dist(APO)
    mdf.vert_dist(RC3)
    mdf.vert_dist(RC6)
    mdf.vert_dist(L3)
    mdf.vert_dist(L6)

    ''' Spectro tests '''
    ges,c3 = mdf.spectro(GES,C3)
    apo,c6 = mdf.spectro(APO,C6)
    rc3,c3_R = mdf.spectro(RC3,C3)
    rc6,c6_R = mdf.spectro(RC6,C6)
    l3, c3_l = mdf.spectro(L3,C3)
    l6, c6_l = mdf.spectro(L6,C6)
    GR3,RG3 = mdf.spectro(GES,RC3)
    AR6,RA6 = mdf.spectro(APO,RC6)

    ''' Gaia-ESO Uncertainty Comps '''
    param = ['mass','rad','age']
    x = [['mass',0.0],['rad',0.0],['age',0.0]]
    err1 = mdf.uncerts(C3,param,x)
    err2 = mdf.uncerts(GES,param,x)
    err3 = mdf.uncerts(ges,param,x)
    err4 = mdf.uncerts(c3,param,x)

    ''' Combine all spectroscopic data for each campaign '''
    c_three = pd.concat([ges,rc3,l3],ignore_index=True)
    c_three.reset_index(drop=True)
    c_six = pd.concat([apo,rc6,l6],ignore_index=True)
    c_six.reset_index(drop=True)
    AS = pd.concat([c_three,c_six],ignore_index=True)
    AS.reset_index(drop=True)

    K2 = pd.concat([C3,C6],ignore_index=True)
    K2.reset_index(drop=True)
    K2.reset_index(drop=True)

    K2_New = pd.concat([C3_New,C6_New],ignore_index=True)
    K2_New.reset_index(drop=True)
    K2_New.reset_index(drop=True)

    TRI = pd.concat([TRI3,TRI6],ignore_index=True)
    TRI.reset_index(drop=True)

    K2['lgs'] = np.log10(27400 * (K2['nmx']/3090) * np.sqrt(K2['Teff']/5777))
    K2_New['lgs'] = np.log10(27400 * (K2_New['nmx']/3090) * np.sqrt(K2_New['Teff']/5777))
    AS['lgs'] = np.log10(27400 * (AS['nmx']/3090) * np.sqrt(AS['Teff']/5777))
    mdf.vert_dist(AS)
    # AS['feh'] = AS['feh'] + 3

    ''' Photometric data for spectroscopic K2 stars '''
    cols_to_use = AS.columns.difference(K2.columns)
    cols_to_use = cols_to_use.union(['#Id'])
    K2_AS = pd.merge(K2,AS[cols_to_use],how='inner',on=['#Id'])
    K2_AS.reset_index(drop=True)

    cols_to_use = c_three.columns.difference(C3.columns)
    cols_to_use = cols_to_use.union(['#Id'])
    C3_AS = pd.merge(C3,c_three[cols_to_use],how='inner',on=['#Id'])
    C3_AS.reset_index(drop=True)

    cols_to_use = c_six.columns.difference(C6.columns)
    cols_to_use = cols_to_use.union(['#Id'])
    C6_AS = pd.merge(C6,c_six[cols_to_use],how='inner',on=['#Id'])
    C6_AS.reset_index(drop=True)

    K2_AS = pd.concat([C3_AS,C6_AS],ignore_index=True)
    K2_AS.reset_index(drop=True)

    K2_dnu = pd.DataFrame()
    K2_dnu = K2[K2['dnu'] > 0]
    K2_dnu.reset_index(drop=True)

    AS_red = AS[AS['logAge'] < 9.8]
    AS_red = AS_red[AS_red['logAge'] > 9.6]
    AS_red.reset_index(drop=True)


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
    # mdf.histo(K2,'mass',bins,r'Mass [M$_{\odot}$]',0,r'K2 PARAM')
    # mdf.histo(alt_sim,'Mass',bins,r'Mass [M$_{\odot}$]',0,r'TRILEGAL')
    # mdf.histo(K2,'m_seismo',bins,r'Mass [M$_{\odot}$]',0,r'K2 S-R')
    # plt.legend(prop={'size':15})
    # plt.show()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'Mass_comp_param_sim_SR.png')

    ''' Parameter distribution plots '''
    # K2_RGB = pd.DataFrame()
    # APK2 = APK2[APK2['radius'] > 0.0]
    # K2_New = K2_New[K2_New['logAge'] > 9.9]
    # AS = AS[AS['logAge'] > 9.9]
    fig, axes = plt.subplots(3,2)
    ax0,ax1,ax2,ax3,ax4,ax5 = axes.flatten()
    # plt.suptitle(r'log$_{10}$(Age) $> 10.0$', fontsize=15)
    # plt.suptitle(r'$\forall$ R', fontsize=15)
    ax0.hist(K2['mass'],bins=np.linspace(0.5,2.5,20),histtype='step',label=r'Photometric',normed=True)
    ax0.hist(K2_New['mass'],bins=np.linspace(0.5,2.5,20),histtype='step',label=r'Spectro. based Met.',normed=True)
    ax0.legend(prop={'size':10})
    ax0.set_xlabel(r'Mass [M$_{\odot}$]')
    ax1.hist(K2['logAge'],bins=np.linspace(8.5,10.5,75),histtype='step',label=r'R $< 9$',normed=True)
    ax1.hist(K2_New['logAge'],bins=np.linspace(8.5,10.5,75),histtype='step',label=r'$\forall$ R',normed=True)
    ax1.set_xlabel(r'log$_{10}$(Age)')
    ax2.hist(K2['rad'],bins=np.linspace(3,20,50),histtype='step',label=r'Spectroscopic')#,normed=1.0)
    ax2.hist(K2_New['rad'],bins=np.linspace(3,20,50),histtype='step',label=r'K2')#,normed=1.0)
    ax2.set_xlabel(r'Radius [R$_{\odot}$]')
    ax3.hist(K2['Z'],bins=np.linspace(-8,8,130),histtype='step',label=r'Spectroscopic',normed=True)
    ax3.hist(K2_New['Z'],bins=np.linspace(-8,8,130),histtype='step',label=r'K2',normed=True)
    ax3.set_xlabel(r'Z [kpc]')
    ax4.hist(K2['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'Spectroscopic',normed=1)
    ax4.hist(K2_New['feh'],bins=np.linspace(-2,0.75,30),histtype='step',label=r'K2',normed=1)
    ax4.set_xlabel(r'[Fe/H]')
    ax5.axis('off')
    # fig.subplots_adjust(top=3.0)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    # if save_out == 1:
    #     plt.savefig(ext_fig+folder_loc+'Kep_K2_age_distr.png')

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
    # ds = mdf.least_squares(c3,ges)
    # dsT = mdf.least_squares2(c3,ges)
    # ds1 = mdf.least_squares(c3_R,rc3)
    # ds1T = mdf.least_squares2(c3_R,rc3)
    # ds2 = mdf.least_squares(c6_R,rc6)
    # ds2T = mdf.least_squares2(c6_R,rc6)
    # ds3 = mdf.least_squares(c6,apo)
    # ds3T = mdf.least_squares2(c6,apo)
    # df = mdf.least_squares2(GR3,RG3)
    # df1 = mdf.least_squares2(AR6,RA6)
    # df2 = mdf.least_squares(GR3,RG3)
    # df3 = mdf.least_squares(AR6,RA6)
    #
    # x = np.linspace(4100,5500,100)
    # plt.figure()
    # plt.scatter(c3['Teff'],ges['Teff'],label=r'GES C3')
    # plt.scatter(c3_R['Teff'],rc3['Teff'],label=r'RAVE C3',color='g')
    # plt.scatter(c6_R['Teff'],rc6['Teff'],label=r'RAVE C6',color='r')
    # plt.scatter(c6['Teff'],apo['Teff'],label=r'APOGEE C6',color='m')
    # plt.xlabel(r'MAST T$_{\rm{eff}}$')
    # plt.ylabel(r'Spectroscopic T$_{\rm{eff}}$')
    # plt.plot([4100,5500],[4100,5500],c='k')
    # plt.plot(x,(x*ds[0])+ds[1],linewidth=2,c='b')
    # plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2,c='g')
    # plt.plot(x,(x*ds2[0])+ds2[1],linewidth=2,c='r')
    # plt.plot(x,(x*ds3[0])+ds3[1],linewidth=2,c='m')
    # plt.xlim(4100,5500)
    # plt.ylim(4100,5500)
    # plt.legend(loc=4)
    # plt.show()
    #
    # plt.figure()
    # plt.subplot(2,2,1)
    # plt.scatter(c3['feh'],ges['feh'],label=r'GES C3')
    # plt.scatter(c3_R['feh'],rc3['feh'],label=r'RAVE C3',color='g')
    # plt.scatter(c6_R['feh'],rc6['feh'],label=r'RAVE C6',color='r')
    # plt.scatter(c6['feh'],apo['feh'],label=r'APOGEE C6',color='m')
    # plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    # plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    # plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    # plt.plot(x,(x*ds[0])+ds[1],linewidth=2,c='b')
    # plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2,c='g')
    # plt.plot(x,(x*ds2[0])+ds2[1],linewidth=2,c='r')
    # plt.plot(x,(x*ds3[0])+ds3[1],linewidth=2,c='m')
    # plt.xlim(-3.0,1.0)
    # plt.ylim(-3.0,1.0)
    # plt.tick_params(labelsize=15)
    # plt.title(r'All Pipelines',fontsize=20)
    #
    # plt.subplot(2,2,2)
    # plt.scatter(c3['feh'],ges['feh'],label=r'C3')
    # plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    # plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    # plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    # plt.plot(x,(x*ds[0])+ds[1],linewidth=2,c='b')
    # plt.xlim(-3.0,1.0)
    # plt.ylim(-3.0,1.0)
    # plt.text(-.1, -2.75,r'Fit: %.6sx $+$ %.6s' %(ds[0][0],ds[1][0]), ha='center', va='center',fontsize=15)
    # plt.tick_params(labelsize=15)
    # plt.title(r'Gaia-ESO',fontsize=20)
    #
    # plt.subplot(2,2,3)
    # plt.scatter(c3_R['feh'],rc3['feh'],label=r'C3',color='g')
    # plt.scatter(c6_R['feh'],rc6['feh'],label=r'C6',color='r')
    # plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    # plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    # plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    # plt.plot(x,(x*ds1[0])+ds1[1],linewidth=2,c='g')
    # plt.plot(x,(x*ds2[0])+ds2[1],linewidth=2,c='r')
    # plt.xlim(-3.0,1.0)
    # plt.ylim(-3.0,1.0)
    # plt.text(-.1, -2.5,r'Fit RC3: %.6sx $+$ %.6s' %(ds1[0][0],ds1[1][0]), ha='center', va='center',fontsize=15)
    # plt.text(-.1, -2.75,r'Fit RC6: %.6sx $+$ %.6s' %(ds2[0][0],ds2[1][0]), ha='center', va='center',fontsize=15)
    # plt.tick_params(labelsize=15)
    # plt.title(r'RAVE',fontsize=20)
    # plt.legend()
    #
    # plt.subplot(2,2,4)
    # plt.scatter(c6['feh'],apo['feh'],label=r'C6',color='m')
    # plt.xlabel(r'MAST [Fe/H]',fontsize=20)
    # plt.ylabel(r'Spectroscopic [Fe/H]',fontsize=20)
    # plt.plot([-3.0,1.0],[-3.0,1.0],c='k')
    # plt.plot(x,(x*ds3[0])+ds3[1],linewidth=2,c='m')
    # plt.xlim(-3.0,1.0)
    # plt.ylim(-3.0,1.0)
    # plt.tick_params(labelsize=15)
    # plt.title(r'APOGEE',fontsize=20)
    # plt.text(-.1, -2.75,r'Fit: %.6sx $+$ %.6s' %(ds3[0][0],ds3[1][0]), ha='center', va='center',fontsize=15)
    #
    # plt.tight_layout()
    # plt.show()
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
    # mdf.histo(AS,'rad',np.linspace(0,20,100),r'Radius [R$_{\odot}$]',0,r'Spectro')
    # mdf.histo(AS_red,'rad',np.linspace(0,20,100),r'Radius [R$_{\odot}$]',0,r'Spectro Age Peak')
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
    # C3['Rs'] = (C3['nmx']/3090) * (C3['dnu']/135.1)**-2 * (C3['Teff']/5777)**0.5
    # C6['Rs'] = (C6['nmx']/3090) * (C6['dnu']/135.1)**-2 * (C6['Teff']/5777)**0.5
    K2_New['Rs'] = (K2['nmx']/3090) * (K2['dnu']/135.1)**-2 * (K2['Teff']/5777)**0.5
    # AS['Rs'] = (AS['nmx']/3090) * (AS['dnu']/135.1)**-2 * (AS['Teff']/5777)**0.5
    # K2_AS['Rs'] = (K2_AS['nmx']/3090) * (K2_AS['dnu']/135.1)**-2 * (K2_AS['Teff']/5777)**0.5
    K2_New['drad'] = (K2_New['rad']-K2_New['Rs'])/K2_New['rad']
    # K2_New = K2_New[K2_New['drad'] > -0.5]
    # K2_New = K2_New[K2_New['drad'] < 0.5]
    plt.figure()
    plt.scatter(K2_New['rad'],K2_New['drad'])
    plt.plot([0,max(K2_New['rad'])+0.1],[0,0])
    plt.xlim(0,max(K2_New['rad'])+0.1)
    plt.plot([0,max(K2_New['rad'])+0.1],[0.5,0.5],color='r',linewidth=3)
    plt.plot([0,max(K2_New['rad'])+0.1],[-0.5,-0.5],color='r',linewidth=3)
    plt.xlabel(r'Radius [R$_{\odot}$]', fontsize=15)
    plt.ylabel(r'R - R$_{sr}$', fontsize=15)
    plt.tight_layout()
    # plt.show()

    K2_New['Ms'] = (K2_New['nmx']/3090)**3 * (K2_New['dnu']/135.1)**-4 * (K2_New['Teff']/5777)**1.5
    plt.figure()
    plt.scatter(K2_New['mass'],(K2_New['mass']-K2_New['Ms'])/K2_New['mass'])
    plt.plot([0,max(K2_New['mass'])+0.1],[0,0])
    plt.xlim(0,max(K2_New['mass'])+0.1)
    plt.plot([0,max(K2_New['mass'])+0.1],[0.5,0.5],color='r',linewidth=3)
    plt.plot([0,max(K2_New['mass'])+0.1],[-0.5,-0.5],color='r',linewidth=3)
    plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    plt.ylabel(r'M - M$_{sr}$', fontsize=15)
    plt.tight_layout()
    plt.show()

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
    # mdf.space_density2(C6)
    # # mdf.space_density2(AS)
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
    f = plt.figure()
    # plt.scatter(C6_New['mass'],C6_New['Z'],label=r'Full Sample')
    x = np.linspace(0,max(C6_New['mass'])+0.05)
    plt.fill_between(x, 0.1, max(APK2['Z']), facecolor='gray', alpha=0.2, interpolate=True,label=r'Kepler Z range')
    plt.scatter(C6_New['mass'],C6_New['Z'],c=C6_New['feh'],cmap=colormaps.parula,label=None)
    # plt.errorbar(2.25, -3, xerr=m_err, yerr=z_err*1e-3,color='k')
    cbar = plt.colorbar()
    cbar.set_label(r'[Fe/H]', rotation=270, fontsize=15, labelpad=25)
    # plt.plot(edges[:-1],mass_Z,color='k',linewidth=2)
    # plt.plot(edges6[:-1],abs(mass_Z6),color='k',linewidth=2)
    # plt.plot([2.25,2.25],[0.1,max(APK2['Z'])],color='k',linewidth=2,label=r'Kepler Z range',alpha=0.2)
    # plt.plot([2.24,2.26],[0.1,0.1],color='k',linewidth=2,label=None,alpha=0.2)
    # plt.plot([2.24,2.26],[max(APK2['Z']),max(APK2['Z'])],color='k',linewidth=2,label=None,alpha=0.2)
    plt.xlabel(r'Mass [M$_{\odot}$]', fontsize=15)
    plt.ylabel(r'Z [kpc]', fontsize=15)
    plt.xlim(min(C6_New['mass'])-0.05,max(C6_New['mass'])+0.05)
    f.set_rasterized(True)
    # plt.title(r'Exploring Spectroscopic Peak Age Properties')
    plt.legend()
    plt.tight_layout()
    f.savefig('/home/bmr135/Dropbox/GES-K2/Ages/figure_1.eps',rasterized=True,dpi=400)
    plt.show()

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
    plt.hist(C6_New['logAge'],bins=np.linspace(8.5,10.5,40),histtype='step',normed=True,label=r'C6 Photometry',linewidth=2)
    plt.hist(RC6['logAge'],bins=np.linspace(8.5,10.5,40),histtype='step',normed=True,label=r'RAVE C6',linewidth=2,alpha=0.65)
    # plt.hist(np.log10(APK2['age']*10**9),bins=np.linspace(8.5,10.5,75),histtype='step')#,normed=True)
    plt.xlabel(r'log$_{10}$(Age)')
    cur_axes = plt.gca()
    cur_axes.axes.get_yaxis().set_ticklabels([])
    cur_axes.axes.get_yaxis().set_ticks([])
    # plt.title(r'Clump cut: R $< 10.0$')
    plt.legend()
    plt.tight_layout()
    plt.show()

    C6_New.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/C6_New',index=False)
    APK2.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/APK2',index=False)
    RC6.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/RC6',index=False)

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
