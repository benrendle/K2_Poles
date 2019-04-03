''' Reading in of data for K2_seismo_comp.py '''

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import sys
from pandas import DataFrame, read_csv
from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
import K2_plot as k2p
import K2_properties as prop
import K2_constants as const
import K2_data as dat
from numbers import Number
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import random

''' Dropbox Path '''
ext_DB = '/home/bmr135/' # Work
# ext_DB = '/home/ben/'   # Laptop
''' GA directory '''
ext_GA = '/media/bmr135/SAMSUNG/' # Work
# ext_GA = '/media/ben/SAMSUNG1/' # Hard-Drive


def TRILEGAL():
    ''' Return TRILEGAL C3/C6 data on demand (Data from Leo Giradi) '''
    TRILEGAL_C3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/k1.6_K16_c3.all.out.txt',delimiter=r'\s+')
    TRILEGAL_C3['Teff'] = 10**(TRILEGAL_C3['logTe'])
    TRILEGAL_C3['g'] = 10**(TRILEGAL_C3['logg'])
    TRILEGAL_C3['L'] = 10**(TRILEGAL_C3['logL'])
    TRILEGAL_C3['radius'] = np.sqrt(TRILEGAL_C3['Mass'] / (TRILEGAL_C3['g']/const.solar_g))
    TRILEGAL_C3['JK'] = TRILEGAL_C3['Jmag'] - TRILEGAL_C3['Kmag']
    TRILEGAL_C3['Vcut'] = TRILEGAL_C3['Kmag'] + 2*(TRILEGAL_C3['JK']+0.14) + 0.382*np.exp(2*(TRILEGAL_C3['JK']-0.2))
    TRILEGAL_C3 = TRILEGAL_C3.dropna(axis=0)

    TRILEGAL_C6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/k1.6_K16_c6.all.out.txt',delimiter=r'\s+')
    TRILEGAL_C6['Teff'] = 10**(TRILEGAL_C6['logTe'])
    TRILEGAL_C6['g'] = 10**(TRILEGAL_C6['logg'])
    TRILEGAL_C6['L'] = 10**(TRILEGAL_C6['logL'])
    TRILEGAL_C6['radius'] = np.sqrt(TRILEGAL_C6['Mass'] / (TRILEGAL_C6['g']/const.solar_g))
    TRILEGAL_C6['JK'] = TRILEGAL_C6['Jmag'] - TRILEGAL_C6['Kmag']
    TRILEGAL_C6['Vcut'] = TRILEGAL_C6['Kmag'] + 2*(TRILEGAL_C6['JK']+0.14) + 0.382*np.exp(2*(TRILEGAL_C6['JK']-0.2))
    TRILEGAL_C6 = TRILEGAL_C6.dropna(axis=0)

    # print(TRILEGAL_C3.columns.values)
    # sys.exit()

    return TRILEGAL_C3, TRILEGAL_C6

# def BESANCON():
#     ''' Return BESANCON C3/C6 fields (Data from Celine Reyle) '''
#     c3 = pd.read_csv(ext_GA+'GA/K2Poles/K2c3.sim1705',delim_whitespace=True)
#     c3 = prop.galactic_coords2(c3)
#     c3['logTe'] = np.log10(c3['Teff'])
#     c3['numax'] = c3['IniMass'] * c3['Radius']**-2 * (c3['Teff']/const.solar_Teff)**-0.5 * const.solar_Numax
#     c3['dnu'] = c3['IniMass']**0.5 * c3['Radius']**-1.5 * const.solar_Dnu
#     c3['imag'] = np.nan
#     c3['Hmag'] = c3['V'] - c3['V-H']
#     c3.rename(columns={'J-K':'JK','V':'Vmag'},inplace=True)
#     c3['Kmag'] = c3['Vmag'] - c3['V-J'] - c3['JK']
#
#     c6 = pd.read_csv(ext_GA+'GA/K2Poles/K2c6.sim1705',delim_whitespace=True)
#     c6['logTe'] = np.log10(c6['Teff'])
#     c6['numax'] = c6['IniMass'] * c6['Radius']**-2 * (c6['Teff']/const.solar_Teff)**-0.5 * const.solar_Numax
#     c6['dnu'] = c6['IniMass']**0.5 * c6['Radius']**-1.5 * const.solar_Dnu
#     c6.rename(columns={'J-K':'JK','V':'Vmag'},inplace=True)
#     c6 = prop.galactic_coords2(c6)
#     c6['imag'] = np.nan
#     c6['Kmag'] = c6['Vmag'] - c6['V-J'] - c6['JK']
#     c6['Vcut'] = c6['Kmag'] + 2*(c6['JK']+0.14) + 0.382*np.exp(2*(c6['JK']-0.2))
#
#     ''' Kepler magnitude calculation: Eq. 4, Huber et al., 2016 '''
#     c3['KepMag'] = 0.314377 + 3.85667*c3['JK'] + 3.176111*c3['JK']**2 - \
#                    25.3126*c3['JK']**3 + 40.7221*c3['JK']**4 - \
#                    19.2112*c3['JK']**5 + c3['Kmag']
#     c6['KepMag'] = 0.314377 + 3.85667*c6['JK'] + 3.176111*c6['JK']**2 - \
#                25.3126*c6['JK']**3 + 40.7221*c6['JK']**4 - \
#                19.2112*c6['JK']**5 + c6['Kmag']
#
#     data = [c3,c6]
#     numax = ['numax','numax']
#     dnu = ['dnu','dnu']
#     Numax = [const.solar_Numax,const.solar_Numax]
#     Dnu = [const.solar_Dnu,const.solar_Dnu]
#     # c3,c6 = prop.lmrl_comps(data,numax,dnu,Numax,Dnu,0)
#
#     return c3, c6

def C3_cat():
    # EPIC C3 Catalogue
    C3_1 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/C3_All_EPICS1.txt')
    C3_2 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/C3_All_EPICS2.txt')
    C3 = pd.concat([C3_1,C3_2],ignore_index=True)
    # C3.to_csv('/home/bmr135/Downloads/C3',index=False)

    return C3

def C6_cat():
    # EPIC C6 Catalogue
    C6_1 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/C6_epic_search1.txt')
    C6_2 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/C6_epic_search2.txt')
    C6_3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/C6_epic_search3.txt')
    C6_4 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/C6_epic_search4.txt')
    C6 = pd.concat([C6_1,C6_2,C6_3,C6_4])
    # C6.to_csv('/home/bmr135/Downloads/C6',index=False)

    return C6

def K2_GAP():
    ''' K2 GAP Targets <- For new files '''
    # GAP3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/GAP3')
    # GAP6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/GAP6')
    # C3_flag = pd.read_csv(ext_DB+'Dropbox/K2Poles/C3_EPIC_param_flags.txt')
    # C6_flag = pd.read_csv(ext_DB+'Dropbox/K2Poles/C6_EPIC_param_flags.txt')
    # GAP3['JK'] = GAP3['Jmag'] - GAP3['Kmag']
    # GAP3['BV'] = GAP3['Bmag'] - GAP3['Vmag']
    # GAP3['sig_Teff'] = (abs(GAP3['ep_teff'])+abs(GAP3['em_teff']))/2
    # GAP3['sig_logg'] = (abs(GAP3['ep_logg'])+abs(GAP3['em_logg']))/2
    # GAP3['sig_feh'] = (abs(GAP3['ep_[Fe/H]'])+abs(GAP3['em_[Fe/H]']))/2
    # GAP3['[Fe/H]'] = GAP3['[Fe/H]'] - 0.2
    # # GAP3['[Fe/H]'] = np.random.normal(-0.294,0.305,len(GAP3)) # mu/std dev. from RAVE (22/02/2018)
    # GAP6['JK'] = GAP6['Jmag'] - GAP6['Kmag']
    # GAP6['BV'] = GAP6['Bmag'] - GAP6['Vmag']
    # GAP6['Vcut'] = GAP6['Kmag'] + 2*(GAP6['JK']+0.14) + 0.382*np.exp(2*(GAP6['JK']-0.2))
    # GAP6['sig_Teff'] = (abs(GAP6['ep_teff'])+abs(GAP6['em_teff']))/2
    # # GAP6['[Fe/H]'] = np.random.normal(-0.405,0.437,len(GAP6)) # mu/std dev. from RAVE (22/02/2018)
    # GAP6['[Fe/H]'] = GAP6['[Fe/H]'] - 0.2
    # ''' [Fe/H] uncertainty threshold using std. dev. of population '''
    # sig3 = np.std(GAP3['[Fe/H]'])
    # sig6 = np.std(GAP6['[Fe/H]'])
    #
    #
    # ''' Minimum threshold uncertainty based on distribution width '''
    # for i in range(len(GAP3['sig_Teff'])):
    #     if GAP3['sig_Teff'][i] < 100:
    #         GAP3['sig_Teff'][i] = 100
    #     if GAP3['sig_feh'][i] < sig3:
    #         GAP3['sig_feh'][i] = sig3
    # GAP6['sig_logg'] = (abs(GAP6['ep_logg'])+abs(GAP6['em_logg']))/2
    # GAP6['sig_feh'] = (abs(GAP6['ep_[Fe/H]'])+abs(GAP6['em_[Fe/H]']))/2
    # for i in range(len(GAP6['sig_Teff'])):
    #     if GAP6['sig_Teff'][i] < 100:
    #         GAP6['sig_Teff'][i] = 100
    #     if GAP6['sig_feh'][i] < sig6:
    #         GAP6['sig_feh'][i] = sig6
    #
    # GAP3 = pd.merge(GAP3,C3_flag,how='inner',on=['EPIC'])
    # GAP3 = GAP3.reset_index(drop=True)
    # GAP6 = pd.merge(GAP6,C6_flag,how='inner',on=['EPIC'])
    # GAP6 = GAP6.reset_index(drop=True)
    # GAP3.to_csv(ext_DB+'Dropbox/K2Poles/GAP3',index=False)
    # GAP6.to_csv(ext_DB+'Dropbox/K2Poles/GAP6',index=False)

    """
    Created on Sun Jul  1 18:13:39 2018
    Bolometric correcion implementation code.
    @author: Saniya Khan
    """
    # K2_camp = pd.concat([GAP3,GAP6],ignore_index=True)
    # K2_camp = K2_camp.reset_index(drop=True)
    #
    # g = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/GAP_Gaia_x_gaiadr2.csv')
    # g = g.join(K2_camp[['2MASS','EPIC','Teff','Kmag','[Fe/H]','logg']].set_index('2MASS'),on='2MASS')
    # g['Kabs'] = g['Kmag'] - 5*np.log10((g['dist_ABJ']*1000)/10)
    # g.to_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/GAP_Gaia_x_gaiadr2.csv')
    # print(g.columns.values)
    # sys.exit()
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # plt.rc('axes', labelsize=12)
    # plt.rc('xtick', labelsize=12)
    # plt.rc('ytick', labelsize=12)
    # plt.rc('legend', fontsize=12)
    # plt.rcParams['savefig.dpi'] = 300
    #
    # pd.set_option('display.max_seq_items', None)
    #
    # # data = pd.read_table('/home/bmr135/K2_Poles/Mass_Distr_In/g_COR_DR14_R7S29.txt', sep=r',', header=0)
    # g['EBV'] = 0
    # cols = ['EPIC', 'logg', '[Fe/H]', 'Teff', 'EBV']
    # g2 = g[cols]
    #
    #
    # g2.to_csv('input.sample.all', header=None, sep=r' ', index=False)
    # p = subprocess.Popen(['./bcall'])
    # p.wait()
    # p = subprocess.Popen(["mv", "output.file.all", "/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/BCs_6.csv"])
    # p.wait()
    # BC = pd.read_table('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/GAP_Gaia_BCs.csv', sep=r'\s+', header=0)
    # BC = BC[['ID', 'BC_1', 'BC_2']]
    # BC = BC.rename(columns={'ID':'EPIC','BC_1':'BC_K','BC_2':'BC_G'})
    #
    # data_BC = pd.merge(g, BC, on=['EPIC'], how='inner')
    # # data_BC = data_BC.drop(['ID'], axis=1)
    # data_BC.to_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/GAP_Gaia_BC_full.csv', index=False)
    # sys.exit()
    ''' Computation of radii from Gaia '''
    # GAP3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/GAP3')
    # GAP3 = GAP3.drop(columns=['Rgaia','radius_val','Kabs','glogg'])
    # GAP6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/GAP6')
    # GAP6 = GAP6.drop(columns=['Rgaia','radius_val','Kabs','glogg'])
    # df = pd.read_csv('/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/GAP_Gaia_BC_full.csv')
    # df['Kabs'] = df['Kmag'] - 5*np.log10((df['dist_ABJ']*1000)/10)
    # # print(df.columns.values)
    # T = (5777/df['Teff'])**4 # Temperature term
    # Mm = (4.75-df['BC_K']-df['Kabs']-df['A_K'])/2.5 # Magnitude term including bolometric correction and extinction
    # df['Rgaia'] = np.sqrt(T * 10**(Mm))
    # df = df.dropna(subset=['Rgaia'])
    # df = df.reset_index(drop=True)
    # GAP3 = pd.merge(GAP3,df[['EPIC','Rgaia','radius_val','Kabs']],on=['EPIC'])
    # GAP3['glogg'] = np.log10((const.G * GAP3['mass']*const.solar_mass)/((GAP3['Rgaia']*const.solar_radius)**2))
    # GAP6 = pd.merge(GAP6,df[['EPIC','Rgaia','radius_val','Kabs']],on=['EPIC'])
    # GAP6['glogg'] = np.log10((const.G * GAP6['mass']*const.solar_mass)/((GAP6['Rgaia']*const.solar_radius)**2))
    # # fig, ax = plt.subplots()
    # # ax.scatter(df['Rgaia'],df['radius_val'])
    # # ax.set_xlabel(r'Radius, Calc.')
    # # ax.set_ylabel(r'Radius, Gaia Catalogue')
    # # plt.show()
    # GAP3.to_csv(ext_DB+'Dropbox/K2Poles/GAP3',index=False)
    # GAP6.to_csv(ext_DB+'Dropbox/K2Poles/GAP6',index=False)
    # sys.exit()

    ''' GAP read in when no updates required '''
    GAP3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/GAP3')
    GAP6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/GAP6')

    return GAP3, GAP6

def KASOC_LC_in():
    ''' Reading in and sorting of KASOC LCs '''
    KASOC6 = pd.read_csv(ext_DB+'/home/bmr135/Drobox/K2Poles/Data0405/KASOC_EPICS')
    KASOC6 = KASOC6.convert_objects(convert_numeric=True)
    KASOC6.columns.to_series().groupby(KASOC6.dtypes).groups    # Converts EPIC IDs from objects to numbers
    KASOC_GAP6 = pd.merge(GAP6,KASOC6,how='inner',on=['EPIC'])
    return KASOC_GAP6

def Yvonne():
    ''' Yvonne Seismo C3 and C6 '''

    Yvonne_C3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Yvonne/c3-YE-results.txt')

    Yvonne_C6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Yvonne/c6-YE-results.txt')

    ''' Everest Light Curves '''
    Yvonne_EC3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Yvonne/YE-C3-Everest-results.txt')
    Yvonne_EC6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Yvonne/YE-resultsC6Everest2017-04-30.csv')

    return Yvonne_C3, Yvonne_C6, Yvonne_EC6, Yvonne_EC3

def Savita():
    ''' Savita C3/C6 '''
    Savita_C3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Savita/Results_A2Z_K2_C3_ok.txt',names=['EPIC','SDnu','e_SDnu','Snumax','e_Snumax'],delim_whitespace=True)
    Savita_C6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Savita/Results_A2Z_K2_C6_ok.txt',names=['EPIC','SDnu','e_SDnu','Snumax','e_Snumax'],delim_whitespace=True)
    Savita_EC6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Savita/A2Z_results_K2_C6_Everest_2017-04-06.csv',skiprows=1,names=['EPIC', 'Snumax', 'e_Snumax', 'SDnu', 'e_SDnu', 'Amax', 'err_Amax'])
    Savita_EC3_1 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Savita/list_K2_C3_Ok_A2Z_A2ZR_Everest.txt',names=['EPIC','SDnu','e_SDnu','Snumax','e_Snumax'],delimiter=r'\s+')
    Savita_EC3_2 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Savita/list_A2Z_C3_redo_ok_05_Everest.txt',names=['EPIC','SDnu','e_SDnu','Snumax','e_Snumax'],delimiter=r'\s+')
    Savita_EC3_3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Savita/list_K2_C3_A2Z_addon1.txt',names=['EPIC','SDnu','e_SDnu','Snumax','e_Snumax'],delimiter=r'\s+')

    Savita_EC3 = pd.concat([Savita_EC3_1,Savita_EC3_2,Savita_EC3_3],ignore_index=True)
    Savita_EC3 = Savita_EC3.drop_duplicates(subset=['EPIC'])
    Savita_EC3 = Savita_EC3.reset_index(drop=True)

    Savita_C3 = Savita_C3[Savita_C3.Snumax > -1]
    Savita_C3 = Savita_C3[Savita_C3.SDnu > -1]
    Savita_C6 = Savita_C6[Savita_C6.Snumax > -1]
    Savita_C6 = Savita_C6[Savita_C6.SDnu > -1]
    Savita_EC3 = Savita_EC3[Savita_EC3.Snumax > -1]
    Savita_EC3 = Savita_EC3[Savita_EC3.SDnu > -1]
    Savita_EC6 = Savita_EC6[Savita_EC6.Snumax > -1]
    Savita_EC6 = Savita_EC6[Savita_EC6.SDnu > -1]

    return Savita_C3, Savita_C6, Savita_EC3, Savita_EC6

def Benoit():
    ''' Benoit Seismo C3/C6 K2P2/C6 Everest '''
    Benoit_C3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Benoit/K2_fin_C3_E_corrected.txt',skiprows=3, \
                            names=['EPIC','Bnumax','e_Bnumax','BDnu','e_BDnu','A','errA'],delim_whitespace=True)

    Benoit_C6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Benoit/K2_fin_C06_01_E.txt',skiprows=3, \
                            names=['EPIC','Bnumax','e_Bnumax','BDnu','e_BDnu','A','errA'],delim_whitespace=True)
    Everest_C3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Benoit/K2_fin_C_3_E.txt',skiprows=3, \
                            names=['EPIC','Enumax','e_Enumax','EDnu','e_EDnu','A','errA'],delim_whitespace=True)
    Everest_C6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/K2_fin_C_6_E.txt',skiprows=3, \
                            names=['EPIC','Enumax','e_Enumax','EDnu','e_EDnu','A','errA'],delim_whitespace=True)

    return Benoit_C3, Benoit_C6, Everest_C3, Everest_C6

def RAVE():
    ''' RAVE C3/C6 Catalogues (Data from Marica) '''
    RAVE3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/RAVE/RAVE-K2C3_TGP_final.csv')
    RAVE3.rename(columns={'EPIC_ID':'EPIC'},inplace=True)
    RAVE3['ALPHA'] = (RAVE3['[Mg/H]']-RAVE3['[Fe/H]_RAVE'] + RAVE3['[Si/H]']-RAVE3['[Fe/H]_RAVE'])/2
    RAVE3['sig_FEH'] = abs(RAVE3['sup.3']-RAVE3['inf.3'])/2
    RAVE6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/RAVE/RAVE-K2C6_TGPfinal.csv')
    RAVE6.rename(columns={'EPIC ID':'EPIC'},inplace=True)
    RAVE6['ALPHA'] = ((RAVE6['[Mg/H]']-RAVE6['[Fe/H]_RAVE']) + (RAVE6['[Si/H]']-RAVE6['[Fe/H]_RAVE']))/2
    RAVE6['sig_FEH'] = abs(RAVE6['sup.3']-RAVE6['inf.3'])/2
    for i in range(len(RAVE3)):
        if (RAVE3['ALPHA'].iloc[i] > 1.5) | (RAVE3['ALPHA'].iloc[i] < -2.5):
            RAVE3['ALPHA'].iloc[i] = -99.9
        if abs(RAVE3['sig_FEH'].iloc[i]) > 0.5:
            RAVE3['sig_FEH'].iloc[i] = 0.2
    for i in range(len(RAVE6)):
        if (RAVE6['ALPHA'].iloc[i] > 1.5) | (RAVE6['ALPHA'].iloc[i] < -2.5):
            RAVE6['ALPHA'].iloc[i] = -99.9
        if abs(RAVE6['sig_FEH'].iloc[i]) > 0.5:
            RAVE6['sig_FEH'].iloc[i] = 0.2
    RAVE3.to_csv(ext_GA+'GA/K2Poles/RAVE3.csv',index=False)
    RAVE6.to_csv(ext_GA+'GA/K2Poles/RAVE6.csv',index=False)
    return RAVE3, RAVE6

def RAVE_merge(K2,RAVE,name):
    ''' Merge K2 data with RAVE and save out file'''
    for i in range(0,len(K2),1):
        a = RAVE[i]
        b = K2[i]
        c = pd.merge(b,a,how='inner',on=['EPIC'])
        K2[i] = c
        # c.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/'+name[i]+'_PARAM_ready.csv',index=False,sep='\t',na_rep='Inf')

    return K2

def Gaia_ESO():
        #Gaia_ESO_C3 = pd.read_csv(ext_GA+'GA/K2Poles/Gaia_ESO/GES_EPICS_Seismic1_SpecParameters1.txt',delim_whitespace=True)
        '''Paula iters'''
        # Gaia_ESO_C3 = pd.read_csv(ext_DB+'Dropbox/GES-K2/Laura_Magrini/epinarbo_ite2_new.txt')
        # Gaia_ESO_C3 = pd.read_csv(ext_DB+'Dropbox/GES-K2/param_tables/params_abund.txt',delim_whitespace=True)
        # Gaia_ESO_C3 = pd.read_csv(ext_DB+'Dropbox/GES-K2/Diane_Feuillet/epinarbo_ite2_new.txt')
        # Gaia_ESO_C3 = pd.read_csv(ext_DB+'Dropbox/GES-K2/Alvin_Gavel/lumba_ite1_uves_new.txt')
        Gaia_ESO_C3 = pd.read_csv(ext_DB+'Dropbox/GES-K2/param_tables/Final_Spec_Params_w_seismo.bin')
        alphas = pd.read_csv('/home/bmr135/Dropbox/GES-K2/param_tables/Alpha_Abundances.dat',delimiter=r'\s+')
        # Gaia_ESO_C3.rename(columns={'OBJECT':'EPIC'},inplace=True)
        Gaia_ESO_C3['EPIC'] = Gaia_ESO_C3['EPIC'].map(lambda x: x.split('_')[-1])
        Gaia_ESO_C3['EPIC'] = Gaia_ESO_C3['EPIC'].convert_objects(convert_numeric=True)
        alphas['EPIC'] = alphas['EPIC'].map(lambda x: x.split('_')[-1])
        alphas['EPIC'] = alphas['EPIC'].convert_objects(convert_numeric=True)
        alphas['ALPHA'] = alphas['MG1'] #+ alphas['SI1'] + alphas['SI2'] + alphas['CA1'] + alphas['CA2'] + alphas['TI1'] + alphas['TI2'])/7.
        Gaia_ESO_C3 = pd.merge(Gaia_ESO_C3,alphas[['EPIC','ALPHA']],how='inner',on=['EPIC'])
        Gaia_ESO_C3['sig_LOGG'] = 0
        for i in range(len(Gaia_ESO_C3)):
            Gaia_ESO_C3['sig_LOGG'].iloc[i] = random.randint(5,25)/100
            if Gaia_ESO_C3['E_TEFF'].iloc[i] < 70:
                Gaia_ESO_C3['E_TEFF'].iloc[i] = Gaia_ESO_C3['E_TEFF'].iloc[i] + 50
            if Gaia_ESO_C3['E_FEH'].iloc[i] < 0.1:
                Gaia_ESO_C3['E_FEH'].iloc[i] = Gaia_ESO_C3['E_FEH'].iloc[i] + 0.1
        # print(Gaia_ESO_C3)
        # sys.exit()
        return Gaia_ESO_C3

def GES_merge(K2,GES,name):
    ''' Merge K2 data with GES and save out file'''
    for i in range(0,len(K2),1):
        b = K2[i]
        cols_to_use = b.columns.difference(GES.columns)
        cols_to_use = cols_to_use.union(['EPIC'])
        c = pd.merge(b[cols_to_use],GES,how='inner',on=['EPIC'])
        c = c.reset_index(drop=True)
        K2[i] = c
        c.to_csv(ext_GA+'/GA/K2Poles/Gaia_ESO/'+name[i]+'_GES_match.csv',index=False,sep='\t',na_rep='Inf')

    return K2

def APOGEE():
    ''' Input for APOGEE data. C6 only as APOGEE does not survey the Southern
        hemisphere. (Data from Marica - need to apply logg cuts)'''
    # APO6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/APOGEE-K2C6.csv')
    ''' Code for new APOGEE date releases '''
    # data = fits.getdata(ext_DB+'Dropbox/K2Poles/Data0405/APOGEE_DR14Exploratory_C3GAP.fits', 1)
    # t = Table(data)
    # cols = ['NINST','STABLERV_CHI2','STABLERV_RCHI2','CHI2_THRESHOLD','STABLERV_CHI2_PROB', \
    #         'PARAM','FPARAM','PARAM_COV','FPARAM_COV','PARAMFLAG','FELEM','FELEM_ERR', \
    #         'X_H','X_H_ERR','X_M','X_M_ERR','ELEM_CHI2','ELEMFLAG','ALL_VISIT_PK', \
    #         'VISIT_PK','FPARAM_CLASS','CHI2_CLASS']
    # cols_corogee = ['PARAM_COV_2','FPARAM_COV_2','FELEM_2','FELEM_ERR_2','X_H','X_H_ERR', \
    #                 'X_M','X_M_ERR','ELEM_CHI2_2','ELEMFLAG_2']
    # t.remove_columns(cols)
    # flag = t['ASPCAPFLAG'] <= 5
    # APO3 = t[flag].to_pandas()
    # APO3.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/APOGEE_DR14_C3_180418',index=False)
    #
    # data = fits.getdata(ext_DB+'Dropbox/K2Poles/Data0405/APOGEE_DR14Exploratory_C6GAP.fits', 1)
    # data = fits.getdata(ext_DB+'Downloads/CoRoGEE_DR14.fits', 1)
    # t = Table(data)
    # print(t.info)
    # t.remove_columns(cols_corogee)
    # flag = t['ASPCAPFLAG'] <= 5
    # APO6 = t[flag].to_pandas()
    # APO6 = t.to_pandas()
    # APO6.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/CoRoGEE_DR14',index=False,na_rep='Inf')
    # sys.exit()
    ''' Input for current APOGEE data '''
    APO3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/APOGEE_DR14_C3_180418')
    APO6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/APOGEE_DR14_C6_180418')
    # print(len(APO3),len(APO6))

    return APO3, APO6

def APO_merge(K2,APO,name):
    ''' Merge K2 data with APOGEE and save out file'''
    for i in range(0,len(K2),1):
        b = K2[i]
        cols_to_use = b.columns.difference(APO.columns)
        cols_to_use = cols_to_use.union(['EPIC'])
        c = pd.merge(b[cols_to_use],APO,how='inner',on=['EPIC'])
        c = c.reset_index(drop=True)
        K2[i] = c
        c.to_csv(ext_GA+'GA/K2Poles/'+name[i]+'_APO_match.csv',index=False,na_rep='Inf')

    return K2

def LAMOST():
    ''' Read in a format LAMOST data '''
    C3 = pd.read_csv(ext_GA+'GA/K2Poles/APO_LAMOST/C3_LAMOST_output.csv',comment='#')
    C3 = C3.convert_objects(convert_numeric=True)
    C3 = C3.drop_duplicates(subset='input_id',keep='first')
    C3 = C3.dropna()
    C3 = C3.reset_index(drop=True)
    C3_list = pd.read_csv(ext_GA+'GA/K2Poles/APO_LAMOST/C3_LAMOST.csv')
    c3_comp = pd.merge(C3,C3_list,how='inner',on=['RA'])
    c3_comp = c3_comp.reset_index(drop=True)

    C6 = pd.read_csv(ext_GA+'GA/K2Poles/APO_LAMOST/C6_LAMOST_output.csv',comment='#')
    C6 = C6.convert_objects(convert_numeric=True)
    C6 = C6.dropna()
    C6 = C6.drop_duplicates(subset='input_id',keep='first')
    C6 = C6.reset_index(drop=True)
    C6_list = pd.read_csv(ext_GA+'GA/K2Poles/APO_LAMOST/C6_LAMOST.csv')
    c6_comp = pd.merge(C6,C6_list,how='inner',on=['RA'])
    c6_comp = c6_comp.reset_index(drop=True)

    return c3_comp, c6_comp

def LAMOST_merge(K2,LAMOST,name):
    ''' Merge K2 data with LAMOST and save out file'''
    for i in range(0,len(K2),1):
        b = K2[i]
        cols_to_use = b.columns.difference(LAMOST.columns)
        cols_to_use = cols_to_use.union(['EPIC'])
        c = pd.merge(b[cols_to_use],LAMOST,how='inner',on=['EPIC'])
        c = c.reset_index(drop=True)
        K2[i] = c
        c.to_csv(ext_GA+'GA/K2Poles/'+name[i]+'_LAMOST_match.csv',index=False,na_rep='Inf')

    return K2[0],K2[1],K2[2],K2[3]

def occurrence():
    ''' Read in occurence data '''
    oc = pd.read_csv(ext_GA+'K2P2_MNL_K2GAP_C6/EPIC_occurrence')

    return oc

def n_epics(df,oc):
    ''' Assign occurrence number to each EPIC value '''
    df = pd.merge(df,oc,how='inner',on=['EPIC'])
    df = df.reset_index(drop=True)

    return df
