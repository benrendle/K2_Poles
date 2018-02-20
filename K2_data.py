''' Reading in of data for K2_seismo_comp.py '''

import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import sys
from pandas import DataFrame, read_csv
from astropy.io import fits
from astropy.table import Table
import numpy as np
import K2_plot as k2p
import K2_properties as prop
import K2_constants as const
import K2_data as dat
from numbers import Number

''' Dropbox Path '''
ext_DB = '/home/bmr135/' # Work
# ext_DB = '/home/ben/'   # Laptop
''' GA directory '''
ext_GA = '/home/bmr135/' # Work
# ext_GA = '/media/ben/SAMSUNG/' # Hard-Drive


def TRILEGAL():
    ''' Return TRILEGAL C3/C6 data on demand (Data from Leo Giradi) '''
    TRILEGAL_C3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/k1.6_K16_c3.all.out.txt',delimiter=r'\s+')
    TRILEGAL_C3['Teff'] = 10**(TRILEGAL_C3['logTe'])
    TRILEGAL_C3['g'] = 10**(TRILEGAL_C3['logg'])
    TRILEGAL_C3['L'] = 10**(TRILEGAL_C3['logL'])
    TRILEGAL_C3['radius'] = np.sqrt(TRILEGAL_C3['Mass'] / (TRILEGAL_C3['g']/const.solar_g))
    TRILEGAL_C3['JK'] = TRILEGAL_C3['Jmag'] - TRILEGAL_C3['Kmag']
    TRILEGAL_C3 = TRILEGAL_C3.dropna(axis=0)

    TRILEGAL_C6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/k1.6_K16_c6.all.out.txt',delimiter=r'\s+')
    TRILEGAL_C6['Teff'] = 10**(TRILEGAL_C6['logTe'])
    TRILEGAL_C6['g'] = 10**(TRILEGAL_C6['logg'])
    TRILEGAL_C6['L'] = 10**(TRILEGAL_C6['logL'])
    TRILEGAL_C6['radius'] = np.sqrt(TRILEGAL_C6['Mass'] / (TRILEGAL_C6['g']/const.solar_g))
    TRILEGAL_C6['JK'] = TRILEGAL_C6['Jmag'] - TRILEGAL_C6['Kmag']
    TRILEGAL_C6['Vcut'] = TRILEGAL_C6['Kmag'] + 2*(TRILEGAL_C6['JK']+0.14) + 0.382*np.exp(2*(TRILEGAL_C6['JK']-0.2))
    TRILEGAL_C6 = TRILEGAL_C6.dropna(axis=0)

    return TRILEGAL_C3, TRILEGAL_C6

def BESANCON():
    ''' Return BESANCON C3/C6 fields (Data from Celine Reyle) '''
    c3 = pd.read_csv(ext_GA+'GA/K2Poles/K2c3.sim1705',delim_whitespace=True)
    c3 = prop.galactic_coords2(c3)
    c3['logTe'] = np.log10(c3['Teff'])
    c3['numax'] = c3['IniMass'] * c3['Radius']**-2 * (c3['Teff']/const.solar_Teff)**-0.5 * const.solar_Numax
    c3['dnu'] = c3['IniMass']**0.5 * c3['Radius']**-1.5 * const.solar_Dnu
    c3['imag'] = np.nan
    c3['Hmag'] = c3['V'] - c3['V-H']
    c3.rename(columns={'J-K':'JK','V':'Vmag'},inplace=True)
    c3['Kmag'] = c3['Vmag'] - c3['V-J'] - c3['JK']

    c6 = pd.read_csv(ext_GA+'GA/K2Poles/K2c6.sim1705',delim_whitespace=True)
    c6['logTe'] = np.log10(c6['Teff'])
    c6['numax'] = c6['IniMass'] * c6['Radius']**-2 * (c6['Teff']/const.solar_Teff)**-0.5 * const.solar_Numax
    c6['dnu'] = c6['IniMass']**0.5 * c6['Radius']**-1.5 * const.solar_Dnu
    c6.rename(columns={'J-K':'JK','V':'Vmag'},inplace=True)
    c6 = prop.galactic_coords2(c6)
    c6['imag'] = np.nan
    c6['Kmag'] = c6['Vmag'] - c6['V-J'] - c6['JK']
    c6['Vcut'] = c6['Kmag'] + 2*(c6['JK']+0.14) + 0.382*np.exp(2*(c6['JK']-0.2))

    ''' Kepler magnitude calculation: Eq. 4, Huber et al., 2016 '''
    c3['KepMag'] = 0.314377 + 3.85667*c3['JK'] + 3.176111*c3['JK']**2 - \
                   25.3126*c3['JK']**3 + 40.7221*c3['JK']**4 - \
                   19.2112*c3['JK']**5 + c3['Kmag']
    c6['KepMag'] = 0.314377 + 3.85667*c6['JK'] + 3.176111*c6['JK']**2 - \
               25.3126*c6['JK']**3 + 40.7221*c6['JK']**4 - \
               19.2112*c6['JK']**5 + c6['Kmag']

    data = [c3,c6]
    numax = ['numax','numax']
    dnu = ['dnu','dnu']
    Numax = [const.solar_Numax,const.solar_Numax]
    Dnu = [const.solar_Dnu,const.solar_Dnu]
    c3,c6 = prop.lmrl_comps(data,numax,dnu,Numax,Dnu,0)

    return c3, c6

def C3_cat():
    # EPIC C3 Catalogue
    C3_1 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/C3_All_EPICS1.txt')
    C3_2 = pd.read_csv(ext_DB'Dropbox/K2Poles/Data0405/C3_All_EPICS2.txt')
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
    # K2 GAP Targets
    GAP3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/C3_epic_search.txt')
    GAP6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/C6_epic_search.txt')
    C3_flag = pd.read_csv(ext_DB+'Dropbox/K2Poles/C3_EPIC_param_flags.txt')
    C6_flag = pd.read_csv(ext_DB+'Dropbox/K2Poles/C6_EPIC_param_flags.txt')
    GAP3['JK'] = GAP3['Jmag'] - GAP3['Kmag']
    GAP6['JK'] = GAP6['Jmag'] - GAP6['Kmag']
    GAP6['Vcut'] = GAP6['Kmag'] + 2*(GAP6['JK']+0.14) + 0.382*np.exp(2*(GAP6['JK']-0.2))
    GAP3['sig_Teff'] = (abs(GAP3['ep_teff'])+abs(GAP3['em_teff']))/2
    for i in range(len(GAP3['sig_Teff'])):
        if GAP3['sig_Teff'][i] < 100:
            GAP3['sig_Teff'][i] = 100
    GAP3['sig_logg'] = (abs(GAP3['ep_logg'])+abs(GAP3['em_logg']))/2
    GAP3['sig_feh'] = (abs(GAP3['ep_[Fe/H]'])+abs(GAP3['em_[Fe/H]']))/2
    GAP6['sig_Teff'] = (abs(GAP6['ep_teff'])+abs(GAP6['em_teff']))/2
    for i in range(len(GAP6['sig_Teff'])):
        if GAP6['sig_Teff'][i] < 100:
            GAP6['sig_Teff'][i] = 100
    GAP6['sig_logg'] = (abs(GAP6['ep_logg'])+abs(GAP6['em_logg']))/2
    GAP6['sig_feh'] = (abs(GAP6['ep_[Fe/H]'])+abs(GAP6['em_[Fe/H]']))/2

    GAP3 = pd.merge(GAP3,C3_flag,how='inner',on=['EPIC'])
    GAP3 = GAP3.reset_index(drop=True)
    GAP6 = pd.merge(GAP6,C6_flag,how='inner',on=['EPIC'])
    GAP6 = GAP6.reset_index(drop=True)
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
    Benoit_C3 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/Benoit/K2_fin_C3_E.txt',skiprows=3, \
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
    RAVE6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/RAVE/RAVE-K2C6_TGPfinal.csv')
    RAVE6.rename(columns={'EPIC ID':'EPIC'},inplace=True)
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
        c.to_csv(ext_DB+'Dropbox/K2Poles/Data0405/'+name[i]+'_PARAM_ready.csv',index=False,sep='\t',na_rep='Inf')

    return K2

def Gaia_ESO():
        #Gaia_ESO_C3 = pd.read_csv(ext_GA+'GA/K2Poles/Gaia_ESO/GES_EPICS_Seismic1_SpecParameters1.txt',delim_whitespace=True)
        '''Paula iters'''
        # Gaia_ESO_C3 = pd.read_csv(ext_DB+'Dropbox/GES-K2/Laura_Magrini/epinarbo_ite2_new.txt')
        Gaia_ESO_C3 = pd.read_csv(ext_DB+'Dropbox/GES-K2/param_tables/params_abund.txt',delim_whitespace=True)
        # Gaia_ESO_C3 = pd.read_csv(ext_DB+'Dropbox/GES-K2/Diane_Feuillet/epinarbo_ite2_new.txt')
        # Gaia_ESO_C3 = pd.read_csv(ext_DB+'Dropbox/GES-K2/Alvin_Gavel/lumba_ite1_uves_new.txt')
        Gaia_ESO_C3.rename(columns={'OBJECT':'EPIC'},inplace=True)
        Gaia_ESO_C3['EPIC'] = Gaia_ESO_C3['EPIC'].map(lambda x: x.split('_')[-1])
        Gaia_ESO_C3['EPIC'] = Gaia_ESO_C3['EPIC'].convert_objects(convert_numeric=True)
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
    APO6 = pd.read_csv(ext_DB+'Dropbox/K2Poles/Data0405/APOGEE-K2C6.csv')
    APO6.rename(columns={'EPIC ID':'EPIC'},inplace=True)

    return APO6

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
    oc = pd.read_csv(ext_DB+'K2P2_MNL_K2GAP_C6/EPIC_occurrence')

    return oc

def n_epics(df,oc):
    ''' Assign occurrence number to each EPIC value '''
    df = pd.merge(df,oc,how='inner',on=['EPIC'])
    df = df.reset_index(drop=True)

    return df
