''' Comparison of spectroscopic results where the surveys have been calibrated
    to one another using spectro_photom_comp.py '''

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import mass_distr_functs as mdf
import seaborn as sns
import sys

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
plt.rcParams["font.family"] = "serif"

ext_load = '/home/bmr135/K2_Poles/Mass_Distr_In/'

''' Plot parameter distributions of originals, calibrated and calibrated to results '''
def param_distr(df,df1,df2,lab):
    '''
    Function for the plotting of the parameter distributions
    - df: Results from original set of input parameters
    - df1: Results to which the input values were calibrated to
    - df2: Results from the calibrated sample
    - lab: labels corresponding to the datasets
    '''
    fig, ((ax,ax1),(ax2,ax3),(ax4,ax5)) = plt.subplots(3,2)
    sns.despine(left=True)
    bins = [np.linspace(0,20,40),\
           np.linspace(0.75,2.5,40),\
           np.linspace(3,25,40),\
           np.linspace(2.,4,40),\
           np.linspace(-2.,0.5,40)]

    sns.distplot(df['age'], ax=ax,bins=bins[0],label=r'Original '+lab[0])
    sns.distplot(df1['age'], ax=ax, bins=bins[0],label=r'Calibrator '+lab[1])
    sns.distplot(df2['age'], ax=ax, bins=bins[0],label=r'Calibrated '+lab[2])
    ax.set_yticks([])
    # ax.set_xlim(0,20)
    ax.set_xlabel(r'Age [Gyr]')
    # ax.legend()

    sns.distplot(df['mass'], ax=ax1,bins=bins[1])
    sns.distplot(df1['mass'], ax=ax1, bins=bins[1])
    sns.distplot(df2['mass'], ax=ax1, bins=bins[1])
    ax1.set_yticks([])
    # ax1.set_xlim(0,20)
    ax1.set_xlabel(r'Mass [M$_{\odot}$]')

    sns.distplot(df['rad'], ax=ax2,bins=bins[2])
    sns.distplot(df1['rad'], ax=ax2, bins=bins[2])
    sns.distplot(df2['rad'], ax=ax2, bins=bins[2])
    ax2.set_yticks([])
    ax2.set_xlabel(r'Radius [R$_{\odot}$]')

    sns.distplot(df['logg'], ax=ax3,bins=bins[3])
    sns.distplot(df1['logg'], ax=ax3, bins=bins[3])
    sns.distplot(df2['logg'], ax=ax3, bins=bins[3])
    ax3.set_yticks([])
    ax3.set_xlabel(r'log$_{10}$(g)')

    sns.distplot(df['feh'], ax=ax4,bins=bins[4])
    sns.distplot(df1['feh'], ax=ax4, bins=bins[4])
    sns.distplot(df2['feh'], ax=ax4, bins=bins[4])
    ax4.set_yticks([])
    ax4.set_xlabel(r'[Fe/H]')

    # sns.distplot(df['Z'], ax=ax5,bins=bins,label=r'')
    # sns.distplot(df1['Z'], ax=ax5, bins=bins,label=r'')
    # sns.distplot(df2['Z'], ax=ax5, bins=bins,label=r'')
    # ax5.set_yticks([])
    # ax5.set_xlabel(r'Age [Gyr]')
    # ax5.legend()
    ax5.axis('off')
    plt.tight_layout()
    # fig.savefig('Age_Z_spec.pdf', bbox_inches='tight')
    # pdf.savefig(fig)
    # pdf.close()
    plt.show()
    # sys.exit()

''' Calibrated Files '''
RC3 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/RC3_GES/RC3_GES_07112018.in.mo',delimiter=r'\s+') # Calibrated with GES
AP3 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/AP3_GES/AP3_GES_07112018.in.mo',delimiter=r'\s+') # Calibrated with RAVE
AP6 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/AP6_RAVE/AP6_RAVE_07112018.in.mo',delimiter=r'\s+') # Calibrated with GES

''' PARAM in files '''
R3 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/RC3_GES/RC3_GES_07112018.in',delimiter=r'\s+')
A3 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/AP3_GES/AP3_GES_07112018.in',delimiter=r'\s+')
A6 = pd.read_csv('/media/bmr135/SAMSUNG/GA/K2Poles/param_outputs/Poles/AP6_RAVE/AP6_RAVE_07112018.in',delimiter=r'\s+')

RC3 = pd.merge(RC3,R3[['#Id','teff','feh','GLON','GLAT']],how='inner',on=['#Id'])
RC3.rename(columns={'GLON':'Glon','GLAT':'Glat'},inplace=True)
AP3 = pd.merge(AP3,A3[['#Id','teff','feh','GLON','GLAT']],how='inner',on=['#Id'])
AP3.rename(columns={'GLON':'Glon','GLAT':'Glat'},inplace=True)
AP6 = pd.merge(AP6,A6[['#Id','teff','feh','GLON','GLAT']],how='inner',on=['#Id'])
AP6.rename(columns={'GLON':'Glon','GLAT':'Glat'},inplace=True)

mdf.vert_dist(RC3)
mdf.vert_dist(AP3)
mdf.vert_dist(AP6)

RC3 = RC3[(RC3['logg']>-90) & (RC3['age']>-90) & (RC3['mass']>-90) & (RC3['rad']>-99.9) & (RC3['age']<19.9)]
AP3 = AP3[(AP3['logg']>-90) & (AP3['age']>-90) & (AP3['mass']>-90) & (AP3['rad']>-99.9) & (AP3['age']<19.9)]
AP6 = AP6[(AP6['logg']>-90) & (AP6['age']>-90) & (AP6['mass']>-90) & (AP6['rad']>-99.9) & (AP6['age']<19.9)]
RC3.reset_index(drop=True)
AP3.reset_index(drop=True)
AP6.reset_index(drop=True)

''' Surveys calibrated to '''
GES = pd.read_csv(ext_load+'Normal/Sep_2018/GES_070918')
GES = GES[GES['age'] < 19.9]
GES.reset_index(drop=True)
RC6 = pd.read_csv(ext_load+'Normal/Sep_2018/RC6_070918')
RC6 = RC6[RC6['age'] < 19.9]
RC6.reset_index(drop=True)


''' Original Results '''
OR3 = pd.read_csv(ext_load+'Normal/Sep_2018/RC3_070918')
OR3 = OR3[OR3['age'] < 19.9]
OR3.reset_index(drop=True)
OA3 = pd.read_csv(ext_load+'Normal/Sep_2018/AP3_070918')
OA3 = OA3[OA3['age'] < 19.9]
OA3.reset_index(drop=True)
OA6 = pd.read_csv(ext_load+'Normal/Sep_2018/AP6_070918')
OA6 = OA6[OA6['age'] < 19.9]
OA6.reset_index(drop=True)

param_distr(OR3,GES,RC3,['RC3','GES','RC3'])
param_distr(OA3,GES,AP3,['AP3','GES','AP3'])
param_distr(OA6,RC6,AP6,['AP6','RAVE','AP6'])
