#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import matplotlib.gridspec as gridspec

if __name__ == "__main__":
    ### Read in the APOKASC data ...
    context_csv_file  ='/home/bmr135/K2_Poles/Mass_Distr_In/Normal/Nov_2018/APOKASC_281118'
    df = pd.read_csv(context_csv_file)
    df = df[df.mass > 0] # Kill stars with no mass

    ### Read in the new sample
    sample_csv_file = '/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/SM_Gaia_BC_full.csv'
    sample = pd.read_csv(sample_csv_file)
    sample6_csv_file = '/home/bmr135/K2_Poles/Mass_Distr_In/Gaia/AS_Gaia_BC_APO.csv'
    sample6 = pd.read_csv(sample6_csv_file)
    sample6 = sample6[(sample6['alpha'] > 0.1) & (sample6.age < 7)]
    sample6 = sample6.reset_index(drop=True)

    ### Set up the plot
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=(2,4))
    ax = fig.add_subplot(gs[0], facecolor='white')
    bx = fig.add_subplot(gs[1], facecolor='white')
    matplotlib.rcParams.update({'font.size': 20})

    ### plot the Context stars
    ax.scatter(df.numax, df.Dnu, label=r'APOKASC', color='gray', alpha=0.2)
    CS = ax.scatter(sample6.nmx, sample6.dnu, c=sample6.mass, cmap='copper', \
                    label='K2 YAR SAMPLE', vmax=1.5)
    ax.set_xlabel(r'$\mathrm{\nu_{max}}$', fontsize=20)
    ax.set_ylabel(r'$\mathrm{\Delta \nu}$', fontsize=20)
    cbar = fig.colorbar(CS, extend='both')
    cbar.ax.set_ylabel('Mass [M$_{\odot}$]', fontsize=20, rotation=270, labelpad=20)
    ax.set_xlim([0,220])
    ax.set_ylim([0,20])

    ### Set up a reference line
    numaxx = np.linspace(0,288)
    alpha = 0.22
    beta = 0.8
    fit = alpha * numaxx **beta
    fitt = alpha * df.numax **beta
    ax.plot(numaxx, fit, 'r--')

    ### Plot the context stars O-C
    bx.scatter(df.numax, df.Dnu - fitt, color='gray',alpha=0.2)#, c=df.mass, cmap='copper', vmax=2,alpha=0.25)
    bx.set_xlim([0,220])
    fig.tight_layout()

    ### Plot the sample stars in the O-C
    oc = sample.Dnu - alpha * sample.numax**beta
    oc6 = sample6.dnu - alpha * sample6.nmx**beta
    bx.scatter(sample6.nmx, oc6, label='SM', alpha=0.73, s=40, c=sample6.mass, cmap='copper', vmax=2,)
    # bx.scatter(sample6.nmx, oc6, edgecolors='m', \
    #            label='AS', alpha=0.73, s=40)
    bx.set_ylim([-2,1.5])
    bx.set_xlabel(r'$\mathrm{\nu_{max}}$', fontsize=20)
    bx.set_ylabel(r'$\mathrm{\Delta \nu}$ - fit', fontsize=20)

    # fig.savefig('context_mass.png')

    plt.show()
