# Comparison of Yvonne's and Savita's C3 numax and dnu values
# Plots styled with the help of GRD
# Last Modified: 1/03/17, bmr135 Rendle

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from astropy.io import fits
from astropy.table import Table
from scipy.stats import norm
import colormaps
import K2_properties as prop
import K2_data as dat
import sys

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams["font.family"] = "serif"

''' Import and organise data '''

TRILEGAL_C3, TRILEGAL_C6 = dat.TRILEGAL()
C3 = dat.C3_cat()
C6 = dat.C6_cat()
GAP3, GAP6 = dat.K2_GAP()
Yvonne_C3, Yvonne_C6, Yvonne_EC6, Yvonne_EC3 = dat.Yvonne()
Savita_C3, Savita_C6, Savita_EC3, Savita_EC6 = dat.Savita()
Benoit_C3, Benoit_C6, Everest_C3, Everest_C6 = dat.Benoit()

YC3 = pd.merge(Yvonne_C3,GAP3,how='inner',on=['EPIC'])
YC3 = YC3.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
YC6 = pd.merge(Yvonne_C6,GAP6,how='inner',on=['EPIC'])
YC6 = YC6.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
SC3 = pd.merge(Savita_C3,GAP3,how='inner',on=['EPIC'])
SC3 = SC3.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
SC6 = pd.merge(Savita_C6,GAP6,how='inner',on=['EPIC'])
SC6 = SC6.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
BC3 = pd.merge(Benoit_C3,GAP3,how='inner',on=['EPIC'])
BC3 = BC3.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
BC6 = pd.merge(Benoit_C6,GAP6,how='inner',on=['EPIC'])
BC6 = BC6.drop_duplicates(subset=['EPIC']).reset_index(drop=True)

''' Number of stars unique ot each pipeline '''
# prop.individ(YC3,BC3,SC3)
# prop.individ(YC6,BC6,SC6)

''' Asteroseismic logg '''
YC3['slogg'] = np.log10(27400*(YC3['nmx']/3135)*np.sqrt(YC3['Teff']/5777))
YC6['slogg'] = np.log10(27400*(YC6['nmx']/3135)*np.sqrt(YC6['Teff']/5777))
SC3['slogg'] = np.log10(27400*(SC3['Snumax']/3097.33)*np.sqrt(SC3['Teff']/5777))
SC6['slogg'] = np.log10(27400*(SC6['Snumax']/3097.33)*np.sqrt(SC6['Teff']/5777))
BC3['slogg'] = np.log10(27400*(BC3['Bnumax']/3104)*np.sqrt(BC3['Teff']/5777))
BC6['slogg'] = np.log10(27400*(BC6['Bnumax']/3104)*np.sqrt(BC6['Teff']/5777))

''' Combining pipelines for overlaps '''
YS_C3 = pd.merge(SC3,YC3,how='inner',on=['EPIC'])
YS_C3 = YS_C3.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
YS_C6 = pd.merge(SC6,YC6,how='inner',on=['EPIC'])
YS_C6 = YS_C6.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
YB_C3 = pd.merge(BC3,YC3,how='inner',on=['EPIC'])
YB_C3 = YB_C3.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
YB_C6 = pd.merge(BC6,YC6,how='inner',on=['EPIC'])
YB_C6 = YB_C6.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
BS_C3 = pd.merge(BC3,SC3,how='inner',on=['EPIC'])
BS_C3 = BS_C3.drop_duplicates(subset=['EPIC']).reset_index(drop=True)
BS_C6 = pd.merge(BC6,SC6,how='inner',on=['EPIC'])
BS_C6 = BS_C6.drop_duplicates(subset=['EPIC']).reset_index(drop=True)


''' numax histrograms - C3 and then C6 '''
# fig, ax = plt.subplots()
# bins = np.linspace(10,280,100)
# ax.hist(YC3['nmx'],bins,histtype='step',linewidth=2,label=r'BHM')
# ax.hist(BC3['Bnumax'],bins,histtype='step',linewidth=2,label=r'COR')
# ax.hist(SC3['Snumax'],bins,histtype='step',linewidth=2,label=r'A2Z')
# ax.set_xlim(10,280)
# ax.set_xlabel(r'$\nu_{\rm{max}}$ [$\mu$Hz]',fontsize=15)
# ax.legend()
# plt.show()
#
# fig, ax = plt.subplots()
# bins = np.linspace(10,280,100)
# ax.hist(YC6['nmx'],bins,histtype='step',linewidth=2,label=r'BHM')
# ax.hist(BC6['Bnumax'],bins,histtype='step',linewidth=2,label=r'COR')
# ax.hist(SC6['Snumax'],bins,histtype='step',linewidth=2,label=r'A2Z')
# ax.set_xlim(10,280)
# ax.set_xlabel(r'$\nu_{\rm{max}}$ [$\mu$Hz]',fontsize=15)
# ax.legend()
# plt.show()


''' dnu histograms - C3 then C6 '''
# fig, ax = plt.subplots()
# bins = np.linspace(0,20,100)
# ax.hist(YC3['dnu'],bins,histtype='step',linewidth=2,label=r'BHM')
# ax.hist(BC3['BDnu'],bins,histtype='step',linewidth=2,label=r'COR')
# ax.hist(SC3['SDnu'],bins,histtype='step',linewidth=2,label=r'A2Z')
# ax.set_xlim(0,20)
# ax.set_xlabel(r'$\Delta\nu$ [$\mu$Hz]',fontsize=15)
# ax.legend()
# plt.show()
#
# fig, ax = plt.subplots()
# bins = np.linspace(0,20,100)
# ax.hist(YC6['dnu'],bins,histtype='step',linewidth=2,label=r'BHM')
# ax.hist(BC6['BDnu'],bins,histtype='step',linewidth=2,label=r'COR')
# ax.hist(SC6['SDnu'],bins,histtype='step',linewidth=2,label=r'A2Z')
# ax.set_xlim(0,20)
# ax.set_xlabel(r'$\Delta\nu$ [$\mu$Hz]',fontsize=15)
# ax.legend()
# plt.show()


''' C3 HR/CM Diagrams - Combined on a single plot '''
# plt.figure()
# plt.subplot(3,2,1)
# # plt.scatter(TRILEGAL_C3['JK'],TRILEGAL_C3['Vcut'],color='0.75',label=r'C3 TRI')
# plt.scatter(GAP3['JK'],GAP3['Vmag'],color='grey',alpha=0.25,label=r'K2 GAP')
# plt.scatter(YC3['JK'],YC3['Vmag'],label=r'Yvonne C3',c=YC3['nmx'],cmap=colormaps.parula)
# plt.plot([0.5,0.5],[max(GAP3['Vmag'])+1,min(GAP3['Vmag'])-1],color='k', \
#         linestyle='--',linewidth=2)
# plt.plot([0.5,1.5],[9.0,9.0],color='k',linestyle='--',linewidth=2)
# plt.plot([0.5,1.5],[15.0,15.0],color='k',linestyle='--',linewidth=2)
# cbar = plt.colorbar()
# cbar.set_label(r'$\nu_{\rm{max}}$ [$\mu$Hz]', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'V',fontsize=20, labelpad=20)
# plt.ylim(max(GAP3['Vmag'])+1,min(GAP3['Vmag'])-1)
# plt.xlim(0.25,1.5)
# plt.xlabel(r'J - K',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,2)
# plt.scatter(GAP3['Teff'],GAP3['logg'],color='grey',alpha=0.25,label=r'C3 GAP')
# plt.scatter(YC3['Teff'],YC3['slogg'],label=r'Yvonne C3',c=YC3['JK'],cmap=colormaps.parula)
# cbar = plt.colorbar()
# cbar.set_label(r'J - K', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'$log(g)$',fontsize=20, labelpad=20)
# plt.ylim(5,0)
# plt.xlim(6000,3500)
# plt.xlabel(r'T$_{\rm{eff}}$',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,3)
# # plt.scatter(TRILEGAL_C3['JK'],TRILEGAL_C3['Vcut'],color='0.75',label=r'C3 TRI')
# plt.scatter(GAP3['JK'],GAP3['Vmag'],color='grey',alpha=0.25,label=r'K2 GAP')
# plt.scatter(SC3['JK'],SC3['Vmag'],label=r'Savita C3',c=SC3['Snumax'],cmap=colormaps.parula)
# plt.plot([0.5,0.5],[max(GAP3['Vmag'])+1,min(GAP3['Vmag'])-1],color='k', \
#         linestyle='--',linewidth=2)
# plt.plot([0.5,1.5],[9.0,9.0],color='k',linestyle='--',linewidth=2)
# plt.plot([0.5,1.5],[15.0,15.0],color='k',linestyle='--',linewidth=2)
# cbar = plt.colorbar()
# cbar.set_label(r'$\nu_{\rm{max}}$ [$\mu$Hz]', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'V',fontsize=20, labelpad=20)
# plt.ylim(max(GAP3['Vmag'])+1,min(GAP3['Vmag'])-1)
# plt.xlim(0.25,1.5)
# plt.xlabel(r'J - K',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,4)
# plt.scatter(GAP3['Teff'],GAP3['logg'],color='grey',alpha=0.25,label=r'C3 GAP')
# plt.scatter(SC3['Teff'],SC3['slogg'],label=r'Savita C3',c=SC3['JK'],cmap=colormaps.parula)
# cbar = plt.colorbar()
# cbar.set_label(r'J - K', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'$log(g)$',fontsize=20, labelpad=20)
# plt.ylim(5,0)
# plt.xlim(6000,3500)
# plt.xlabel(r'T$_{\rm{eff}}$',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,5)
# # plt.scatter(TRILEGAL_C3['JK'],TRILEGAL_C3['Vcut'],color='0.75',label=r'C3 TRI')
# plt.scatter(GAP3['JK'],GAP3['Vmag'],color='grey',alpha=0.25,label=r'K2 GAP')
# plt.scatter(BC3['JK'],BC3['Vmag'],label=r'Benoit C3',c=BC3['Bnumax'],cmap=colormaps.parula)
# plt.plot([0.5,0.5],[max(GAP3['Vmag'])+1,min(GAP3['Vmag'])-1],color='k', \
#         linestyle='--',linewidth=2)
# plt.plot([0.5,1.5],[9.0,9.0],color='k',linestyle='--',linewidth=2)
# plt.plot([0.5,1.5],[15.0,15.0],color='k',linestyle='--',linewidth=2)
# cbar = plt.colorbar()
# cbar.set_label(r'$\nu_{\rm{max}}$ [$\mu$Hz]', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'V',fontsize=20, labelpad=20)
# plt.ylim(max(GAP3['Vmag'])+1,min(GAP3['Vmag'])-1)
# plt.xlim(0.25,1.5)
# plt.xlabel(r'J - K',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,6)
# plt.scatter(GAP3['Teff'],GAP3['logg'],color='grey',alpha=0.25,label=r'C3 GAP')
# plt.scatter(BC3['Teff'],BC3['slogg'],label=r'Benoit C3',c=BC3['JK'],cmap=colormaps.parula)
# cbar = plt.colorbar()
# cbar.set_label(r'J - K', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'$log(g)$',fontsize=20, labelpad=20)
# plt.ylim(5,0)
# plt.xlim(6000,3500)
# plt.xlabel(r'T$_{\rm{eff}}$',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.tight_layout(pad=0.2, w_pad=0.125, h_pad=0.15)
# # plt.show()
# # sys.exit()

''' C6 CMD/HRD - combined plot '''
# plt.figure()
# plt.subplot(3,2,1)
# # plt.scatter(TRILEGAL_C6['JK'],TRILEGAL_C6['Vcut'],color='0.75',label=r'C6 TRI')
# plt.scatter(GAP6['JK'],GAP6['Vmag'],color='grey',alpha=0.25,label=r'K2 GAP')
# plt.scatter(YC6['JK'],YC6['Vmag'],label=r'Yvonne C6',c=YC6['nmx'],cmap=colormaps.parula)
# plt.plot([0.5,0.5],[max(GAP6['Vmag'])+1,min(GAP6['Vmag'])-1],color='k', \
#         linestyle='--',linewidth=2)
# plt.plot([0.5,2.0],[9.0,9.0],color='k',linestyle='--',linewidth=2)
# plt.plot([0.5,2.0],[15.0,15.0],color='k',linestyle='--',linewidth=2)
# plt.xlim(-0.35,2.0)
# cbar = plt.colorbar()
# cbar.set_label(r'$\nu_{\rm{max}}$ [$\mu$Hz]', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'V',fontsize=20, labelpad=20)
# plt.ylim(max(GAP6['Vmag'])+1,min(GAP6['Vmag'])-1)
# plt.xlabel(r'J - K',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,2)
# plt.scatter(GAP6['Teff'],GAP6['logg'],color='grey',alpha=0.25,label=r'K2 GAP')
# plt.scatter(YC6['Teff'],YC6['slogg'],label=r'Yvonne C6',c=YC6['JK'],cmap=colormaps.parula)
# cbar = plt.colorbar()
# cbar.set_label(r'J - K', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'$log(g)$',fontsize=20, labelpad=20)
# plt.ylim(5,0)
# plt.xlim(6000,3500)
# plt.xlabel(r'T$_{\rm{eff}}$',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,3)
# # plt.scatter(TRILEGAL_C6['JK'],TRILEGAL_C6['Vcut'],color='0.75',label=r'C6 TRI')
# plt.scatter(GAP6['JK'],GAP6['Vmag'],color='grey',alpha=0.25,label=r'K2 GAP')
# plt.scatter(SC6['JK'],SC6['Vmag'],label=r'Savita C6',c=SC6['Snumax'],cmap=colormaps.parula)
# plt.plot([0.5,0.5],[max(GAP6['Vmag'])+1,min(GAP6['Vmag'])-1],color='k', \
#         linestyle='--',linewidth=2)
# plt.plot([0.5,2.0],[9.0,9.0],color='k',linestyle='--',linewidth=2)
# plt.plot([0.5,2.0],[15.0,15.0],color='k',linestyle='--',linewidth=2)
# plt.xlim(-0.35,2.0)
# cbar = plt.colorbar()
# cbar.set_label(r'$\nu_{\rm{max}}$ [$\mu$Hz]', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'V',fontsize=20, labelpad=20)
# plt.ylim(max(GAP6['Vmag'])+1,min(GAP6['Vmag'])-1)
# plt.xlabel(r'J - K',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,4)
# plt.scatter(GAP6['Teff'],GAP6['logg'],color='grey',alpha=0.25,label=r'C3 GAP')
# plt.scatter(SC6['Teff'],SC6['slogg'],label=r'Savita C6',c=SC6['JK'],cmap=colormaps.parula)
# cbar = plt.colorbar()
# cbar.set_label(r'J - K', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'$log(g)$',fontsize=20, labelpad=20)
# plt.xlabel(r'T$_{\rm{eff}}$',fontsize=20, labelpad=10)
# plt.ylim(5,0)
# plt.xlim(6000,3500)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,5)
# # plt.scatter(TRILEGAL_C6['JK'],TRILEGAL_C6['Vcut'],color='0.75',label=r'C6 TRI')
# plt.scatter(GAP6['JK'],GAP6['Vmag'],color='grey',alpha=0.25,label=r'K2 GAP')
# plt.scatter(BC6['JK'],BC6['Vmag'],label=r'Benoit C6',c=BC6['Bnumax'],cmap=colormaps.parula)
# plt.plot([0.5,0.5],[max(GAP6['Vmag'])+1,min(GAP6['Vmag'])-1],color='k', \
#         linestyle='--',linewidth=2)
# plt.plot([0.5,2.0],[9.0,9.0],color='k',linestyle='--',linewidth=2)
# plt.plot([0.5,2.0],[15.0,15.0],color='k',linestyle='--',linewidth=2)
# plt.xlim(-0.35,2.0)
# cbar = plt.colorbar()
# cbar.set_label(r'$\nu_{\rm{max}}$ [$\mu$Hz]', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'V',fontsize=20, labelpad=20)
# plt.ylim(max(GAP6['Vmag'])+1,min(GAP6['Vmag'])-1)
# plt.xlabel(r'J - K',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.subplot(3,2,6)
# plt.scatter(GAP6['Teff'],GAP6['logg'],color='grey',alpha=0.25,label=r'C3 GAP')
# plt.scatter(BC6['Teff'],BC6['slogg'],label=r'Benoit C6',c=BC6['JK'],cmap=colormaps.parula)
# cbar = plt.colorbar()
# cbar.set_label(r'J - K', rotation=270, fontsize=20, labelpad=25)
# cbar.ax.tick_params(labelsize=20)
# plt.ylabel(r'$log(g)$',fontsize=20, labelpad=20)
# plt.ylim(5,0)
# plt.xlim(6000,3500)
# plt.xlabel(r'T$_{\rm{eff}}$',fontsize=20, labelpad=10)
# plt.tick_params(labelsize=15)
# plt.legend(prop={'size':10})
#
# plt.tight_layout(pad=0.2, w_pad=0.125, h_pad=0.15)
# plt.show()
# sys.exit()

''' Difference plots dnu --> Guy style '''
def comp_seismic_plot(df,comp_err,params,uncerts,labels,seismo):
    '''
    Plot to compare different asteroseimic pipelines. Scatter plot of the
    difference between the the two pipelines in terms of the compound uncertainty
    and plotted against one set of the seismic parameter. Histogram of the
    distribution is included, as well as makrings of the mean and std. dev. A
    zero point line is included for guidance. Any values with a sigma difference
    of greater than 10 are removed as the sources are likely outliers.

    - Inputs:
        - Combined dataset to compare (df)
        - Compound uncertainty type: dnu or numax (comp_err)
        - Seismic constraint uncertainty (uncerts, [])
        - Seismic constraint value (params, [])
        - Plot labels -> Names of pipelines/authors (labels, [])
        - Latex description of the input seismic parameter (seismo)

    - Outputs:
        - Plot.
        - Could be modified to save or show accordingly.

    '''
    df[comp_err] = np.sqrt(df[uncerts[0]]**2+df[uncerts[1]]**2)
    df['Diff'] = ((df[params[0]]-df[params[1]])/df[comp_err])
    df = df[abs(df['Diff']) < 10]
    mu = np.mean(df['Diff'])
    std = np.std(df['Diff'])

    fig = plt.figure(figsize=(8,4.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4,1])

    ax1 = fig.add_subplot(gs[0])
    ax1.scatter(df[params[0]], df['Diff'], s=25, facecolors='none', edgecolors='b', label='__nolegend__',alpha=0.5)
    ax1.set_ylabel(r'('+labels[0]+r' - '+labels[1]+r')/$\sigma_{comb.}$ ',fontsize=15,labelpad=10)
    ax1.set_xlabel(labels[0]+r', '+seismo+r' [$\mu$Hz]',fontsize=15, labelpad=10)
    ax1.plot([0,(max(df[params[0]])+1)],[mu,mu],color='k',linewidth=2,label=r'$\mu$: %.5s'%(mu))
    # ax1.plot([0,(max(df[params[0]])+1)],[mu+std,mu+std],color='r',linewidth=2,alpha=0.25,label=r'$\sigma$: %.5s'%(std))
    # ax1.plot([0,(max(df[params[0]])+1)],[mu-std,mu-std],color='r',linewidth=2,alpha=0.25,label='__nolegend__')
    ax1.fill_between([0,(max(df[params[0]])+1)],(mu-std),(mu+std),facecolor='grey', alpha=0.25,label=r'$\sigma$: %.5s'%(std))
    ax1.tick_params(labelsize=10)
    ax1.set_xlim(0,(max(df[params[0]])+1))
    ax1.set_ylim([-5,5])
    # ax1.plot([0,(max(df[params[0]])+1)],[0,0], 'k--',label='__nolegend__')
    ax1.legend(prop={'size':10})

    ax2 = fig.add_subplot(gs[1], sharey=ax1)
    n,bins,patches = ax2.hist(df['Diff'], 50, normed=1, orientation='horizontal', \
             facecolor='b',alpha=0.5, label='__nolegend__')
    x = np.linspace(min(df['Diff']),max(df['Diff']),100)
    y = mlab.normpdf(x, mu, std)
    print(mu, std)
    ref = mlab.normpdf(x, 0.0, 1.0)
    ax2.plot(y,x,'r--',linewidth=2, label='Fit')
    ax2.plot(ref,x,'k--',linewidth=2,label='N(0,1)')
    ax2.set_ylim([-5,5])
    ax2.set_xticks([])
    ax2.set_xlabel('PDF',fontsize=15)
    ax2.legend(prop={'size':10})
    ax2.tick_params(labelsize=10)
    plt.setp(ax2.get_yticklabels(), visible=False)
    fig.subplots_adjust(wspace=0.0)
    plt.tight_layout()
    # sys.exit()
    plt.savefig('/home/bmr135/Dropbox/K2Poles/pop_trends/seismo_comp/C6_'+labels[0]+'_'+labels[1]+'_dnu.png')
    df = df[df[params[0]]>7.5]
    # df[['EPIC',params[0],uncerts[0],params[1],uncerts[1]]].to_csv('/home/bmr135/seismo_comp/C6_'+labels[0]+'_'+labels[1]+'_dnu.txt',index=False)
    # plt.show()
    # sys.exit()

''' C3 comparisons '''
# comp_seismic_plot(YS_C3,'comp_err_nmx',['nmx','Snumax'],['nmx_err','e_Snumax'],[r'BHM',r'A2Z'],r'$\nu_{\rm{max}}$')
# comp_seismic_plot(YB_C3,'comp_err_nmx',['nmx','Bnumax'],['nmx_err','e_Bnumax'],[r'BHM',r'COR'],r'$\nu_{\rm{max}}$')
# comp_seismic_plot(BS_C3,'comp_err_nmx',['Bnumax','Snumax'],['e_Bnumax','e_Snumax'],[r'COR',r'A2Z'],r'$\nu_{\rm{max}}$')
# comp_seismic_plot(YS_C3,'comp_err_Dnu',['dnu','SDnu'],['dnu_err','e_SDnu'],[r'BHM',r'A2Z'],r'$\Delta\nu$')
# comp_seismic_plot(YB_C3,'comp_err_Dnu',['dnu','BDnu'],['dnu_err','e_BDnu'],[r'BHM',r'COR'],r'$\Delta\nu$')
# comp_seismic_plot(BS_C3,'comp_err_Dnu',['BDnu','SDnu'],['e_BDnu','e_SDnu'],[r'COR',r'A2Z'],r'$\Delta\nu$')

''' C6 comparisons '''
# comp_seismic_plot(YS_C6,'comp_err_nmx',['nmx','Snumax'],['nmx_err','e_Snumax'],[r'BHM',r'A2Z'],r'$\nu_{\rm{max}}$')
# comp_seismic_plot(YB_C6,'comp_err_nmx',['nmx','Bnumax'],['nmx_err','e_Bnumax'],[r'BHM',r'COR'],r'$\nu_{\rm{max}}$')
# comp_seismic_plot(BS_C6,'comp_err_nmx',['Bnumax','Snumax'],['e_Bnumax','e_Snumax'],[r'COR',r'A2Z'],r'$\nu_{\rm{max}}$')
comp_seismic_plot(YS_C6,'comp_err_Dnu',['dnu','SDnu'],['dnu_err','e_SDnu'],[r'BHM',r'A2Z'],r'$\Delta\nu$')
comp_seismic_plot(YB_C6,'comp_err_Dnu',['dnu','BDnu'],['dnu_err','e_BDnu'],[r'BHM',r'COR'],r'$\Delta\nu$')
comp_seismic_plot(BS_C6,'comp_err_Dnu',['BDnu','SDnu'],['e_BDnu','e_SDnu'],[r'COR',r'A2Z'],r'$\Delta\nu$')

sys.exit()
