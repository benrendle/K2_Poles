''' Plotting program for K2_seismo_comp.py '''

import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.mlab as mlab
import pandas as pd
import sys
from pandas import DataFrame, read_csv
import numpy as np
from scipy.stats import norm
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde
from pyqt_fit import kde
import colormaps

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

def plot_CMD(X,Y,Z,lab,a):

    # X, Y, Z = input data arrays
    # lab = array of labels for the plots
    # a = colorbar parameter

    plt.scatter(X['JK'],X['Vmag'],color='0.75',label=lab[0])
    plt.scatter(Y['JK'],Y['Vmag'],color='r',alpha=0.25,label=lab[1])
    plt.scatter(Z['JK'],Z['Vmag'],label=lab[2],c=Z[a],cmap=colormaps.parula)
    plt.plot([0.5,0.5],[4,18],color='k', \
            linestyle='--',linewidth=2)
    plt.plot([0.5,1.5],[9.0,9.0],color='k',linestyle='--',linewidth=2)
    plt.plot([0.5,1.5],[15.0,15.0],color='k',linestyle='--',linewidth=2)
    cbar = plt.colorbar()
    cbar.set_label(r'$\nu_{\rm{max}}$ [$\mu$Hz]', rotation=270, fontsize=20, labelpad=25)
    cbar.ax.tick_params(labelsize=20)
    plt.ylabel(r'Vmag',fontsize=20, labelpad=20)
    plt.ylim(18,4)
    plt.xlim(-0.5,1.5)
    plt.xlabel(r'J - K',fontsize=20, labelpad=10)
    plt.tick_params(labelsize=15)
    plt.legend(prop={'size':10})    # Currently takes 3 data sets   # CMD diagram

def plot_HRD(X,lab,b):
    ''' HRD without simulated data '''
    if b == 1:
        if len(lab) > 0:
            plt.scatter(X['Teff'],X['slogg'],label=lab[0],c=X['JK'],cmap=colormaps.parula)
            cbar = plt.colorbar()
            cbar.set_label(r'J - K', rotation=270, fontsize=20, labelpad=25)
            cbar.ax.tick_params(labelsize=20)
            plt.legend(prop={'size':10})
        elif len(lab) == 0:
            plt.scatter(X['Teff'],X['slogg'],c=X['JK'],cmap=colormaps.parula)
            cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize=20)
    elif b ==0:
        if len(lab) > 0:
            plt.scatter(X['Teff'],X['slogg'],label=lab[0],color='b',alpha=0.1)#,c=X['JK'],cmap=colormaps.parula)
            plt.legend(prop={'size':10})
        elif len(lab) == 0:
            plt.scatter(X['Teff'],X['slogg'],color='b',alpha=0.1)#c=X['JK'],cmap=colormaps.parula)
    plt.xlim(max(X['Teff'])+500,min(X['Teff'])-500)
    plt.ylim(max(X['slogg'])+0.1,min(X['slogg'])-0.1)
    plt.ylabel(r'$log(g)$',fontsize=20, labelpad=20)
    plt.xlabel(r'T$_{\rm{eff}}$',fontsize=20, labelpad=10)
    plt.tick_params(labelsize=15)

def plot_simHRD(X,model,lab,b):

    # HRD diagram
    # X = input data set
    # TRILEGAL = relevant simulation from TRILEGAL to the data set
    # lab = label for legends
    # b = 0: no colourbar, 1: colourbar

    if b == 1:
        if len(lab) > 0:
            plt.scatter(10**(model['logTe']),model['logg'],color='0.75',label=lab[0])
            plt.scatter(X['Teff'],X['slogg'],label=lab[1],c=X['JK'],cmap=colormaps.parula)
            cbar = plt.colorbar()
            cbar.set_label(r'J - K', rotation=270, fontsize=20, labelpad=25)
            cbar.ax.tick_params(labelsize=20)
            plt.legend(prop={'size':10})
        elif len(lab) == 0:
            plt.scatter(10**(model['logTe']),model['logg'],color='0.75')
            plt.scatter(X['Teff'],X['slogg'],c=X['JK'],cmap=colormaps.parula)
            cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize=20)
    elif b ==0:
        if len(lab) > 0:
            plt.scatter(10**(model['logTe']),model['logg'],color='0.75',label=lab[0])
            plt.scatter(X['Teff'],X['slogg'],label=lab[1],color='b',alpha=0.1)#,c=X['JK'],cmap=colormaps.parula)
            plt.legend(prop={'size':10})
        elif len(lab) == 0:
            plt.scatter(10**(model['logTe']),model['logg'],color='0.75')
            plt.scatter(X['Teff'],X['slogg'],color='b',alpha=0.1)#c=X['JK'],cmap=colormaps.parula)
    plt.ylabel(r'$log(g)$',fontsize=20, labelpad=20)
    plt.ylim(6,min(model['logg'])-1)
    plt.xlim(8000,min(10**model['logTe'])-1000)
    plt.xlabel(r'T$_{\rm{eff}}$',fontsize=20, labelpad=10)
    plt.tick_params(labelsize=15)

def plt_comp_scat(X,a,b,a1,b1,lab,p,u,d):
    ''' Turn off truncation of combined params in K2_seismo_comp before
    using this function for desired results. Otherwise toggle off the limits
    imposed on the data here unless they are to be more stringent than
    initial cuts imposed. '''

    # X = input data set
    # a/b = input Dnu/numax from different pipelines
    # a1/b1 = uncertainties on dnu/numax
    # lab = legend labels
    # p = height of std dev makers
    # u = upper std dev marker
    # d = lower std dev marker

    fig = plt.figure(figsize=(8,4.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4,1])

    ax1 = fig.add_subplot(gs[0])
    # Difference in numax for comp plot
    X['comp_err_Dnu'] = np.sqrt(X[a1]**2+X[b1]**2)
    X['Dnu_Diff'] = ((X[a]-X[b])/X['comp_err_Dnu'])
    Y = pd.DataFrame([])
    Y = X[(X['Dnu_Diff'] > -5)]# & (X['Dnu_Diff'] > 10)]
    Y = Y[(Y['Dnu_Diff'] < 5)]# & (X['Dnu_Diff'] > 10)]
    X = X[X['Dnu_Diff'] > -10]
    X = X[X['Dnu_Diff'] < 10]
    x1 = x2 = x3 = pd.DataFrame()
    x1 = X[X['occ'] == 1.0]
    x2 = X[X['occ'] == 2.0]
    x3 = X[X['occ'] == 3.0]
    mu = np.mean(Y['Dnu_Diff'])
    std = np.std(Y['Dnu_Diff'])

    ax1.scatter(x1[a], x1['Dnu_Diff'], s=25, facecolors='none', edgecolors='b', label='Single EPIC')
    ax1.scatter(x2[a], x2['Dnu_Diff'], s=25, facecolors='none', edgecolors='r', label='Double EPIC')
    ax1.scatter(x3[a], x3['Dnu_Diff'], s=25, facecolors='none', edgecolors='g', label='Triple EPIC')

    ax1.set_ylabel(lab[0],fontsize=20,labelpad=10)
    ax1.set_xlabel(lab[1],fontsize=20, labelpad=10)
    ax1.plot([0,(max(X[a])+1)],[mu,mu],color='k',linewidth=2)
    ax1.fill_between([0,(max(X[a])+1)],(mu-std),(mu+std),facecolor='red', alpha=0.25)
    ax1.tick_params(labelsize=15)
    ax1.set_xlim(0,(max(X[a])+1))
    # ax1.set_ylim([-3,3])
    ax1.plot([0,(max(X[a])+1)],[0,0], 'k--')
    ax1.legend(loc=2,prop={'size':15})

    ax2 = fig.add_subplot(gs[1], sharey=ax1)
    n,bins,patches = ax2.hist(X['Dnu_Diff'], 50, normed=1, orientation='horizontal', \
             facecolor='b',alpha=0.25, label='Data')
    y = mlab.normpdf(bins, mu, std)
    ref = mlab.normpdf(bins, 0.0, 1.0)
    ax2.plot(y,bins,'r--',linewidth=2, label='Fit')
    ax2.plot(ref, bins, 'k--', linewidth=2, label='N(0,1)')
    bins1,y1 = std_dev_histo(mu,std,bins,y,p,u,d)
    ax2.fill(y1,bins1, facecolor='r',edgecolor='red', alpha=0.25, label=r'Std. Dev.: %s' %(std))
    ax2.plot([0,max(y)],[mu,mu],color='k',linewidth=2,label=r'$\mu$: %s' %(mu))
    # ax2.set_ylim([-3,3])
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax2.set_xticks([])
    ax2.set_xlabel('PDF',fontsize=20)
    ax2.legend(prop={'size':15})
    ax2.tick_params(labelsize=15)
    fig.subplots_adjust(wspace=0.0)
    plt.tight_layout()

def std_dev_histo(mu,std,bins,y,p,u,d):

    # Calculate the area to fill between for the vertical histogram
    # y = input data from fitted Gaussian
    # p = intercept with fitted curve
    # u = upper limit to the std dev area
    # d = lower limit to the std dev area

    x = mu+std
    z = mu-std
    diff=bins-mu
    i = np.argmin(abs(diff))
    j = i+u
    k = i-d
    h = j-k-1
    y1 = y[k:j]
    # print( y)
    # print( y1)
    bins1 = bins[k:j]
    y1[0]=0
    y1[1]=p
    y1[h-1]=p
    y1[h]=0
    bins1[0]=z
    bins1[1]=z
    bins1[h-1]=x
    bins1[h]=x

    return bins1, y1

def comp_histo(b1,b2,b3,dat,a1,a2,lab):

    # Plot comparative histograms
    # b1,b2,b3 = min bin value, max bin value, number of bins
    # dat = data to be compared - use length for number to be plotted
    # a1, a2 = variables to be plotted from different pipelines

    bin_list = np.linspace(b1,b2,b3)
    plt.figure()
    if len(dat) == 2:
        X = dat[0]
        Y = dat[1]
        n1,bins,patches=plt.hist(X[a1], bins=bin_list, facecolor='g',alpha=0.5,label=lab[0])
        n2,bins,patches=plt.hist(Y[a1], bins=bin_list, facecolor='b',alpha=0.5,label=lab[1])
        plt.xlabel(lab[2],fontsize=20, labelpad=10)
        plt.ylabel(lab[3], fontsize=20, labelpad=10)
        plt.tick_params(labelsize=15)
        plt.legend(prop={'size':15})
    elif len(dat) == 3:
        X = dat[0]
        Y = dat[1]
        Z = dat[2]
        n1,bins,patches=plt.hist(X[a1], bins=bin_list, facecolor='g',alpha=0.5,label=lab[0])
        n2,bins,patches=plt.hist(Y[a1], bins=bin_list, facecolor='b',alpha=0.5,label=lab[1])
        plt.hist(Z[a2], bins=bin_list, facecolor='r',alpha=0.25,label=lab[2])
        plt.xlabel(lab[3],fontsize=20, labelpad=10)
        plt.ylabel(lab[4], fontsize=20, labelpad=10)
        plt.tick_params(labelsize=15)
        plt.legend(prop={'size':15})
    # n=n2-n1
    # center = (bins[:-1] + bins[1:]) / 2
    # plt.figure()
    # plt.plot(center,abs(n))

def plot_KDE(a,param,label,colour,linestyle):
    density = kde.KDE1D(a[param])
    x = np.r_[min(a[param]):max(a[param]):1024j]
    plt.plot(x,density(x),label=label,color=colour,linewidth=3,linestyle=linestyle)
    plt.tick_params(labelsize=15)
    plt.ticklabel_format(useOffset=False)	# Should turn off the scientific label formatting (+n*10^x notations)
    # plt.xlabel(label,fontsize=20)
