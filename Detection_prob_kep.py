''' Highly selective deconstruction of Mat's detection probability code'''

import numpy as np
from scipy import stats
from scipy.interpolate import interp1d
import pandas as pd
from pandas import DataFrame, read_csv
from astropy.io import fits
from astropy.table import Table
import sys
import os
# sys.path.insert(0, '/home/mxs191/Dropbox/01.12')
import granulation
import timeit
import warnings

# for TESS observation time
# from astropysics.coords.coordsys import EquatorialCoordinatesEquinox, \
# EclipticCoordinatesEquinox, GalacticCoordinates
# from itertools import groupby
# from operator import itemgetter

# glat and glon removed from globalDetections as not needed for K2

def properties(X, constants_only=False):
    # solar parameters
    # teff_solar = 5777.0 # Kelvin
    teffred_solar = 8907.0 #in Kelvin
    # numax_solar = 3104 # in micro Hz
    # dnu_solar = 138.8 # in micro Hz
    # lum = (rad**2)*((teff/teff_solar)**4)

    X['Tred'] = teffred_solar*(X['L']**-0.093) # from (3) eqn 8. red-edge temp

    return X['Tred']


def globalDetections(imag, Kp, lum, rad, teff, \
                    numax, max_T, teffred, teff_solar, teffred_solar, \
                    numax_solar, dnu_solar, sys_limit, dilution, vnyq, \
                    cadence, vary_beta=True):

    dnu = dnu_solar*(rad**-1.42)*((teff/teff_solar)**0.71) # from (14) eqn 21
    beta = 1.0-np.exp(-(teffred-teff)/1550.0) # beta correction for hot solar-like stars from (6) eqn 9.
    # if isinstance(teff, float):  # for only 1 star
    #     if (teff>=teffred):
    #         beta = 0.0
    # else:
    #     beta[teff>=teffred] = 0.0


    # to remove the beta correction, set Beta=1
    if vary_beta == False:
        beta = 1.0 # added on 04.08.16 after 02.08.16 Bill+Tiago meeting

    amp = 2.5*beta*(rad**2)*((teff/teff_solar)**0.5) # from (6) eqn 11

    env_width = 0.66 * numax**0.88 # From (5) table 2 values for delta nu_{env}. env_width is defined as +/- some value.
    # if isinstance(teff, float):  # for only 1 star
    #     if (numax>=100):
    #         env_width = numax/2
    # else:
    #     env_width[numax>100]=numax[numax>100]/2 # from (6) p12

    total, per_cam, pix_limit, npix_aper = pixel_cost(imag)

    # noise = calc_noise(imag=imag, teff=teff, exptime=cadence, e_lng=e_lng, e_lat=e_lat, \
    # g_lng=g_lng, g_lat=g_lat, sys_limit=sys_limit, npix_aper=npix_aper)
    # noise = noise*10.0**6 # total noise in units of ppm

    noise = keplerNoise(Kp)

    a_nomass = 3382*numax**-0.609 # multiply by 0.85 to convert to redder TESS bandpass.
    b1 = 0.317 * numax**0.970
    b2 = 0.948 * numax**0.992

    # call the function for the real and aliased components (above and below vnyq) of the granulation
    # the order of the stars is different for the aliases so fun the function in a loop
    Pgran, eta = granulation.granulation(numax, dilution, a_nomass, b1, b2, vnyq)
    Pgranalias = np.zeros(len(Pgran))
    etaalias = np.zeros(len(eta))

    # if vnyq is 1 fixed value
    if isinstance(vnyq, float):
        for i in range(len(numax)):

        	if numax[i] > vnyq:
        		Pgranalias[i], etaalias[i] = granulation.granulation((vnyq - (numax[i] - vnyq)), \
        		dilution, a_nomass[i], b1[i], b2[i], vnyq)

        	elif numax[i] < vnyq:
        		Pgranalias[i], etaalias[i] = granulation.granulation((vnyq + (vnyq - numax[i])), \
        		dilution, a_nomass[i], b1[i], b2[i], vnyq)

    # if vnyq varies for each star
    else:
        for i in range(len(numax)):

        	if numax[i] > vnyq[i]:
        		Pgranalias[i], etaalias[i] = granulation.granulation((vnyq[i] - (numax[i] - vnyq[i])), \
        		dilution, a_nomass[i], b1[i], b2[i], vnyq[i])

        	elif numax[i] < vnyq[i]:
        		Pgranalias[i], etaalias[i] = granulation.granulation((vnyq[i] + (vnyq[i] - numax[i])), \
        		dilution, a_nomass[i], b1[i], b2[i], vnyq[i])

    Pgrantotal = Pgran + Pgranalias

    ptot = (0.5*3.15*amp**2*((2*env_width)/dnu)*eta**2) / (dilution**2)
    Binstr = 2.0 * (noise)**2 * cadence*10**-6.0 # from (6) eqn 18
    bgtot = ((Binstr + Pgrantotal) * 2*env_width) # units are ppm**2

    snr = ptot/bgtot # global signal to noise ratio from (11)
    fap = 0.05 # false alarm probability
    pdet = 1.0 - fap
    pfinal = np.full(rad.shape[0], -99)
    # idx = np.where(max_T != 0) # calculate the indexes where T is not 0
    tlen=max_T*80*86400.0 # the length of the K2 observations in seconds

    bw=1.0 * (10.0**6.0)/tlen
    nbins=(2*env_width/bw).astype(int) # from (11)
    snrthresh = stats.chi2.ppf(pdet, 2.0*nbins) / (2.0*nbins) - 1.0
    pfinal = stats.chi2.sf((snrthresh+1.0) / (snr+1.0)*2.0*nbins, 2*nbins)

    return pfinal, snr # snr is needed in TESS_telecon2.py

def calc_noise(imag, exptime, teff, e_lng = 0, e_lat = 30, g_lng = 96, g_lat = -30, subexptime = 2.0, npix_aper = 10, \
frac_aper = 0.76, e_pix_ro = 10, geom_area = 60.0, pix_scale = 21.1, sys_limit = 0):

    omega_pix = pix_scale**2.0
    n_exposures = exptime/subexptime

    # electrons from the star
    megaph_s_cm2_0mag = 1.6301336 + 0.14733937*(teff-5000.0)/5000.0
    e_star = 10.0**(-0.4*imag) * 10.0**6 * megaph_s_cm2_0mag * geom_area * exptime * frac_aper
    e_star_sub = e_star*subexptime/exptime

    # e/pix from zodi
    dlat = (abs(e_lat)-90.0)/90.0
    vmag_zodi = 23.345 - (1.148*dlat**2.0)
    e_pix_zodi = 10.0**(-0.4*(vmag_zodi-22.8)) * (2.39*10.0**-3) * geom_area * omega_pix * exptime

    # e/pix from background stars
    dlat = abs(g_lat)/40.0*10.0**0

    dlon = g_lng
    q = np.where(dlon>180.0)
    if len(q[0])>0:
    	dlon[q] = 360.0-dlon[q]

    dlon = abs(dlon)/180.0*10.0**0
    p = [18.97338*10.0**0, 8.833*10.0**0, 4.007*10.0**0, 0.805*10.0**0]
    imag_bgstars = p[0] + p[1]*dlat + p[2]*dlon**(p[3])
    e_pix_bgstars = 10.0**(-0.4*imag_bgstars) * 1.7*10.0**6 * geom_area * omega_pix * exptime

    # compute noise sources
    noise_star = np.sqrt(e_star) / e_star
    noise_sky  = np.sqrt(npix_aper*(e_pix_zodi + e_pix_bgstars)) / e_star
    noise_ro   = np.sqrt(npix_aper*n_exposures)*e_pix_ro / e_star
    noise_sys  = 0.0*noise_star + sys_limit/(1*10.0**6)/np.sqrt(exptime/3600.0)

    noise1 = np.sqrt(noise_star**2.0 + noise_sky**2.0 + noise_ro**2.0)
    noise2 = np.sqrt(noise_star**2.0 + noise_sky**2.0 + noise_ro**2.0 + noise_sys**2.0)

    return noise2

# calculate the noise for a source from Kepler
def keplerNoise(Kp):

     c = 1.28 * 10**(0.4*(12.-Kp) + 7.)  # detections per cadence, from (4)
     noise = 1e6/c * np.sqrt(c + 9.5 * 1e5*(14./Kp)**5) # from (4) eqn 17. in ppm
     return noise

# the total number of pixels used by the highest ranked x number of targets in the tCTL
def pixel_cost(x):
    N = np.ceil(10.0**-5.0 * 10.0**(0.4*(20.0-x)))
    N_tot = 10*(N+10)
    total = np.cumsum(N_tot)

    # want to find: the number of ranked tCTL stars (from highest to lowest rank) that correspond to a pixel cost of 1.4Mpix at a given time
    per_cam = 26*4 # to get from the total pixel cost to the cost per camera at a given time, divide by this
    pix_limit = 1.4e6 # the pixel limit per camera at a given time
    #print (total/per_cam)[np.where((total/per_cam)<pix_limit)] # the total pixel cost values per camera that are less than the limit
    #print np.where((total/per_cam)<pix_limit)[0][-1] # the number of tCTL stars that can be observed before the pixel limit per cam is exceeded
    #print len(total[total/per_cam<pix_limit]) # the number of tCTL stars that can be observed before the pixel limit per cam is exceeded


    #print per_cam, pix_limit, N_tot

    # for the tCTL stars only. Not to be used when globalDetections() is called
    """
    sname1 = 'Pixel_Cost_for_top500k_ranked_targets.pdf'
    sname2 = 'Pixel_Cost_for_ALL_ranked_targets.pdf'

    fig1, ax, width, size = generalPlot()
    counter = np.arange(0, len(total), 1)
    plt.plot(counter, total/per_cam, linewidth=2)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    plt.xlabel('Number of Ranked tCTL Stars')
    plt.ylabel('Total Pixel Cost per Camera per Observation')

    # if plotting ALL ranked tCTL targets, add these lines to show the limits per cam
    #plt.axhline(y=pix_limit, c='k', linestyle='--', linewidth=2)
    #plt.axvline(x=(len(total[total/per_cam<pix_limit])), c='k', linestyle='--', linewidth=2)
    #plt.show()
    #fig1.savefig(loc+sname2)
    #plt.close()

    plt.show()
    fig1.savefig(loc+sname1)
    plt.close()
    """
    # print total.iloc[-1]
    # print per_cam
    # print pix_limit
    # print N_tot

    return total.iloc[-1], per_cam, pix_limit, N_tot
