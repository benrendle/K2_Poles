import numpy as np

# calculate the granulation at a set of frequencies from (8) eqn 2 model F
def granulation(nu0, dilution, a_nomass, b1, b2, vnyq):

     # Divide by dilution squared as it affects stars in the time series.
     # The units of dilution change from ppm to ppm^2 microHz^-1 when
    #  going from the
     # time series to frequency. p6: c=4 and zeta = 2*sqrt(2)/pi
     Pgran = (((2*np.sqrt(2))/np.pi) * (a_nomass**2/b1) / (1 + ((nu0/b1)**4)) \
     + ((2*np.sqrt(2))/np.pi) * (a_nomass**2/b2) / (1 + ((nu0/b2)**4))) / (dilution**2)

     # From (9). the amplitude suppression factor. Normalised sinc with pi (area=1)
     eta = np.sinc((nu0/(2*vnyq)))

     # the granulation after attenuation
     Pgran = Pgran * eta**2

     return Pgran, eta


# calculate the total granulation from real and aliased components
def TotalGranulation(bins, dilution, numax, vnyq, bpass):

     # multiply by 'bpass' to convert to the redder TESS bandpass. Eqns from (8)
     a_nomass = bpass * 3382*numax**-0.609
     b1 = bpass * 0.317 * numax**0.970
     b2 = bpass * 0.948 * numax**0.992

     Pgran, eta = granulation(bins[:,0], dilution, a_nomass, b1, b2, vnyq)

     if numax >= vnyq:
         Pgranalias, etaalias = granulation((vnyq - (bins[:,0] - vnyq)), \
         dilution, a_nomass, b1, b2, vnyq)
     elif numax < vnyq:
         Pgranalias, etaalias = granulation((vnyq + (vnyq - bins[:,0])), \
         dilution, a_nomass, b1, b2, vnyq)

     # The total granulation = the sum of the real and aliased components
     Pgrantotal = Pgran + Pgranalias
     return Pgrantotal
