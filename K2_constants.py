''' Constants for K2_seismo_comp.py '''

solar_Teff = 5777 # Kelvin
solar_Numax = 3104 # muHz, Mosser 2013
solar_Dnu = 138.8 # muHz, Mosser 2013
solar_Numax_y = 3135 # pm 9, muHz
solar_Dnu_y = 135.045 # pm 0.013, muHz
solar_Numax_s = 3097.33 # pm 115.97, muHz
solar_Dnu_s = 135.2 # pm 1.55, muHz

solar_g = 27400 # cm/s^2
teffred_solar = 8907.0 # Solar red edge value, Kelvin
cadence = 29.4*60 # 58.5 sec or 29.4 mins depending on choice of SC or LC data
vnyq = 1.0/(2*cadence) * 1e6 # Nyquist frequency muHz
sys_limit = 0
dilution = 1.0
G = 6.67408E-11 # Gravitational Constant m^3 kg^-1 s^-2
