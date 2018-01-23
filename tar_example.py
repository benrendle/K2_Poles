import numpy as np
import pandas as pd
import tar_test as tt

filename = '/home/bmr135/GA/K2Poles/Benoit.tar.gz'
sf = ''
data = 'Benoit_nodet_K2P2_C3'
delim = ','

tar = tt.TAR(filename, sf, data, delim)
z = pd.DataFrame()
z = tar()
print(z)
