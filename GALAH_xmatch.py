''' Cross match expected GALAH observations with asteroseismic detections
    found within C3.
'''

import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.table import Table
import sys

''' Read in GALAH data from fits file '''
data = fits.getdata('/media/ben/SAMSUNG1/GA/K2Poles/GALAH/GALAH_DR2_Public_Release.fits', 1)
t = Table(data)
df = t.to_pandas()
df.rename(columns={'star_id':'2MASS'},inplace=True)
df = df[df['flag_cannon'] == 0]
df = df.reset_index(drop=True)
print(df['2MASS'])
# df['2MASS'] = df['star_id'].map(lambda x: x.split('-')[0])
# df['2MASSv2'] = df['star_id'].map(lambda x: x.split('-')[-1])
# df['2MASS'] = df['2MASS'].convert_objects(convert_numeric=True)
# df['2MASSv2'] = df['2MASSv2'].convert_objects(convert_numeric=True)
# GALAH.to_csv('/media/ben/SAMSUNG1/GA/K2Poles/GALAH/GALAH_DR2',index=False)
# sys.exit()

# df = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/GALAH/GALAH_DR2.csv')
# df1 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/matlab_in/15_03_2018/C6_15032018_191820.csv')
df1 = pd.read_csv('/home/ben/Dropbox/K2Poles/GAP6')
df2 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/YC3_TL')
df3 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/SC3_TL')
df4 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/BC3_TL')
df5 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/matlab_in/15_03_2018/GES_15032018_191821.csv')
df6 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/matlab_in/15_03_2018/RC3_15032018_191821.csv')


df1 = df1[df1['2MASS'].notnull()]
df1['star_ID'] = df1['2MASS']
# df1['2MASS'] = df1['star_ID'].map(lambda y: y.split('-')[0])
# df1['2MASS2'] = df1['star_ID'].map(lambda y: y.split('-')[-1])
# df1['2MASS'] = df1['2MASS'].convert_objects(convert_numeric=True)
# df1['2MASS2'] = df1['2MASS2'].convert_objects(convert_numeric=True)
# print(df1['2MASS'][0],df1['2MASS2'][0])

GALAH_C6 = pd.merge(df,df1,how='inner',on=['2MASS'])

# for i in GALAH_C6['2MASS2']:
#     for j in GALAH_C6['2MASSv2']:
#         # print(i,j)
#         if i == j:
#             print(GALAH_C6['star_id'])

# GALAH_YC3 = pd.merge(df,df2,how='inner',on=['2MASS'])
# GALAH_SC3 = pd.merge(df,df3,how='inner',on=['2MASS'])
# GALAH_BC3 = pd.merge(df,df4,how='inner',on=['2MASS'])
# GALAH_GES = pd.merge(df,df5,how='inner',on=['2MASS'])
# GALAH_R3 = pd.merge(df,df6,how='inner',on=['2MASS'])

# G1 = pd.concat([GALAH_YC3,GALAH_BC3,GALAH_SC3],ignore_index=True)
# print(len(G1))
# G1 = G1.drop_duplicates(subset=['EPIC'])
# print(len(G1))
# G2 = pd.concat([GALAH_GES,GALAH_R3],ignore_index=True)
# print(len(G2))
# G2 = G2.drop_duplicates(subset=['EPIC'])
# print(len(G2))

print(len(GALAH_C6))
# print(len(GALAH_YC3))
# print(len(GALAH_SC3))
# print(len(GALAH_BC3))
# print(len(GALAH_GES))
# print(len(GALAH_R3))
