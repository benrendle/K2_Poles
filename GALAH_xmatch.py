''' Cross match expected GALAH observations with asteroseismic detections
    found within C3.
'''

import pandas as pd
import numpy as np

df = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/GALAH_C3_TL',delimiter=r'\s+')
df1 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/matlab_in/15_03_2018/C3_15032018_191820.csv')
df2 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/YC3_TL')
df3 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/SC3_TL')
df4 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/BC3_TL')
df5 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/matlab_in/15_03_2018/GES_15032018_191821.csv')
df6 = pd.read_csv('/media/ben/SAMSUNG1/GA/K2Poles/matlab_in/15_03_2018/RC3_15032018_191821.csv')

df = df[['EPIC','GALAH']]
df = df[df['GALAH'] == 1]
df = df.reset_index(drop=True)

GALAH_K2 = pd.merge(df,df1,how='inner',on=['EPIC'])
GALAH_YC3 = pd.merge(df,df2,how='inner',on=['EPIC'])
GALAH_SC3 = pd.merge(df,df3,how='inner',on=['EPIC'])
GALAH_BC3 = pd.merge(df,df4,how='inner',on=['EPIC'])
GALAH_GES = pd.merge(df,df5,how='inner',on=['EPIC'])
GALAH_R3 = pd.merge(df,df6,how='inner',on=['EPIC'])

G1 = pd.concat([GALAH_YC3,GALAH_BC3,GALAH_SC3],ignore_index=True)
print(len(G1))
G1 = G1.drop_duplicates(subset=['EPIC'],keep=False)
G2 = pd.concat([GALAH_GES,GALAH_R3],ignore_index=True)
print(len(G2))
G2 = G2.drop_duplicates(subset=['EPIC'],keep=False)

# print(len(GALAH_K2))
# print(len(GALAH_YC3))
# print(len(GALAH_SC3))
# print(len(GALAH_BC3))
# print(len(GALAH_GES))
# print(len(GALAH_R3))
print(len(G1))
print(len(G2))
