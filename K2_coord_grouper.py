''' Group K2 stars based on the channel they were observed in '''

import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon, Point, LinearRing
import sys
import colormaps

footprint_dictionary = json.load(open("/media/bmr135/SAMSUNG/GA/K2Poles/k2-footprint.json"))
df = pd.read_csv('/home/bmr135/Dropbox/GES-K2/Ages/C3_TL')
df1 = pd.read_csv('/home/bmr135/Dropbox/GES-K2/Ages/C3_New')
df2 = pd.read_csv('/home/bmr135/Dropbox/GES-K2/Ages/APK2')
df3 = pd.read_csv('/home/bmr135/Dropbox/GES-K2/Ages/C3_spec')


print(len(df1))
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(18,18))

x = np.linspace(0,max(df1['mass'])+0.05)
# ax1.fill_between(x, 0.1, max(df2['Z']), facecolor='gray', interpolate=True,label=r'Kepler Z range')
# ax1.plot([2.25,2.25],[0.1,1.5],color='k',linewidth=2,label=r'Kepler Z range')
# ax1.plot([2.24,2.26],[0.1,0.1],color='k',linewidth=2,label=None)
# ax1.plot([2.24,2.26],[1.5,1.5],color='k',linewidth=2,label=None)
im = ax1.scatter(df1['mass'],df1['Z'],c=df1['feh'],cmap=colormaps.parula,label=None)
ax1.scatter(df2['mass'],df2['Z'],c=df2['feh'],cmap=colormaps.parula,label=None)
cbar = fig.colorbar(im,ax=ax1)
cbar.set_label(r'[Fe/H]', rotation=270, labelpad=25, fontsize=20)
ax1.set_xlabel(r'\begin{center}Mass [M$_{\odot}$]\\Fig. 1 - Mass-Z plot for the C3(Z$<0$) and \textit{Kepler} APOKASC (Z$>$0)\\ photometric samples, with photometric $[\rm{Fe/H}]$ colourbar\end{center}',fontsize=20)
ax1.set_ylabel(r'Z [kpc]',fontsize=20)
# ax1.set_xlim(min(df1['mass'])-0.05,max(df1['mass'])+0.05)
# ax1.set_rasterized(True)
ax1.legend()


ax2.hist(df1['age'],bins=np.linspace(0,20,100),histtype='step',normed=True,label=r'C3 Photometry',linewidth=2)
ax2.hist(df3['age'],bins=np.linspace(0,20,100),histtype='step',normed=True,label=r'RAVE C3',linewidth=2)#,alpha=0.65)
ax2.set_xlabel(r'\begin{center}log$_{10}$(Age)\\Fig. 2 - Normalised log$_{10}$(Age) distribution for the C3\\ photometric (blue) and RAVE/Gaia-ESO C3\\ (orange) samples.\end{center}',fontsize=20)
# cur_axes = fig.gca()
ax2.axes.get_yaxis().set_ticklabels([])
ax2.axes.get_yaxis().set_ticks([])
ax2.legend()

df = df[df['Vmag'] >= 9]
df = df.reset_index()

ax3.hist(df['Vmag'],bins=[9,10,11,12,13,14,15])
ax3.set_xlabel(r'\begin{center}V\\Fig. 3 - Distribution of C3 V-band magnitudes.\end{center}',fontsize=20)

df = df[df['Vmag'] >= 13]
df = df.reset_index()

v13 = df[(df['Vmag'] >= 13) & (df['Vmag'] < 14)]
v14 = df[(df['Vmag'] >= 14) & (df['Vmag'] <= 15)]
ax4.scatter(v13['RA'],v13['Dec'],label=r'V = 13-14',color='c')#,alpha=0.3)
ax4.scatter(v14['RA'],v14['Dec'],label=r'V = 14-15',color='m')#,alpha=0.3)
ax4.set_xlabel(r'\begin{center}RA\\Fig. 4 - Sky map in RA and DEC of primary targets in the C3 field.\\ Orange boxes show the nearest neighbour observation groupings.\end{center}',fontsize=20)
ax4.set_ylabel(r'DEC',fontsize=20)

for j in footprint_dictionary["c3"]["channels"]:
    df2 = pd.DataFrame()
    mychannel = footprint_dictionary["c3"]["channels"][j]

    ax4.plot(mychannel["corners_ra"] + mychannel["corners_ra"][:1], \
            mychannel["corners_dec"] + mychannel["corners_dec"][:1], \
            color='orange',label=None)

    polygon = Polygon([(mychannel['corners_ra'][0], mychannel['corners_dec'][0]), \
                    (mychannel['corners_ra'][1], mychannel['corners_dec'][1]), \
                    (mychannel['corners_ra'][2], mychannel['corners_dec'][2]), \
                    (mychannel['corners_ra'][3], mychannel['corners_dec'][3])])
    linearring = LinearRing(list(polygon.exterior.coords))
    plt.text(polygon.centroid.coords[0][0],polygon.centroid.coords[0][1],str(j),fontsize=15)
    for i in range(len(df)):
        point = Point(df['RA'][i],df['Dec'][i])
        if polygon.contains(point):
            df2.append(df.iloc[i], ignore_index=True)
        if polygon.touches(point):
            df2.append(df.iloc[i], ignore_index=True)
    # df1.to_csv('/home/bmr135/Dropbox/GES-K2/Ages/C6_Channels/C6_channel_'+j,index=False)
ax4.legend()

fig.subplots_adjust(hspace=0.52,top=0.97,bottom=0.14)
fig.set_rasterized(True)
plt.show()
# fig.savefig('/home/bmr135/Dropbox/GES-K2/Ages/C6_combined.eps',reasterized=True,dpi=550)
