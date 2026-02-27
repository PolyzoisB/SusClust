import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.io import loadmat
from datetime import datetime
from sklearn.linear_model import LinearRegression

#################### LOAD output catalog #####################
file_path = '../results/output_catalog.txt'

df = pd.read_csv(file_path,sep=r'\s+',engine='python',header=None)
df.columns=['time', 'lat', 'lon','mag','bg_id']
df = pd.DataFrame(df).astype(float)
bg = df[df['bg_id']==0]
main = df[df['mag']>=6]

# Convert date and time columns to datetime objects
dt3 = bg['time']/(3600*24*365)
quake3 = list(range(1, len(dt3) + 1))

dt4 = df['time']/(3600*24*365)
quake4 = list(range(1, len(dt4) + 1))

prc_bg = np.round((len(bg) / len(df))*100,2)
################# FIGURE ##########################################################
plt.rcParams.update({'font.size': 16})

fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True, dpi=400)
fig.patch.set_facecolor('white')

ax3=axes[0]
ax3.set_facecolor((0.7, 0.7, 0.7))
ax3.plot(dt4, df['lat'], 'o', markersize=2, markerfacecolor='cyan')
ax3.plot(main['time']/(3600*24*365), main['lat'], 'o', markersize=8,
         markerfacecolor='red',markeredgecolor='black')
ax3.set_ylim([min(df['lat']), max(df['lat'])])
ax3.set_xlim([min(dt4), max(dt4)])
ax3.set_ylabel('Latitude', fontsize=14)
ax3.set_title('Initial Catalog')

ax3_right = ax3.twinx()
ax3_right.plot(dt4, quake4, linewidth=1.5,color='m')
ax3_right.set_ylabel('Cum. Num.', fontsize=14)

ax4=axes[1]
ax4.set_facecolor((0.7, 0.7, 0.7))
ax4.plot(dt3, bg['lat'], 'o', markersize=1, color='k')
ax4.plot(main['time']/(3600*24*365), main['lat'], 'o', markersize=8,
         markerfacecolor='red',markeredgecolor='black')

ax4.set_ylim([min(bg['lat']), max(bg['lat'])])
ax4.set_xlim([min(dt4), max(dt4)])
ax4.set_title('Susceptibility-Based Declustering 'f'{prc_bg}'r'$\%$' ' bg events')
ax4.set_xlabel('Time since origin (in years)', fontsize=14)
ax4.set_ylabel('Latitude', fontsize=14)

ax4_right = ax4.twinx()
ax4_right.plot(dt3, quake3, linewidth=1.5,color='m')
ax4_right.set_ylabel('Cum. Num.', fontsize=14)

plt.savefig('Fig_2.png',bbox_inches='tight',dpi=400)
