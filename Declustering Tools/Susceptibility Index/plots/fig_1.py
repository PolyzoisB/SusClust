import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde

########### Functions #############################################
# Define the Gaussian kernel function
def gauss(x, y, sigma_x, sigma_y):
    return np.exp(-((x - sigma_x) ** 2 + (y - sigma_y) ** 2) / 2.0) / (2.0 * np.pi * sigma_x * sigma_y)

# Generate Gaussian kernel and normalize
w = gauss(30, 30, 15, 15)
w = w / np.sum(w)

# Define xymap function
def xymap(x, y, v, type, Nx, Ny, XLim, YLim):
    if v is None:
        v = np.ones_like(x)
    
    # set boundaries
    minX, maxX = XLim
    minY, maxY = YLim
    
    # keep those within boundaries
    valid_indices = (x >= minX) & (x <= maxX) & (y >= minY) & (y <= maxY)
    x = x[valid_indices]
    y = y[valid_indices]
    v = v[valid_indices]
    
    # set the grid based on Nx, Ny
    stepx = (maxX - minX) / Nx
    stepy = (maxY - minY) / Ny
    xn = np.ceil((x - minX) / stepx).astype(int)
    yn = np.ceil((y - minY) / stepy).astype(int)
    
    xa = np.arange(minX, maxX, stepx)
    ya = np.arange(minY, maxY, stepy)
    
    M = np.zeros((Nx, Ny))
    
    for i in range(Nx):
        for j in range(Ny):
            indices = (xn == i + 1) & (yn == j + 1)
            if type == 's':
                M[i, j] = np.nansum(v[indices])
            elif type == 'm':
                M[i, j] = np.nanmean(v[indices])
            else:
                raise ValueError("type must be either 's' or 'm'")
    
    return M, xa, ya, xn, yn

########### LOAD Susceptibility Index ##################################
file_path = '../results/susc_index.txt'

df = pd.read_csv(file_path,sep=r'\s+',engine='python',header=None)
df.columns=['N-clusters', 'SI', 'Freq-NN','Threshold','Threshold-Index','N-links']
df = pd.DataFrame(df).astype(float)
Kn = df.iloc[:,0]                  # Number of clusters
ssm_index = df.iloc[:,1]           # SI
nd = df.iloc[:,2]                  # NN counter
eta = df.iloc[:,3]                 # Similarity threshold
log_eta = df.iloc[:,4].astype(int) # Similarity threshold index

################ LOAD Final Results #######################################
file_path = '../results/final_results.txt'

# Read the data
with open(file_path,'r') as file:
        lines = file.readlines()
        data = [line.split() for line in lines]
data = pd.DataFrame(data).astype(float)
data.columns=['k_best', 'SI_best', 'Thr(K_best)','K_min'\
                 ,'SI_min','Thr(K_min)','kmax','Thr(kmax)', 'kmax2'\
                    , 'Thr(kmax2)']
K_best = int(data['k_best'].values[0]) 
ind = df[df.iloc[:,0] == K_best].index[0]
log_eta_best = log_eta[ind]
eta_best = eta[ind]

K_min = int(data['K_min'].values[0]) 
ind = df[df.iloc[:,0] == K_min].index[0]
log_eta_min = log_eta[ind]
eta_min = eta[ind]

thr_best = data['Thr(K_best)'].values[0] 
thr_min = data['Thr(K_min)'].values[0] 

############ LOAD Rescaled distances #########################
file_path = '../results/rescaled_distances.txt'

# Read the data
with open(file_path,'r') as file:
        lines = file.readlines()
        data = [line.split() for line in lines]
data = pd.DataFrame(data).astype(float)
data.columns=['T', 'R']

#######################################################################################
TN = data.iloc[:,0]
RN = data.iloc[:,1]

# Convert to logarithmic scale
log_TN = np.log10(TN)
log_RN = np.log10(RN)

# Compute the 2D histogram
K, Ta, Xa, xn, yn = xymap(log_TN, log_RN, None, 's', 200, 200, [-6, 7], [-9, 1])

# Smooth the histogram using a Gaussian filter
K1 = gaussian_filter(K.T, sigma=1.5)
K1 = K1 / np.sum(K1)

# Cut lower 5% of the cumulative distribution
Ksort = np.sort(K1.flatten())
Ksum = np.cumsum(Ksort)
ii = np.min(np.where(Ksum >= 0.05))
K1[K1 < Ksort[ii]] = np.nan
K1 = K1 / (np.ptp(Ta) / K1.shape[0]) / (np.ptp(Xa) / K1.shape[1])  # Convert to density

a = np.log10(1/thr_min)
line_1 =  a * np.ones_like(Ta)-Ta
a = np.log10(1/thr_best)
line_2 =a* np.ones_like(Ta)-Ta

# Filter out the non-NaN values in K1 to determine the bounds
non_nan_indices = ~np.isnan(K1)  # Mask of non-NaN values in K1
non_nan_Ta = Ta[non_nan_indices.any(axis=0)]  # Select Ta values where any K1 column has non-NaN values
non_nan_Xa = Xa[non_nan_indices.any(axis=1)]  # Select Xa values where any K1 row has non-NaN values

# Define the x and y limits based on filtered Ta and Xa
x_min, x_max = 10**np.min(non_nan_Ta), 10**np.max(non_nan_Ta)
y_min, y_max = 10**np.min(non_nan_Xa), 10**np.max(non_nan_Xa)
######################### FIGURE ##########################################################

plt.figure(figsize=(21, 8),dpi=400)
plt.rcParams.update({'font.size': 16})  # Change 14 to the desired font size

plt.subplot(1,2,1)
h = plt.pcolormesh(10**Ta, 10**Xa, K1, shading='auto')
# plt.plot(10**Ta,10**line_1,linewidth=2,linestyle='--',color='g',label=r'$K_{min}$')
plt.plot(10**Ta,10**line_2,linewidth=3,linestyle='--',color='m',label=r'$K_{best}$')
#a = np.log10(1/(thr_best-1.0))
#line_4 =a* np.ones_like(Ta)-Ta
#line_1 = line_1-0.3*np.ones(len(line_1))
#plt.fill_between(10**Ta, 10**line_1, 10**line_4, color='grey', alpha=0.3)

plt.colorbar(h, label='Density')
plt.title('Rescaled distances', fontsize=14)
plt.xlabel('Rescaled time, T', fontsize=14)
plt.ylabel('Rescaled distance, R', fontsize=14)
#plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlim([x_min, x_max])
plt.ylim([y_min, y_max])
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.subplot(1,2,2)
log_dist = ((np.log(1/eta))/0.05).astype(int)
plt.plot(Kn, ssm_index/np.sum(ssm_index), linewidth=1.5,marker='o',color='b',label='SI')
K_full = Kn[nd>0]
nd_full = nd[nd>0]
plt.bar(K_full, nd_full/np.sum(nd_full), width=1, align='edge', facecolor='green',
        edgecolor='green',linewidth=1,label='H(K)', alpha=0.5)
plt.axvline(x=K_best, linewidth=3,color='m', linestyle='--',label=r'$K_{best}=$'f'{K_best}')
#plt.axvline(x=K_min,linewidth=2, color='g', linestyle='--',label=r'$K_{min}=$'f'{K_min}')
plt.legend(fontsize=16)
plt.title('NN Histogram',fontsize=14)
plt.xlabel('K',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.savefig('Fig_1_SI.png', bbox_inches='tight',dpi=400)
