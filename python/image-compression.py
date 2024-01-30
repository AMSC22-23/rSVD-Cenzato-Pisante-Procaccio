import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def readTXT(filename):
    in_file=open(filename,"r")
    args=list(map(lambda x: int(x),in_file.readline().split()))
    T=np.zeros(args)
    count=0
    for x in in_file:
        numbers=list(map(lambda x:float(x),x.split()))
        for i in range(len(numbers)):
            T[count,i]=numbers[i]
        count=count+1
    in_file.close()
    return T

filename_R = "/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/redS.txt"
s_R = readTXT(filename_R)
filename_G = "/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/redG.txt"
s_G = readTXT(filename_G)
filename_B = "/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/redB.txt"
s_B = readTXT(filename_B)

fig, axs =plt.subplots(2,3, figsize=(18,12))

axs[0][0].plot(s_R, 'o-')
axs[0][0].set_title('singular values of R')
plt.grid(True)

axs[1][0].plot(np.cumsum(s_R**2)/np.sum(s_R**2), '*-')
axs[1][0].set_title('explained variance of R') 
plt.grid(True)

axs[0][1].plot(s_G, 'o-')
axs[0][1].set_title('singular values of G')
plt.grid(True)

axs[1][1].plot(np.cumsum(s_G**2)/np.sum(s_G**2), '*-')
axs[1][1].set_title('explained variance of G') 
plt.grid(True)

axs[0][2].plot(s_B, 'o-')
axs[0][2].set_title('singular values of B')
plt.grid(True)

axs[1][2].plot(np.cumsum(s_B**2)/np.sum(s_B**2), '*-')
axs[1][2].set_title('explained variance of B') 
plt.grid(True)

output_path = '/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/compression_variance.txt'
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)