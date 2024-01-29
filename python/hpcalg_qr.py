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
    in_file.close();
    return T;

filename = "/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/dur_givens.txt"
s_R = readTXT(filename)
filename = "/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/dur_householder.txt"
s_B = readTXT(filename)

filename="/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/zeros.txt"
dim = readTXT(filename)
# Plot delle curve
plt.plot(dim, s_R, label='dur_givens')
plt.plot(dim, s_B, label='dur_householder')

# Aggiunta di etichette e legenda
plt.xlabel('Dimensione')
plt.ylabel('Durata')
plt.legend()

plt.savefig('grafico.pdf')