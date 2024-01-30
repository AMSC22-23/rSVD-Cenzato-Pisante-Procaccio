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

filename = "/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/dur_g.txt"
s_R = readTXT(filename)
filename = "/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/dur_mpi_2.txt"
s_B = readTXT(filename)
filename = "/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/dur_mpi_3.txt"
s_D = readTXT(filename)
filename = "/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/dur_mpi_4.txt"
s_F = readTXT(filename)

filename="/home/giuseppepisante/HPC/AMSC-Proj/rSVD-Cenzato-Pisante-Procaccio/python/sizem.txt"
dim = readTXT(filename)
# Plot delle curve
plt.loglog(dim, s_R, label='Serial')
plt.loglog(dim, s_B, label='MPI 2 proc')
plt.loglog(dim, s_D, label='MPI 3 proc')
plt.loglog(dim, s_F, label='MPI 4 proc')


plt.grid(True)

# Aggiunta di etichette e legenda
plt.xlabel('Space')
plt.ylabel('Time')
plt.legend()

plt.savefig('mpi.pdf')