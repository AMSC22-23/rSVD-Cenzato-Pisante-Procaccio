import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

ovariancancer_obs_path = r'/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/data/ovariancancer_obs.csv'
#each row contains a patient
#each coloumn a gene of that patient
ovariancancer_grp_path = r'/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/data/ovariancancer_grp.csv'
#contains labels: 1 has cancer, 0 no

A = np.genfromtxt(ovariancancer_obs_path, delimiter=',').transpose() 
f = open(ovariancancer_grp_path)
grp = np.array(f.read().split("\n"))
grp = grp[grp != ''] 
n_patients=A.shape[1]
n_cancer=np.sum(grp=='Cancer')

# Read T
filename = '/home/ari/rSVD-Cenzato-Pisante-Procaccio/src/T.txt'
in_file=open(filename,"r")
# The first line is the number of rows and cols
args=list(map(lambda x: int(x),in_file.readline().split()))
T=np.zeros(args)
count=0
for x in in_file:
	numbers=list(map(lambda x:float(x),x.split()))
	for i in range(len(numbers)):
		T[count,i]=numbers[i]
	count=count+1
in_file.close();

fig = plt.figure()
for i in range(n_patients) : 
    if i <= n_cancer:
        color='r'
    else:
        color='b'
    plt.scatter(T[0,i], T[1,i],color=color)

output_path = '/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/output/pca_c++_2.png'
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(n_patients) : 
    if grp[i] == 'Cancer':
        color='r'
    else:
        color='b'
    plt.scatter(-T[0,i], T[1,i], T[2,i],color=color, marker='x')
    
    ax.view_init(25,20)

output_path = '/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/output/pca_c++_3.png'
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)