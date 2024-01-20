import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

ovariancancer_obs_path = r'/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/data/ovariancancer_obs.csv'
#each row contains a patient
#each coloumn a gene of that patient
ovariancancer_grp_path = r'/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/data/ovariancancer_grp.csv'
#contains labels: 1 has cancer, 0 no

A = np.genfromtxt(ovariancancer_obs_path, delimiter=',').transpose()  #generates a matrix from the data in the file
f = open(ovariancancer_grp_path)
grp = np.array(f.read().split("\n"))
grp = grp[grp != ''] #remove empty entreces

n_features=A.shape[0]
n_patients=A.shape[1]
n_cancer=np.sum(grp=='Cancer')
n_noCanc=n_patients-n_cancer

A_mean=np.mean(A,axis=1) #mean over rows 
U,s,Vt=np.linalg.svd(A-A_mean[:,np.newaxis],full_matrices=False)

fig, axs =plt.subplots(1,3, figsize=(18,6))

axs[0].loglog(s[:-1], 'o-') #tolgo ultimo ele perché è quasi 0
axs[0].set_title('singular values')

axs[1].loglog(np.cumsum(s)/np.sum(s), 'o-')
axs[1].set_title('comulative fraction')

axs[2].loglog(np.cumsum(s**2)/np.sum(s**2), 'o-')
axs[2].set_title('explained variance') 
#with 1 principal component we can explain 85% of variability

output_path = '/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/output/pca_explained_variance.png'
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)

fig = plt.figure()
for i in range(n_patients) : 
    if grp[i] == 'Cancer':
        color='r'
    else:
        color='b'
    plt.scatter(np.inner(U[:,0], A[:,i]- A_mean),
                np.inner(U[:,1], A[:,i]- A_mean),color=color)

output_path = '/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/output/pca_scatterplot_2_components.png'
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(n_patients) : 
    if grp[i] == 'Cancer':
        color='r'
    else:
        color='b'
    plt.scatter(np.inner(U[:,0], A[:,i]- A_mean),
                np.inner(U[:,1], A[:,i]- A_mean),
                np.inner(U[:,2], A[:,i]- A_mean),color=color, marker='x')
    
    ax.view_init(25,20)

output_path = '/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/output/pca_scatterplot_3_components.png'
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)