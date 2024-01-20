from matplotlib.image import imread
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.figsize'] = [8, 8]

# write here the import path of the image
image_path = r"/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/data/TarantulaNebula.jpg"

A = imread(image_path)

img = plt.imshow(A)
plt.axis('off')

X = np.mean(A, axis=2) #this converts it in gray scale
img = plt.imshow(X)
plt.axis('off')
img.set_cmap('gray') #just for visualization

output_path = '/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/output/grayscale_image.png'
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)

U, s, Vt = np.linalg.svd(X, full_matrices=False)

fig, axs =plt.subplots(1,3, figsize=(18,6))

#plt.plot(s, 'o-')
#plt.semilogy(s, 'o-') #to visualize it better
axs[0].loglog(s, 'o-') #even better
axs[0].set_title('singular values')

axs[1].semilogx(np.cumsum(s)/np.sum(s), 'o-')
axs[1].set_title('comulative fraction')

axs[2].semilogx(np.cumsum(s**2)/np.sum(s**2), 'o-')
axs[2].set_title('explained variance')

output_path = '/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/output/explained_variance.png'
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)

#some ranks to see differences
fig, axs = plt.subplots(5,6,figsize=(18,12))
axs=axs.flatten() 
idxs=[1,5,10,50,100,500]

for i in range(len(idxs)):
    k=idxs[i]
    Xk = U[:,:k] @ np.diag(s[:k]) @ Vt[:k,:]
    plt.set_cmap('gray') 
    axs[i].imshow(Xk)
    axs[i].axis('off')
    axs[i].set_title('k = %d' %k)

output_path = '/home/ari/rSVD-Cenzato-Pisante-Procaccio/python/output/compressions_image.png'
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)