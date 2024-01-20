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

print(X.shape)
print(A.shape)