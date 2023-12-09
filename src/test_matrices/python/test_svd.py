# To execute it:
# python test_svd.py

import scipy.linalg as la
import numpy as np

A = np.array([[2,-1,0,0,0],
        [-1,2,-1,0,0],
        [0,-1,2,-1,0],
        [0,0,-1,2,-1],
        [0,0,0,-1,2]])
 
U, s, VT = np.linalg.svd(A)
print('U shape: ', U.shape)
print('s shape: ', s.shape)
print('VT shape: ', VT.shape)

print('U : ', U)
print('s : ', s)
print('V : ', VT.transpose())
