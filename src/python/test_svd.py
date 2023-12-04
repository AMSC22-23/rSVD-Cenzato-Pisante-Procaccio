# To execute it:
# python test_svd.py

import scipy.linalg as la
import numpy as np

A = np.array([[0.8147, 0.0975, 0.1576],
        [0.9058, 0.2785, 0.9706],
        [0.1270, 0.5469, 0.9572],
        [0.9134, 0.9575, 0.4854],
        [0.6324, 0.9649, 0.8003]])

U, s, VT = np.linalg.svd(A)
print('U shape: ', U.shape)
print('s shape: ', s.shape)
print('VT shape: ', VT.shape)

print('U : ', U)
print('s : ', s)
print('V : ', VT.transpose())
