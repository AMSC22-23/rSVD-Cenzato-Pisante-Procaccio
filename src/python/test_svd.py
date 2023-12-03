# To execute it:
# python test_svd.py

import scipy.linalg as la
import numpy as np

A = np.array([[4, 3, 2, 1],
            [3, 4, 3, 2],
            [2, 3, 4, 3],
            [1, 2, 3, 4]])

#B = A.T @ A
B=A
print(B)

#U = np.identity(4)
for i in range(20):
    [Q,R] = la.qr(B)
    B = R @ Q
    #U = U @ Q
    print(R)
    print(Q)
    print(B)
    print("===============")
