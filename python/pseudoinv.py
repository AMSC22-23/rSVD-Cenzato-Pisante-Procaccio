#------------------------------------
#
#	Code to test Moore-Penrose pseudoinverse
#
# Command to execute: python pseudoinv.py
#
#------------------------------------

from numpy import zeros
from numpy import shape
import numpy.linalg as linalg
import sys

def exportmatrix(filename,matrix):
	out_file=open(filename,"w")
	rows,cols=shape(matrix)
	out_file.write(f'{rows} {cols}\n')
	for row in matrix:
		out_file.write(' '.join(f'{num:.4f}' for num in row))
		out_file.write('\n')
	out_file.close()
	return

def readmatrix(filename):
    in_file=open(filename,"r")
    # The first line is the number of rows and cols (some magic with lambdas happens)
    args=list(map(lambda x: int(x),in_file.readline().split()))
    # Create a numpy vector of zeros
    A=zeros(args)
    count=0
    # Read each line of the file (it is a row)
    for x in in_file:
        # Extract the numbers as a list
        numbers=list(map(lambda x:float(x),x.split()))
        # Append to the matrix the numbers
        for i in range(len(numbers)):
            A[count,i]=numbers[i]
        count=count+1
    # Close everything at the end for safety
    in_file.close()
    return A

A = readmatrix("../../A.txt")
Ainv_mia = readmatrix("../../A_inv.txt")
A_inv = linalg.pinv(A)
exportmatrix("A_inv_py.txt",A_inv)

print(linalg.norm(A_inv-Ainv_mia))