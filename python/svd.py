#------------------------------------
#
#	Code to test SVD with numpy
#
# Command to execute: python svd.py
#
#------------------------------------

# Python imports
from numpy import zeros
from numpy import shape
import numpy.linalg as linalg
import sys

#Utility functions to write to a file
def print_vector_to_file(filename,vector):
	out_file=open(filename,"w")
	out_file.write(f'1 {len(vector)}\n')
	out_file.write(' '.join(f'{num:.4f}' for num in vector))
	out_file.close()
	return 

def print_matrix_to_file(filename,matrix):
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

filename = "../src/A.txt"
A = readmatrix(filename)
# Compute the svd using numpy
U,s,Vt=linalg.svd(A,full_matrices=False) 

s_mio=readmatrix("../src/s_pm.txt")
print("Difference eigenvalues :")
print(linalg.norm(s-s_mio))

# Print the results to a file
# Split the filename to add the extensions
#print_matrix_to_file(filename[0]+"_U"+ext,U)
#print_vector_to_file("s_py.txt",s)
#print_matrix_to_file(filename[0]+"_V"+ext,Vt.transpose())
#print("Successfully printed everything, check your folder :)")