#------------------------------------
#
#	Code to generate SVD and QR factorization with numpy
#
# Command to execute: python python_test_complete.py filename.txt
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
	out_file.write(' '.join(str(num) for num in vector))
	out_file.close()
	return 

def print_matrix_to_file(filename,matrix):
	out_file=open(filename,"w")
	rows,cols=shape(matrix)
	out_file.write(f'{rows} {cols}\n')
	for row in matrix:
		out_file.write(' '.join(str(num) for num in row))
		out_file.write('\n')
	out_file.close()
	return

#	Help command
if len(sys.argv)==2 and sys.argv[1]=='-h':
	print("For a correct usage of the script execute:")
	print("	> python python_test_complete.py [filename.txt]")
	exit(0)

# Open the file in input (if present)
if len(sys.argv)!=2:
	print("Invalid number of arguments >:(")
	exit(1)
filename=sys.argv[1]
in_file=open(filename,"r")
# The first line is the number of rows and cols (some magic with lambdas happens)
args=list(map(lambda x: int(x),in_file.readline().split()))
# Create a numpy vector of zeros
a=zeros(args)
count=0
# Read each line of the file (it is a row)
for x in in_file:
	# Extract the numbers as a list
	numbers=list(map(lambda x:float(x),x.split()))
	# Append to the matrix the numbers
	for i in range(len(numbers)):
		a[count,i]=numbers[i]
	count=count+1
# Close everything at the end for safety
in_file.close();

print("This is the matrix:")
print(a)
# Compute the svd using numpy
U,S,Vt=linalg.svd(a) # Can be also full_matrices=False

Q,R=linalg.qr(a,'reduced') # Can be also 'complete'

# Split the filename to add the extensions
filename=filename.split(".")

# Print the results to a file
print_matrix_to_file(filename[0]+"_U."+filename[1],U)
print_vector_to_file(filename[0]+"_S."+filename[1],S)
print_matrix_to_file(filename[0]+"_Vt."+filename[1],Vt)
print_matrix_to_file(filename[0]+"_Q."+filename[1],Q)
print_matrix_to_file(filename[0]+"_R."+filename[1],R)
print("Successfully printed everything, check your folder :)")
