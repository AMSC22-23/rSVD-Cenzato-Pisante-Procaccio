#------------------------------------------
#
# Python code to check if two matrices are approximately equal
#
# Command to execute: python check_matrix_equality.py file1.txt file2.txt
#
#------------------------------------------

# Python import
import sys

# Help command
if len(sys.argv)==2 and sys.argv[1]=='-h':
	print("For a correct usage of the script execute:")
	print("	> python check_matrix_equality.py [filename1.txt] [filename2.txt]")
	exit(0)

# If the user did not put the right arguments adios
if len(sys.argv)!=3:
	print("Invalid number of arguments >:(")
	exit(1)

filename1=sys.argv[1]
filename2=sys.argv[2]

# Open the files
file1=open(filename1,"r")
file2=open(filename2,"r")

# Extract first row
args1=list(map(lambda x: int(x),file1.readline().split()))
args2=list(map(lambda x: int(x),file2.readline().split()))

# If the first row is different then no need to check anything
if args1[0]!=args2[0] or args1[1]!=args2[1]:
	print("Incompatible matrices, no check needed.")
	exit(1)

row_err=0.
tot_err=0.
#	Read each row of the file
for row in range(args1[0]):

	row_err=0.

	numbers1=list(map(lambda x: float(x),file1.readline().split()))
	numbers2=list(map(lambda x: float(x),file2.readline().split()))

	#	Compute the error of each row and print it
	for idx in range(len(numbers1)):
		row_err+=abs(numbers1[idx]-numbers2[idx])
	
	print(f'Row {idx} error: {row_err}')
	tot_err+=row_err

print(f'Total error: {tot_err}')
