# Python testing

This folder contains two scripts in python that perform some testing on matrices using NumPy.

In order to execute the files you need of course [Python](https://www.makeuseof.com/install-python-ubuntu/) and [NumPy](https://numpy.org/install/) installed.

### General format

The format of the files should be a *.txt* file where the first line is the dimension of the matrix (**rows x cols**) and the columns are spaced using only a single ' '.

### Calculating SVD and QR

To calculate the factorizations we can execute the [file](python_test_complete.py):

```
python python_test_complete.py test_matrix.txt
```

This will produce some other files in the same directory as the file provided in input.

In particular the files generated are:

+  `test_matrix_U.txt` : matrix of left singular vectors
+  `test_matrix_S.txt` : array of singular values
+  `test_matrix_Vt.txt`: matrix transposed of right singular vectors
+  `test_matrix_Q.txt` : orthogonal matrix reduced
+  `test_matrix_R.txt` : upper triangular matrix

### Checking equality between matrices

To check if two matrices are almost identical (roundoff problems) we can execute the other [file](check_matrix_equality.py):

```
python python_test_complete.py test_matrix_1.txt test_matrix_2.txt
```

It will print if the matrices are compatible (rows and columns are the same number) and for each row it prints the difference value by value.

