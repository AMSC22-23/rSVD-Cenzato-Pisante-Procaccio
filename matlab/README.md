# Matlab testing

This folder contains two scripts in matlab that perform some testing on matrices qr function provided by the matlab library.

In order to execute the files you need of course matlab installed [https://it.mathworks.com/help/install/ug/install-products-with-internet-connection.html].

To calculate the factorizations we can execute the [file](Test1.m), which uses the function svd_qr_test provided in the script [file](svd_qr_test.m). 

This will produce some other files in the same directory as the files provided in input.

In particular the files generated are:

+  `test_matrix_U.txt` : matrix of left singular vectors
+  `test_matrix_S.txt` : array of singular values
+  `test_matrix_V.txt` : matrix of right singular vectors
+  `test_matrix_Q.txt` : orthogonal matrix reduced
+  `test_matrix_R.txt` : upper triangular matrix
