clc;
clear;
close all;

%Initialize matrix A
A=[	1, 2, 3, 7, 19;
	4, 5, 6, 7, 10;
    9, 10, 13, 12, 18;
	29, 35, 42, 15, 2;
	1, 3, 13, 71, 98];

%Set the tolerance
tol=1.e-10;

%Compute U, S, V
[U,S,V] = svd_qr_test(A,tol)

%Export the matrixes on .txt files
writematrix(U, 'test_matrix_U.txt')
writematrix(S, 'test_matrix_S.txt')
writematrix(V, 'test_matrix_V.txt')
