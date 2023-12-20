clc;
clear;
close all;

%Initialize matrix A
A=[
4,3,2,1;
3,4,3,2;
2,3,4,3;
1,2,3,4];

%Set the tolerance
tol=1.e-10;

%Compute U, S, V
[U,S,V] = svd_qr_test(A,tol)

%Export the matrixes on .txt files
writematrix(U, 'test_matrix_U.txt')
writematrix(S, 'test_matrix_S.txt')
writematrix(V, 'test_matrix_V.txt')
