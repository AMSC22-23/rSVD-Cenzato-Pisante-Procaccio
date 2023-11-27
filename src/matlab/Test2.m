clc
clear
close all

A=[0.8147, 0.0975, 0.1576;
0.9058, 0.2785, 0.9706;
0.1270, 0.5469, 0.9572;
0.9134, 0.9575, 0.4854;
0.6324, 0.9649, 0.8003];

[Q R]=qr(A)

writematrix(Q, 'test_matrix_Q.txt')
writematrix(R, 'test_matrix_R.txt')