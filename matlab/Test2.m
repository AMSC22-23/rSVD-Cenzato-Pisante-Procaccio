clc
clear
close all
m=10;
n=10;

for  j = 1: n
    
        for i = 1:m
        
            A(i, j) = 2 * (i-1) - (j-1);
        
    
        end
end

B=[
4,3,2,1;
3,4,3,2;
2,3,4,3;
1,2,3,4];

C=[12, -51, 4;
6, 167, -68;
-4, 24, -41];

[Q R]=qr(A);

writematrix(Q, 'test_matrix_Q.txt')
writematrix(R, 'test_matrix_R.txt')
