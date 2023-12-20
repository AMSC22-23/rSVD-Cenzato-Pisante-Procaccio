function [U,S,V] = svd_qr_test(A,tol)
%SVDSIM  simple SVD program
%
% usage: [U,S,V]= svd_qr_test(A)
%     or      S = svd_qr_test(A)
%
% Input:
% A : matrix of size m x n
% tol . tolerance for the error
%
% Output:
% U : matrix of left singular vectors
% S : matrix of singular values
% V : matrix of right singular vectors
%
%
% The idea is to use the QR decomposition on A to gradually "pull" U out from
% the left and then use QR on A transposed to "pull" V out from the right.
% This process makes A lower triangular and then upper triangular alternately.
% Eventually, A becomes both upper and lower triangular at the same time,
% (i.e. Diagonal) with the singular values on the diagonal.

if ~exist('tol','var')
   tol=eps*1024;
end

%reserve space in advance
sizeA=size(A);
iter_max=100*max(sizeA);
iter=0;

% or use Bidiag(A) to initialize U, S, and V
U=eye(sizeA(1));
S=A';
V=eye(sizeA(2));

Err=realmax;
while Err>tol & iter<iter_max ;

    [Q,S]=qr(S'); 
    U=U*Q;
    [Q,S]=qr(S'); 
    V=V*Q;

%Compute the error
    E=triu(S,1);
    E=norm(E(:));
    F=norm(diag(S));
    if F==0, F=1;end
    Err=E/F;
    iter=iter+1;
end
% [Err/tol loopcount/loopmax]

%fix the signs in S
ss=diag(S);
S=zeros(sizeA);
for n=1:length(ss)
    ssn=ss(n);
    S(n,n)=abs(ssn);
    if ssn<0
       U(:,n)=-U(:,n);
    end
end

if nargout<=1
   U=diag(S);
end

return
