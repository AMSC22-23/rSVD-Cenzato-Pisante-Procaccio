function [u,s,v] = svd_qr_test(a,tol)
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
sizea=size(a);
loopmax=100*max(sizea);
loopcount=0;

% or use Bidiag(A) to initialize U, S, and V
u=eye(sizea(1));
s=a';
v=eye(sizea(2));

Err=realmax;
while Err>tol & loopcount<loopmax ;
%   log10([Err tol loopcount loopmax]); pause
    [q,s]=qr(s'); u=u*q;
    [q,s]=qr(s'); v=v*q;

% exit when we get "close"
    e=triu(s,1);
    E=norm(e(:));
    F=norm(diag(s));
    if F==0, F=1;end
    Err=E/F;
    loopcount=loopcount+1;
end
% [Err/tol loopcount/loopmax]

%fix the signs in S
ss=diag(s);
s=zeros(sizea);
for n=1:length(ss)
    ssn=ss(n);
    s(n,n)=abs(ssn);
    if ssn<0
       u(:,n)=-u(:,n);
    end
end

if nargout<=1
   u=diag(s);
end

return
