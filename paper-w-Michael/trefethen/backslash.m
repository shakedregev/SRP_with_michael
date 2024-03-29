%% load matrix
clear all;
load('Trefethen_200000.mat');
%spy(Problem.A);
%% problem setup
%A=Problem.A;
A=top;
tol=10^-8;
nmax=length(A);
%% normalize A
C=diag(sparse(1./sqrt(diag(A))));
A=C*tril(A,-1)*C;
A=A+A'+speye(nmax);
%% solution setup
b=sparse(A*(1:nmax)'/nmax);
%b=sparse(A*ones(nmax,1));
%b=sparse(ones(nmax,1));
%% solution
tic;
x=A\b;
toc;
err=norm(A*x-b)
merr=max(abs((1:nmax)'/nmax-x));