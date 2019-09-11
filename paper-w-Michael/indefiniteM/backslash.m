clear all;
% ill-conditioned
%load('bcsstk13.mat');
%load('bcsstk18.mat')
%load('pdb1HYS.mat');
%load('hood.mat');
% load('cvxbqp1.mat');
% load('Fault_639.mat');
% load('cfd2.mat');
%load('sts4098.mat');
%load('StocF-1465.mat');
%load('Queen_4147.mat');
%load('Emilia_923.mat');
%% problem setup
A=Problem.A;
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