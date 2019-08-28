%% load matrix
clear all;
load('FEM3D_2.mat');
%load('FEM3D_3.mat');
%load('FEM3D_4.mat');
%load('FEM3D_5.mat');
%% problem setup
A=K;
p=symamd(A);
A=A(p,p);
tol=10^-8;
nmax=length(A);
%% solution setup
b=sparse(A*(1:nmax)'/nmax);
%b=sparse(A*ones(nmax,1));
%b=sparse(ones(nmax,1));
x=spalloc(nmax,1,nmax);
%% sparse inverse
lfil=ceil(nnz(A)/nmax);
tic
M=entire_r_sparse_inverse(A,nmax,lfil);
%M=alt_r_sparse_inverse(A,nmax,lfil);
M=(M+M')/2;
toc;
%R=chol(M);
%% PCG
tic;
r=b-A*x;
z=M*r;
p=z;
rho_new=z'*r;
for niter=1:nmax
    q=A*p;
    beta=p'*q;
    if beta<=0
        disp('Matrix A is not positive definite!');
        break;
    end
    rho=rho_new;
    alpha=rho/beta;
    x=x+alpha*p;
    r=r-alpha*q;
    if norm(r)<tol
        break;
    end
    z=M*r;
    rho_new=z'*r;
    if rho_new<=0
        disp('Matrix M is not positive definite!');
        break;
    end
    p=z+(rho_new/rho)*p;
end
toc;
err=norm(A*x-b);
merr=max(abs((1:nmax)'/nmax-x));