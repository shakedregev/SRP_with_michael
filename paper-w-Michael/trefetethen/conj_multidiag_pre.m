%% load matrix
clear all;
load('Trefethen_200000.mat');
%% problem setup
A=top;
%A=Problem.A;
p=symamd(A);
A=A(p,p);
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
resn=b-A*x;
zn=M*resn;
p=zn;
resnzn=zn'*resn;
for niter=1:nmax
    Ap=A*p;
    resz=resnzn;
    alpha=resz/(p'*Ap);
    x=x+alpha*p;
    resn=resn-alpha*Ap;
    if norm(resn)<tol
        break;
    end
    zn=M*resn;
    resnzn=zn'*resn;
    p=zn+(resnzn/resz)*p;
end
toc;
err=norm(A*x-b);
merr=max(abs((1:nmax)'/nmax-x));