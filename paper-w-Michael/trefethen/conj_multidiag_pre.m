%% load matrix
clear all;
load('Trefethen_2000.mat');
%% problem setup
A=Problem.A;
%A=top;
tol=10^-6;
tol2=10*eps;
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
r=b-A*x;
z=M*r;
p=z;
rho_new=z'*r;
for niter=1:nmax
    q=A*p;
    beta=real(p'*q);
    if beta<tol2
        if beta<=0
            disp('Matrix A is not positive definite!');
            break;
        end
        disp('Matrix A is ill conditioned');
    end
    rho=rho_new;
    alpha=rho/beta;
    x=x+alpha*p;
    r=r-alpha*q;
    if norm(r)<tol
        break;
    end
    z=M*r;
    rho_new=real(z'*r);
    if rho_new<tol2
        if rho_new<=0
            disp('Matrix M is not positive definite!');
            break;
        end
        disp('Matrix M is ill conditioned');
    end
    p=z+(rho_new/rho)*p;
end
toc;
err=norm(A*x-b)
disp(niter)
merr=max(abs((1:nmax)'/nmax-x));