%% load matrix
clear all;
load('FEM3D_2.mat');
%load('FEM3D_3.mat');
%load('FEM3D_4.mat');
%load('FEM3D_5.mat');
%% problem setup
A=K;
tol=10^-6;
tol2=10*eps;
nmax=length(A);
%% other normalize
C=diag(sparse(1./sqrt(diag(A))));
A=C*tril(A,-1)*C;
A=A+A'+speye(nmax);
%% solution setup
b=sparse(A*(1:nmax)'/nmax);
%b=sparse(A*ones(nmax,1));
%b=sparse(ones(nmax,1));
x=spalloc(nmax,1,nmax);
%% ichol
tic;
opts.michol = 'off';
L=ichol(A,opts);
toc;
%% PCG
tic;
Lt=L';
r=b-A*x;
y=L\r;
z=Lt\y;
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
    y=L\r;
    z=Lt\y;
    rho_new=z'*r;
    if rho_new<=0
        disp('Matrix M is not positive definite!');
        break;
    end
    p=z+(rho_new/rho)*p;
end
toc;
err=norm(A*x-b)
disp(niter);
merr=max(abs((1:nmax)'/nmax-x));