%% load matrix
clear all;
% ill-conditioned
load('bcsstk13.mat');
%load('bcsstk18.mat')
%load('pdb1HYS.mat');
%load('hood.mat');
% load('cvxbqp1.mat');
% load('Fault_639.mat');
% load('cfd2.mat');
%load('sts4098.mat');
%load('Pres_Poisson.mat');
%load('StocF-1465.mat');
%load('Queen_4147.mat');
%% problem setup
A=Problem.A;
tol=1e-6;
tol2=tol^2;
nmax=length(A);
%% normalize A
C=diag(sparse(1./sqrt(diag(A))));
A=C*tril(A,-1)*C;
A=A+A'+speye(nmax);
%% solution setup
b=sparse(A*(1:nmax)'/nmax);
%b=sparse(A*ones(nmax,1));
%b=sparse(ones(nmax,1));
x_0=spalloc(nmax,1,nmax);
%% sparse inverse
lfil=ceil(nnz(A)/nmax);
tic
M=entire_r_sparse_inverse(A,nmax,lfil);
%M=alt_r_sparse_inverse(A,nmax,lfil);
M=(M+M')/2;
toc;
%R=chol(M);
%% PCG
delta=2;
waste_iter=0;
tic;
flag=1;
count=0;
while(flag==1)
    dx=spalloc(nmax,1,nmax);
    r=b-A*x_0;
    z=M*r;
    p=z;
    rho_new=real(z'*r);
    for niter=1:nmax
        q=A*p;
        beta=real(p'*q);
        if beta/(q'*q)<tol2
            if beta<=0
                disp('Matrix A is not positive definite!');
                flag=0;
                break;
            end
            disp('Matrix A is ill conditioned');
        end
        rho=rho_new;
        alpha=rho/beta;
        dx=dx+alpha*p;
        r=r-alpha*q;
        norm_r2=r'*r;
        if sqrt(norm_r2)<tol
            flag=0;
            break;
        end
        z=M*r;
        rho_new=real(z'*r);
        if rho_new/norm_r2<tol2
            if rho_new<=0
                waste_iter=waste_iter+niter;
                M=M + (delta*abs(rho)/norm_r2)*speye(nmax);
                x_0=x_0+dx;
                count=count+1;
%                 disp('M is being restarted');
%                 disp(niter);
                break;
            end
            disp('Matrix M is ill conditioned');
        end
        p=z+(rho_new/rho)*p;
    end
end
x_0=x_0+dx;
toc;
err=norm(A*x_0-b)
disp(niter+waste_iter);
disp(waste_iter);
disp(count);
merr=max(abs((1:nmax)'/nmax-x_0));
