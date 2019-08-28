%% load matrix
clear all;
load('hood.mat');
%% problem setup
A=Problem.A;
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
%% changing part
fac=[91:100];
for kk=1:length(fac)
    x=spalloc(nmax,1,nmax);
    lfil=fac(kk);
    M=entire_r_sparse_inverse(A,nmax,lfil);
    M=(M+M')/2;
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
    disp(kk)
    disp(lfil)
    disp(niter)
end