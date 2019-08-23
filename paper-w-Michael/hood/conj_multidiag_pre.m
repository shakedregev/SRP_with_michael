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
fac=[56 57 58 59 61 62 63];
for kk=1:length(fac)
    x=spalloc(nmax,1,nmax);
    lfil=fac(kk);
    M=entire_r_sparse_inverse(A,nmax,lfil);
    M=(M+M')/2;
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
    disp(kk)
    disp(lfil)
    disp(niter)
end