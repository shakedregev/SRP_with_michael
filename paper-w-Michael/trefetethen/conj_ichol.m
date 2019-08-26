%% load matrix
clear all;
load('Trefethen_2000000.mat');
%% problem setup
A=top;
%A=Problem.A;
% p=symamd(A);
% A=A(p,p);
tol=10^-8;
nmax=length(A);
% %% normalize A
% %C=diag(sparse(1./sqrt(diag(A))));
% tic
% C=diag(sparse(sqrt(1./diag(A))));
% A=C*A*C;
% A=(A+A')/2;
% toc;
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
resn=b-A*x;
y=L\resn;
zn=Lt\y;
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
    y=L\resn;
    zn=Lt\y;
    resnzn=zn'*resn;
    p=zn+(resnzn)/(resz)*p;
end
toc;
err=norm(A*x-b);
merr=max(abs((1:nmax)'/nmax-x));