%% load matrix
clear all;
load('crystm03.mat');
% load('finan512.mat');
%load('bcsstk13.mat');
%load('aft01.mat');
%load('1138_bus.mat');
%load('bcsstk18.mat');
%load('Chem97ZtZ.mat');
%load('mhd4800b.mat');
%load('pdb1HYS.mat');
% load('ted_B_unscaled.mat');
% load('torsion1.mat');
% load('Trefethen_20000.mat');
% load('wathen100.mat');
%load('G2_circuit.mat');
%load('thermomech_TC.mat');
%load('ecology2.mat');
%load('parabolic_fem.mat');
%load('G3_circuit.mat');
%load('Bump_2911.mat');
%load('thermal2.mat');
%load('bone010.mat');
% load('bundle1.mat');
%load('qa8fm.mat');
%load('mhd1280b.mat');
%load('hood.mat');
% load('gridgena.mat');
% load('boneS01.mat');
% load('tmt_sym.mat');
% load('Trefethen_2000.mat');
% load('wathen120.mat');
% load('cvxbqp1.mat');
%load('Pres_Poisson.mat');
% load('shallow_water2.mat');
%load('2cubes_sphere.mat');
% load('Kuu.mat');
% load('bodyy6.mat');
% load('s2rmq4m1.mat');
%load('boneS10.mat');
% load('Fault_639.mat');
% load('Flan_1565.mat');
 load('StocF-1465.mat');
%load('msdoor.mat');
%% problem setup
A=Problem.A;
p=symamd(A);
A=A(p,p);
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