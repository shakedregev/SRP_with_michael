clear all;
% well-conditioned
load('crystm03.mat');
% load('finan512.mat');
%load('aft01.mat');
%load('1138_bus.mat');
%load('Chem97ZtZ.mat');
%load('mhd4800b.mat');
% load('ted_B_unscaled.mat');
% load('torsion1.mat');
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
% load('gridgena.mat');
% load('boneS01.mat');
% load('tmt_sym.mat');
% load('wathen120.mat');
% load('shallow_water2.mat');
%load('2cubes_sphere.mat');
% load('Kuu.mat');
% load('bodyy6.mat');
% load('s2rmq4m1.mat');
%load('boneS10.mat');
% load('Flan_1565.mat');
% load('StocF-1465.mat');
%load('Andrews.mat');
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