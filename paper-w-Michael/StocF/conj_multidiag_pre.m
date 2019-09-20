%% load matrix
clear all;
load('StocF-1465.mat');
%% problem setup
A=Problem.A;
tol=10^-6;
tol2=10^-2;
nmax=length(A);
%% normalize A
C=diag(sparse(1./sqrt(diag(A))));
A=C*tril(A,-1)*C;
A=A+A'+speye(nmax);
%% solution setup
b=sparse(A*(1:nmax)'/nmax);
%b=sparse(A*ones(nmax,1));
%b=sparse(ones(nmax,1));
%% sparse inverse
delta=10;
fac=[20:5:100];
for kk=1:length(fac)
    lfil=fac(kk);
    tic
    M=entire_r_sparse_inverse(A,nmax,lfil);
    %M=alt_r_sparse_inverse(A,nmax,lfil);
    M=(M+M')/2;
    toc;
    %R=chol(M);
    %% PCG
    tic;
    x_0=spalloc(nmax,1,nmax);
    waste_iter=0;
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
            rho_norm=rho_new/norm_r2;
            if rho_norm<tol2
                cur_iter=niter+waste_iter;
                if abs(rho_norm)>tol2
                    fprintf('M is indefinite,      restarting. Iteration = %d \n',cur_iter);
                else
                    fprintf('M is nearly singular, restarting. Iteration = %d \n',cur_iter);
                end
                gamma = delta*(tol2-rho_norm);
                waste_iter=waste_iter+niter;
                %M=M + (delta*abs(rho)/norm_r2)*speye(nmax);
                M=M + gamma*speye(nmax);
                x_0=x_0+dx;
                count=count+1;
                %           disp(niter);
                break;
            end
            p=z+(rho_new/rho)*p;
        end
    end
    x_0=x_0+dx;
    toc;
    err=norm(A*x_0-b)
    disp(lfil);
    disp(niter+waste_iter);
    disp(waste_iter);
    disp(count);
    disp(nnz(M));
    merr=max(abs((1:nmax)'/nmax-x_0));
end