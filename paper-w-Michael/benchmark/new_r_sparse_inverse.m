function [m]=new_r_sparse_inverse(A,j,tol,n,lfil)
m=sparse(n,1,lfil);
for k=1:2*lfil
    r=sparse(j,1,1,n,1)-A*m;
    if norm(r)<tol
        break;
    end
    [~, i]=max(abs(r));
    %alpha=r(i)/A(i,i);
    %m=m+r(i)/A(i,i)*sparse(i,1,1,n,1);
    m(i)=m(i)+r(i)/A(i,i);
    if nnz(m)>=lfil
        break;
    end
end
end