function [M]=alt_r_sparse_inverse(A,n,lfil)
M=spalloc(n,n,lfil*n);
%d=diag(A);
for j=1:n
    r=sparse(j,1,1,n,1);
    m=M*r;
    r=r-A*m;
    for k=1:2*lfil
        [~, i]=max(abs(r));
        del=r(i)/A(i,i);
        m(i)=m(i)+del;
        if nnz(m)>=lfil
            break;
        end
        r=r-del*A(:,i);
    end
    M(:,j)=m;
end
end