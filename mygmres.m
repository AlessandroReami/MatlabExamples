function [x, iter, resvec, flag] = mygmres (A, b, tol, maxit, x0)


r=b-A*x0;
k=0;
rho=norm(r);   
beta=rho;
v=(1/beta)*r;
norm_b=norm(b);
V=[v, NaN(size(A,1),maxit)];
resvec=[rho;NaN(maxit,1)];
flag=0;
H=zeros(maxit+1,maxit);

while rho>tol*norm_b    &&    k<maxit
    k=k+1;
    V(:,k+1)=A*V(:,k);
    for j=1:k
        H(j,k)=V(:,k+1).'*V(:,j);
        V(:,k+1)=V(:,k+1)-H(j,k)*V(:,j);
    end
    H(k+1,k)=norm(V(:,k+1));
    if H(k+1,k)==0
        flag=-1;
        break
    else
        V(:,k+1)=(1/H(k+1,k))*V(:,k+1);
        [Q,R]=qr(H(1:k+1,1:k));
        rho=abs(beta*Q(1,k+1));
        resvec(k+1)=rho;
    end
end
y=R(1:k,:)\(beta*Q(1,1:k).');
x=x0+V(:,1:k)*y;

resvec=resvec(1:k+1);
iter=k;


