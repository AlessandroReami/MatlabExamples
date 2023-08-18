function [x, resvec, iter]=mypcg (A, b,tol, maxit, L)

% resvec is a vector with the norm of the (absolute) residual at each iteration.


x= zeros(size(A,2),1);
r=b-A*x;
p=L.'\(L\r);
i=0;
rho= r.'*p;

condition=tol*norm(b);   %we will use the 2-norm as in MATLAB pcg function
resvec=NaN(maxit+1,1);
resvec(1)=norm(r);

while (   i<maxit   &&     resvec(i+1)>condition    )
    z=A*p;
    alfa=rho/(z.'*p);
    x=x+alfa*p;
    r=r-alfa*z;
    temp=L\r;
    g=L.'\temp;
    rho_next=r.'*g;
    beta=rho_next/rho;
    rho=rho_next;
    p=g+beta*p;
    resvec(i+2)=norm(r);

    i=i+1;
end

iter=i;
resvec=resvec(1:i+1);  % notice that length(resvec) is iter+1

% A=delsq(numgrid('S',102));
% L=ichol(A);
% n=size(A,1);
% b=A*ones(n,1);
% tol=1e-8;
% maxit=50;



