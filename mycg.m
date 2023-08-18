function [x, resvec, iter]=mycg (A, b,tol, maxit)

% resvec is a vector with the norm of the (absolute) residual at each iteration.

x= zeros(size(A,2),1);
r=b-A*x;
p=r;
i=0;
rr=r.'*r;

%condition=tol*norm(b,"inf");
condition=tol*norm(b);
resvec=NaN(maxit+1,1);
%resvec(1)=norm(r,"inf");
resvec(1)=norm(r);

while (   i<maxit   &&     resvec(i+1)>condition    )
    z=A*p;
    alfa=rr/(z.'*p);
    x=x+alfa*p;
    r=r-alfa*z;
    rr_next=r.'*r;
    beta=rr_next/rr;
    rr=rr_next;
    p=r+beta*p;
    %resvec(i+2)=norm(r,"inf");
    resvec(i+2)=norm(r);

    i=i+1;
end

iter=i;
resvec=resvec(1:i+1);

% A=delsq(numgrid('S',102));
% L=ichol(A);
% n=size(A,1);
% b=A*ones(n,1);
% tol=1e-8;
% maxit=50;



