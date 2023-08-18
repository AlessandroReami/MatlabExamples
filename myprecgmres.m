function [x, iter, resvec, flag] = myprecgmres (A, b, tol, maxit, x0, L, U)

% we want to solve Ax=b
%                  (LU)^(-1)Ax=b
%                   .....
%                   L^(-1)AU^(-1) (Ux)=L^(-1)b


% r=L\(b-A*x0);
% k=0;
% rho=norm(r);   
% beta=rho;
% v=(1/beta)*r;
% norm_b=norm(U\(L\b));
% 
% V=[v, NaN(size(A,1),maxit)];
% resvec=[rho;NaN(maxit,1)];
% flag=0;
% H=zeros(maxit+1,maxit);
% 
% while rho>tol*norm_b    &&    k<maxit % MODIFICA LA RICHIESTA SULLA
%                                       % TOL!!!!!!!!!!
%                                       %
%                                       %
%                                       %
%     k=k+1;
% 
%     %temp1=U\v_vec(:,k); 
%     %temp2=A*temp1;      
%     %v_vec(:,k+1)=L\temp2;
%     V(:,k+1)=L\(A*(U\V(:,k)));
% 
%     for j=1:k
%         H(j,k)=V(:,k+1).'*V(:,j);
%         V(:,k+1)=V(:,k+1)-H(j,k)*V(:,j);
%     end
%     H(k+1,k)=norm(V(:,k+1));
%     if H(k+1,k)==0
%         flag=-1;
%         break
%     else
%         V(:,k+1)=(1/H(k+1,k))*V(:,k+1);
%         [Q,R]=qr(H(1:k+1,1:k));
%         rho=abs(beta*Q(1,k+1));
%         resvec(k+1)=rho; %commentami
% 
%         %calcola a mano norm(r)
%     end
% end
% y=R(1:k,:)\(beta*Q(1,1:k).');
% x=x0+U\(V(:,1:k)*y);
% 
% resvec=resvec(1:k+1);
% iter=k;
% 



r=L\(b-A*x0);
k=0;
rho=norm(r);   
beta=rho;
v=(1/beta)*r;
norm_b=norm(U\(L\b));

V=[v, NaN(size(A,1),maxit)];
resvec=[rho;NaN(maxit,1)];
flag=0;
H=zeros(maxit+1,maxit);

while rho>tol*norm_b    &&    k<maxit 
    k=k+1;

    %temp1=U\v_vec(:,k); 
    %temp2=A*temp1;      
    %v_vec(:,k+1)=L\temp2;
    V(:,k+1)=L\(A*(U\V(:,k)));

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
        %rho=abs(beta*Q(1,k+1));
        %resvec(k+1)=rho; 

        y=R(1:k,:)\(beta*Q(1,1:k).');
        x=x0+U\(V(:,1:k)*y);
        r=b-A*x;
        resvec(k+1)=norm(L\(r));
        rho=norm(U\(L\r));

    end
end
%y=R(1:k,:)\(beta*Q(1,1:k).');
%x=x0+U\(v_vec(:,1:k)*y);

resvec=resvec(1:k+1);
iter=k;

