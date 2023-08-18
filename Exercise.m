close all
clear

tic

% load input for "mesh0"
func=@(x) -4+2*x(1)^2+2*x(2)^2;
exactsol=@(x) x(1)^2+x(2)^2-x(1)^2*x(2)^2-1;
tol=10^(-8);
maxit=250;

%preallocate some variables
resvecJ=NaN(maxit+1,5);
resvecC=NaN(maxit+1,5);
iterj=NaN(5,1);
iterc=NaN(5,1);
epsilon=zeros(5,2);
time=NaN(5,1);

% n=0 is mesh0, n=1 is mesh=1, ... n=4 is mesh=4
for n=0:4
    tic
    if n==0
        bound=load("mesh0\mesh0.bound");
        coord=load("mesh0\mesh0.coord");
        topol=load("mesh0\mesh0.topol");
         %ams= averagemeshsize(coord,topol);    % average mesh size of the mesh0
    elseif n==1
            bound=load("mesh1\mesh1.bound");
            coord=load("mesh1\mesh1.coord");
            topol=load("mesh1\mesh1.topol");
    elseif n==2
            bound=load("mesh2\mesh2.bound");
            coord=load("mesh2\mesh2.coord");
            topol=load("mesh2\mesh2.topol");
    elseif n==3
            bound=load("mesh3\mesh3.bound");
            coord=load("mesh3\mesh3.coord");
            topol=load("mesh3\mesh3.topol");
    else
        bound=load("mesh4\mesh4.bound");
        coord=load("mesh4\mesh4.coord");
        topol=load("mesh4\mesh4.topol");
    end
    
    [H,f, deltai]=stiffmatrix(coord, topol, func);
    
    %impose boundary conditions
    R=10^15;
    for i=1:size(bound,1)
        c=bound(i,1);
        H(c,c) =H(c,c)*R;
        f(c)=0;
    end
    
    % solve linear system with Jacobi and Choleski 
    M=sparse(diag(diag(H)));
    [uj,flagj, relresj,iterj(n+1),resvec]=pcg(H,f,tol,maxit,M);  % x0=zeros by default setting
    resvecJ(n+1,1:iterj(n+1)+1)=resvec';

    L=ichol(H);
    [uc,flagc, relresc,iterc(n+1),resvec]=pcg(H,f,tol,maxit,L,L');  % x0=zeros by default setting
    resvecC(n+1,1:iterc(n+1)+1)=resvec';

    for i=1:length(uj)
        epsilon(n+1,1)=epsilon(n+1,1)+(uj(i)-exactsol(coord(i,:)))^2*deltai(i);
        epsilon(n+1,2)=epsilon(n+1,2)+(uc(i)-exactsol(coord(i,:)))^2*deltai(i);
    end
    epsilon(n+1,:)=[sqrt(epsilon(n+1,1)),sqrt(epsilon(n+1,2))];
    time(n+1)=toc;
end

text={'mesh0','mesh1','mesh2','mesh3','mesh4'};
epsilonj=epsilon(:,1);
epsilonc=epsilon(:,2);
table(text', time, epsilonj, iterj ,epsilonc, iterc)

timetot = toc;

for n=1:5

    semilogy(0:iterj(n),resvecJ(n,1:iterj(n)+1),'m*')
    title('PCG convergence profile using Jacobi preconditioning')
    xlabel('iterations')
    ylabel('residul norm')
    pause
    semilogy(0:iterc(n),resvecC(n,1:iterc(n)+1),'b*')
    title('PCG convergence profile using Cholesky preconditioning')
    xlabel('iterations')
    ylabel('residul norm')
    pause
end

%v=[ams,ams/2,ams/4,ams/8,ams/16]; %look also at line 27
v=[1,1/2,1/4,1/8,1/16];

loglog( v ,epsilon(:,1),'m*', v, epsilon(1,1)*(v/v(1)).^(2) )
legend({'error norm','slope 2'},'Location','northeast')
title('Convergence profile using Jacobi preconditioner')
xlabel('Average mesh size')
ylabel('Euclidian norm of the error')
pause
loglog( v ,epsilon(:,2),'b*', v, epsilon(1,2)*(v/v(1)).^(2) )
legend({'error norm','slope 2'},'Location','northeast')
title('Convergence profile using Cholesky preconditioner')
xlabel('Average mesh size')
ylabel('Euclidian norm of the error')







