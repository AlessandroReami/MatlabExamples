%ex5c Homework1
% solve   y'=-Ay;
%         y(0)=ones(n,1)

%close all
%clear all

nx=100;
G=numgrid('S',nx);
A=delsq(G)*(nx-1)^2;
y0=ones(size(A,2),1);
odefun=@(y) -A*y;
y_exact = load('accurate_solution.m');
%lambda=-eigs(A,1,'lm');
t0=0;
T=0.1;
time=NaN(1,3);
error=NaN(1,3);
steps=[100,1000,10000];% altrimenti viene un errore nell'ultimo passaggio del for
counter=0; 

for h=[10^-3,10^-4,10^-5]
    counter=counter+1;
    tspan=linspace(t0,T,steps(counter)+1);
    yout=[y0,NaN(length(y0),steps(counter))];
    
    tic;
    

    c=10;
    k=h/c;
    matrixtemp=speye(size(A))+k/2*A;
    ytemp=[y0,NaN(length(y0),2*c)];
    for n=1:2*c
        %yout(:,n+1)=(eye(size(A))+h/2*A)\(yout(:,n)-h/2*A*yout(:,n));
        [ytemp(:,n+1),flag]=pcg(matrixtemp, ytemp(:,n)-k/2*A*ytemp(:,n),k^(3),1000);
    end
    yout(:,2)=ytemp(:,c+1);
    yout(:,3)=ytemp(:,2*c+1);
    
    matrix= speye(size(A))+6/11*h*A;
    for n=1:steps(counter)-2
        [yout(:,n+3),flag]=pcg(matrix, 18/11*yout(:,n+2)-9/11*yout(:,n+1)+2/11*yout(:,n), h^3, 1000);
    end
    time(counter)=toc;
    error(counter)=norm(yout(:,end)-y_exact,inf);
end
time, error, steps,




