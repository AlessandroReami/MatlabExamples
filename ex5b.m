%ex5b Homework1
% solve   y'=-Ay;
%         y(0)=ones(n,1)

close all
%clear all

nx=100;
G=numgrid('S',nx);
A=delsq(G)*(nx-1)^2;
y0=ones(size(A,2),1);
odefun=@(t,y) -A*y;
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
    matrix=speye(size(A))+h/2*A;
    for n=1:steps(counter)
        %yout(:,n+1)=(eye(size(A))+h/2*A)\(yout(:,n)-h/2*A*yout(:,n));
        [yout(:,n+1),flag]=pcg(matrix, yout(:,n)-h/2*A*yout(:,n),h^(3),1000);
    end
    time(counter)=toc;
    error(counter)=norm(yout(:,end)-y_exact,inf);
end
time, error, steps,


%runna ma ci mette 300 seconti =5 min per il primo caso, per fare gli altri
%servirebbero 50minuti e 500 min=8h e 20 min