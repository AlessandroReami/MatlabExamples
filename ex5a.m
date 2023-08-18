%ex5a Homework1
% solve   y'=-Ay;
%         y(0)=ones(n,1)
close all
clear all

nx=100;
G=numgrid('S',nx);
A=delsq(G)*(nx-1)^2;
y0=ones(size(A,2),1);
odefun=@(t,y) -A*y;
y_exact = load('accurate_solution.m');
%lambda=-eigs(A,1,'lm')
t0=0;
T=0.1;
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
tic;
[tout,yout]=ode45(odefun,[t0,T],y0,options);
time=toc
yout=yout.';
error=norm(yout(:,end)-y_exact,inf)
 size(yout,2)



% Vsol = load('accurate_solution.m');
% nx = 100;
% A = -delsq(numgrid('S',nx))*(nx-1)^2;
% rhsfun = @(t,y) A*y;
% tspan = [0,.1];
% init = ones(size(A,1),1);
% options = odeset('RelTol',1e-06,'AbsTol',1e-08);
% tic;
% [Times,Y] = ode45(rhsfun,tspan,init,options);
% toc
% norm(Y(end,:)'-Vsol,inf)/norm(Vsol,inf)