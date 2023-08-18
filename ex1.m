% Ex 1 from Homework 1
% We will use the Simpson's method
% y_(n+2)=y_n +(h/3)*(f_n+4f_(n+1)+f_(n+2))
% initial data y'=-5y    y(0)=1 
close all
clear all

odefun=@(x) -5*x;
exactsol=@(x) exp(-5*x);
y0=1;
h=0.02;
tN=10; %so we want to do 500 steps8

%compute y1 with FE
y1FE= y0 + h*odefun(y0);

%compute y1 with RK4 tableau 0  |
%                           1/2 | 1/2
%                           1/2 |  0  1/2
%                            1  |  0   0   1
%                                 1/6 1/3 1/3 1/6            
stadi=NaN(1,4);
stadi(1)=y0;
stadi(2)=y0+h/2*odefun(stadi(1));
stadi(3)=y0+h/2*odefun(stadi(2));
stadi(4)=y0+h*odefun(stadi(3));
y1RK=y0+h/6*(odefun(stadi(1))+2*odefun(stadi(2))+2*odefun(stadi(3))+odefun(stadi(4)));

%Simpson Method
[youtFE, tspan, errFE]=simpson1(odefun,y0,y1FE,h, tN);
plot(tspan,errFE,'r*')
pause
[youtRK, tspan, errRK]=simpson1(odefun,y0,y1RK,h, tN);
plot(tspan,errRK,'b*')
pause


plot(tspan,youtFE,'rx',tspan,exactsol(tspan), 'm')
legend({'FE approximation','exact solution'},'Location','northeast')
pause
plot(tspan,youtRK,'bx',tspan,exactsol(tspan), 'm')
legend({'RK approximation','exact solution'},'Location','northeast')




