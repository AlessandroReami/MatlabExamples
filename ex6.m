close all
clear all 
%
% Lotka-Volterra equations
% 
%

alfa=0.2;
beta=0.01;
gamma=0.07;
delta=0.004;
odefun=@(y) [y(1)*(alfa-beta*y(2));
             y(2)*(delta*y(1)-gamma)];
y0=[19;22];
T=300;
t0=0;
h=10^(-3);
step=(T-t0)/h;
tspan=linspace(t0,T,step+1);
yout=[y0, NaN(length(y0),step)];
for n=1:step
    stadi=NaN(length(y0),4);
    stadi(:,1)=yout(:,n);
    stadi(:,2)=yout(:,n)+h/2*odefun(stadi(:,1));
    stadi(:,3)=yout(:,n)+h/2*odefun(stadi(:,2));
    stadi(:,4)=yout(:,n)+h*odefun(stadi(:,3));
    yout(:,n+1)=yout(:,n)+h/6*(odefun(stadi(:,1))+2*odefun(stadi(:,2))+2*odefun(stadi(:,3))+odefun(stadi(:,4)));
end

plot(tspan,yout(1,:),'b*',tspan,yout(2,:),'ro')
legend({'x(t)','y(t)'},'Location','northeast')
% MANCANO LE ETICHETTE!!!
pause 
plot(yout(1,:),yout(2,:),'r')




%ts=300000;
%tspanFE=linspace(0,300,ts+1);
%hFE=300/ts;
%youtFE=[y0,NaN(2,ts)];
%for n=1:ts
%    youtFE(:,n+1)=yout(:,n)+hFE*odefun(youtFE(:,n));
%end
%plot(tspan,yout(1,:),'*',tspan,yout(2,:),'o',tspanFE,youtFE(1,:),tspanFE,youtFE(2,:));














