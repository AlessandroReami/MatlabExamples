function [yout,tout, err] = simpson1(odefun,y0,y1,h, tN)
%
%     just for SCALAR ODE!!! 
%     based on y'=-5y for solving the implicit method
%     
step=(tN-h)/h;
tspan=linspace(0,tN,step+2);
yout=[y0,y1,NaN(1,step)];
exactsol=@(x) exp(-5*x);
err=NaN(1,step+2);
err(1)=abs(y0-exactsol(0));
err(2)=abs(y1-exactsol(h));

for n=1:step
    yout(n+2)=1/(1+5/3*h)*(yout(n)+h/3*(odefun(yout(n))+4*odefun(yout(n+1))));
    err(n+2)=abs(yout(n+2)-exactsol(tspan(n+2)));
end
tout=tspan;

end

