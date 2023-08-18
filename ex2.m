
clear all
close all
%    ex2 homework1
%   y'=-10y^2
%   y0=1

odefun=@(x) -10*x^2;
exactsol= @(t) 1./(10*t+1);
y0=1;
T=2;
hvalues=[2^-5,2^-6,2^-7,2^-8,2^-9,2^-10];
err=NaN(1,6);

%odefun=@(t,x) -10*x^2;
%[inutile,solesatta]=ode45(odefun,[0,2],y0);
%plot(tspan,yout,'o',inutile,solesatta,'*-');

%subplot(2,3,1)
%plot(tspan,yout,'o', tspan,exactsol(tspan));

counter=0;
for h=hvalues
    counter=counter+1;
    step=T/h;
    tspan=linspace(0,T,step+1);
    yout=[y0,NaN(1,step)];
    
    for n=1:step
        stadi=NaN(1,4);
        stadi(1)=yout(n);
        stadi(2)=yout(n)+h/2*odefun(stadi(1));
        stadi(3)=yout(n)+h/2*odefun(stadi(2));
        stadi(4)=yout(n)+h*odefun(stadi(3));
        yout(n+1)=yout(n)+h/6*(odefun(stadi(1))+2*odefun(stadi(2))+2*odefun(stadi(3))+odefun(stadi(4)));
    end
    err(counter)=abs(exactsol(tspan(end))-yout(end));
    subplot(2,3,counter);
    plot(tspan,yout,'o', tspan,exactsol(tspan));
end
pause
subplot(1,1,1);
stepvalues=T./hvalues;
loglog(stepvalues,err,'*',stepvalues,err(end)*(stepvalues/stepvalues(end)).^(-4));
legend({'error','slope m= -4'},'Location','northeast')







