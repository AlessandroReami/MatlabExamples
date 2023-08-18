function [tout,yout]=rkes7(odefun,tspan,y0)
%   function [tout,yout]=rkes7(odefun,tspan,y0)
%   applica un Runge kutta di tableau:
%     0| 
%   1/2| 1/2
%   1/2|  0  1/2
%     1|  0   0   1
%        1/6 1/3 1/3 1/6

yout=NaN( length(y0), length(tspan));
yout(:,1)=y0;
m=length(tspan)-1;

a=[ 0, 0, 0, 0; 0.5, 0, 0, 0; 0, 0.5, 0, 0; 0, 0, 1, 0];
b=[1/6, 1/3, 1/3, 1/6].';
c=[0, 1/2, 1/2, 1].';
stadi=NaN(length(y0),length(c));
funvalue=zeros(length(y0),length(c));

for n=1:m  %ciclo per avere yn+1
    k=tspan(n+1)-tspan(n);
    yout(:,n+1)=yout(:,n);
    for i=1:length(c)     %ciclo per avere gli stadi
        stadi(:,i)=yout(:,n);
        for j=1:i-1       %ciclo per calcolare la sommatoria delle valutazioni di f che producono uno stadio
            stadi(:,i)=stadi(:,i)+a(i,j)*funvalue(:,j);
        end
        funvalue(:,i)=odefun(tspan(n)+k*c(i),stadi(:,i));
        yout(:,n+1)=yout(:,n+1)+k*b(i)*funvalue(:,i);
    end
end
yout=yout.';
tout=tspan;
   