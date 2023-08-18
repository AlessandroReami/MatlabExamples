clear all
close all
tstar = 1;
m = 200;
x = linspace(0,pi/2,m)';
u0 = sin(x);
h = pi/2/(m-1);
D2 = spdiags(ones(m,1)*[1,-2,1]/h^2,-1:1,m,m);
D2(1,1:2) = [0,0];
D2(m,m-1:m) = [2,-2]/h^2;
b = @(t) 2*exp(t)*[0;sin(x(2:m))];
ts = 100;
k = tstar/ts;
[L,U,P] = lu(speye(m)-k/2*D2);
u = u0;
t = 0;
for n = 1:ts
    u = U\(L\(P*(u+k/2*D2*u+k/2*b(t)+k/2*b(t+k))));
    t = t+k;
end
tsrange = (5:5:50);
counter = 0;
for ts = tsrange
    counter = counter+1;
    k = tstar/ts;
    [L,U,P] = lu(speye(m)-k/2*D2);
    u = u0;
    t = 0;
    for n = 1:ts
        u = U\(L\(P*(u+k/2*D2*u+k/2*b(t)+k/2*b(t+k))));
        t = t+k;
    end
    err(counter) = norm(u-exp(tstar)*sin(x),inf);
end
loglog(tsrange,err,'*',tsrange,tsrange.^(-2))