clear all
close all
c = 0.8;
d = 0.01;
rho = 50;
tstar = 1;
m = 500;
x = linspace(0,1,m)';
h = 1/(m-1);
u0 = 5*x.*(1-x).^2;
D2 = spdiags(ones(m,1)*[1,-2,1]/h^2,-1:1,m,m);
%D2(m,m-1:m) = [2,-2]/h^2;
D1 = spdiags(ones(m,1)*[-1,0,1]/(2*h),-1:1,m,m);
%D1(m,m-1:m) = [0,0];
DT = d*D2-c*D1;
Peclet = c*h/(2*d)
DT(1,1:2) = [0,0];
DT(m,m-1:m) = [0,0];
g = @(u) rho*u.*(1-u).*(u-1/2).*[0;ones(m-2,1);0]; % Dirichlet omogenee
dg = @(u) rho*((1-u).*(u-1/2)-u.*(u-1/2)+u.*(1-u)).*[0;ones(m-2,1);0];
ts = 100;
k = tstar/ts;
u = u0;
tol = k/100;
plot(x,u0)
axis([0,1,0,1])
pause(0.5)
for n = 1:ts
    F = @(y) y-u-k*DT*y-k*g(y);
    JF = @(y) speye(m)-k*DT-k*spdiags(dg(y),0,m,m);
    y = u;
    delta = -JF(y)\F(y);
    while (norm(delta,inf) > tol)
        y = y+delta;
        delta = -JF(y)\F(y);
    end
    u = y;
    plot(x,u)
    axis([0,1,0,1])
    pause(0.5)
end