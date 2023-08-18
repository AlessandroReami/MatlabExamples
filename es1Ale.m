%es1 da solo
clear all
close all

%discretizzo nello spazio
m=500;
h=pi/2/(m-1);
x=linspace(0,pi/2,m).';
A=toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,m));
b=@(t) 2*exp(t)*[0;sin(x(2:m))];
u0=sin(x);  %soddisfa le condizioni al bordo
A(1,1:2)=[0,0];
A(m,m-1)=2/h^2;


%ora sto risolvendo il sistema u'(t)=Au(t)+b(t)
%uso trapezi
tstar=1;
ts=100;   %timesteps
t=0;      %istante di tempo corrente
k= (tstar-t)/ts;  %lunghezza del passo

%ricavo una fattorizzazione della matrice che useremo nel metodo implicito
[L,U,P]=lu(speye(m)-k/2*A);   %LU=PB

u=u0;
for n=1:ts
    u=U\(L\(P*(u+k/2*A*u+k/2*b(t)+k/2*b(t+k))));
    t=t+k;
end

plot(x,u,'*',x,exp(t)*sin(x))
pause

tsrange=[10,20,30,40];
err=NaN(1,length(tsrange));
counter=0;
for ts=tsrange
    counter=counter+1;
    t=0;      %istante di tempo corrente
    k= (tstar-t)/ts;  %lunghezza del passo
    [L,U,P]=lu(speye(m)-k/2*A);
    u=u0;
    for n=1:ts
        u=U\(L\(P*(u+k/2*A*u+k/2*b(t)+k/2*b(t+k))));
        t=t+k;
    end
    err(counter)=norm(u-exp(tstar)*sin(x),inf);
end

loglog(tsrange,err,'*',tsrange,err(end)*(tsrange/tsrange(end)).^(-2))