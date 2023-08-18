%es5
close all
clear all

mrif=3^8+1;
m=mrif;
h=1/(m-1);
A=toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
B=toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));
A(1,2)=0;
A(m,m-1)=0;
B(1,2)=0;
B(m,m-1)=0;
b=[1,zeros(1,m-2),1].';
x=linspace(0,1,m).';

M=(ones(m,1)+x).*A+B;                      %e fa quello che voglio perche' sto punto-moltiplicando un vettore colonna per una matrice e questa operazione 
%M=A(1,:)*(1+x)+B;  ma serve un ciclo!!     %consiste nel moltiplicare ogni riga della matrice
                                           % per corrispettivo elemento del vettore
u=M\(b-ones(m,1));
plot(x,u,'*')
hold on
solesatta=@(x) -x+1/log(2)*log(1+x);
plot(x,solesatta(x))
pause

%guardiamo ora la convergenza

mrange=2:1:7;
mrange=(3.^mrange)+1;
counter=0;
err=NaN(1,length(mrange)+1);

%ricavo l'ultimo valore che ho gia calcolato
err(1,7)=norm(u-solesatta(x),inf);
for n=mrange
    m=n;
    counter=counter+1;
    h=1/(m-1);
    A=toeplitz(sparse([1,1],[1,2],[-2,1]/(h^2),1,m));
    B=toeplitz(sparse(1,2,-1/(2*h),1,m),sparse(1,2,1/(2*h),1,m));
    A(1,2)=0;
    A(m,m-1)=0;
    B(1,2)=0;
    B(m,m-1)=0;
    b=[1,zeros(1,m-2),1].';
    x=linspace(0,1,m).';
    M=(ones(m,1)+x).*A+B;
    u=M\(b-ones(m,1));
    err(counter)=norm(u-solesatta(x),inf);
end
close all
mrange=[mrange,mrif];
loglog(mrange,err,'*',mrange,err(end)*(mrange/mrange(end)).^(-2));
    