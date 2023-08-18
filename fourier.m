%
% equazione
% -u''(x)+u(x) = g(x), -pi <= x < pi
% u(a) = u(b)
% u'(a)=u'(b)
%
close all
clear all
mrif = 256;
m = mrif;
a = -pi;
b = pi;
x = linspace (a, b, m + 1)';
% Butto via l'ultimo punto di "quadratura" (Escluso da ogni calcolo)
x = x(1:m);
g = @(x) 1./(sin (x) + 2);
ghat = fftshift (fft (g(x))) * sqrt (b - a) / m;
%                              aggiungiamo noi l'intervallo
lambda = 1i * 2 * pi * (-m/2:m/2-1)' / (b - a);
% Usare 1i per non sovrascrivere
utilde = ghat./(1-lambda.^2);
u = ifft (ifftshift (utilde)) * m / sqrt(b - a);
% Nota: se mi fossi dimenticato m/sqrt sia sopra che sotto non cambiava niente: l'intervallo di def
% per la valutazione della funzione non importava
% Invece per l'analisi dei coefficienti bisogna ricordarselo, anche per Id di Parseval
% Funzione u già calcolata nei nodi di quadratura
uhatrif = utilde;
% Calcolo dell'ordine
urif = u;
plot(x,u,'*')
pause
%Nota che u(1) =\ u(end): lo avevamo tolto prima
% Possiamo rimetterlo se volessimo
mrange = [8:8:mrif-16];
counter = 0;
for m = mrange
  counter = counter+1;
  x = linspace (a, b, m + 1)';
  x = x(1:m);
  ghat = fftshift (fft (g(x))) * sqrt (b - a) / m;
  lambda = 1i * 2 * pi * (-m/2:m/2-1)' / (b - a);
  utilde = ghat ./ (1 - lambda .^ 2);
  u = ifft (ifftshift (utilde)) * m / sqrt (b - a);
% Corrispondenza fra i coefficienti: devo confrontare i centrati (freq neg e freq positive)
  utilde = [zeros((mrif - m) / 2, 1); utilde; zeros((mrif - m) / 2, 1)];
% Parseval: norma due della differenza dei coefficienti
  errore(counter) = norm (uhatrif-utilde);
end
% Mostriamo la convergenza spettrale: vediamo che converge più velocemente
% di ogni ordine
% con 48 punti raggiungo la precisione di macchina!
% Perciò il problema avrà più derivate periodiche 
loglog(mrange, errore(1) * (mrange / mrange(1)) .^ (-1),'k',...
mrange,errore(1) * (mrange / mrange(1)) .^ (-2),'r',...
mrange,errore(1) * (mrange / mrange(1)) .^ (-4),'m',...
mrange,errore(1) * (mrange / mrange(1)) .^ (-8),'g',mrange, errore, '*b')
xlabel('m')
ylabel('errore in norma-2')
legend('ordine 1','ordine 2','ordine 4','ordine 8','errore',...
'location','NorthWest')
%!demo
%! m = 64;
%! a = -pi;
%! b = pi;
%! g = @(x) 1 ./ (sin (x) + 2);
%! x = linspace (a, b, m + 1)';
%! x = x(1:m);
%! ghat = fftshift (fft (g(x))) * sqrt (b - a) / m;
%! phi = @(j, x) exp(1i*(j-1-m/2)*2*pi*(x-a)/(b-a)) / sqrt(b-a);
%! ghatx = zeros (size (x));
%! for j = 1:m
%!   ghatx = ghatx + ghat(j) * phi(j,x);
%! end
%! ghathat = ifft (ifftshift (ghat)) / sqrt (b - a) * m;
%! plot(x,g(x),x,ghatx,'*',x,ghathat,'o')
%! legend('g(x)','ghat(x)','ghathat')
%! disp('ghathat interpola: massimo errore')
%! norm(ghathat-g(x),inf)
%!demo
%! m = 32;
%! a = -pi;
%! b = pi;
%! g = @(x) 1 ./ (sin (x) + 2);
%! x = linspace (a, b, m + 1)';
%! x = x(1:m);
%! ghat = fftshift (fft (g(x))) * sqrt (b - a) / m;
%! phi = @(j, x) exp(1i*(j-1-m/2)*2*pi*(x-a)/(b-a)) / sqrt(b-a);
%! xtarget = (b - a) * rand (30, 1) + a;
%! ghatxtarget = zeros (size (xtarget));
%! for j = 1:m
%!   ghatxtarget = ghatxtarget + ghat(j) * phi(j,xtarget);
%! end
%! plot(x,g(x),xtarget,real(ghatxtarget),'*')
%! legend('g(x)','ghat(xtarget)')
%! disp('ghat(xtarget) NON interpola: massimo errore')
%! norm(ghatxtarget-g(xtarget),inf)
%! disp('ghat(xtarget) NON reale')
%! ghatxtarget
