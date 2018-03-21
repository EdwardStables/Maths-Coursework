r=1000;
c=100*10^(-9);
w=20000*pi;
qi = 500*10^(-9);
ti=0.000; tf=0.008;
T = 100*10^(-6); %
f = 1/T;
Vin =@(t) 5*cos(2*pi*f*t);
func=@(t,q) ((Vin(t)-q/c)/r);
x=[4:1:8];
h= 10.^(-x);
error = zeros(1, size(h,2)); 
for j=1:size(h,2)
    h(j)
    [t,y] =RK2(func,qi,[ti tf], h(j));
    Vout = y/c;
    exact= ((5*c*cos(w*t)+5*w*r*c^2*sin(w*t)+(qi*(1+4*pi^2)-5*c)*exp(-t/(r*c)))/((1+(w*r*c)^2)))/(100*10^-9);
    error(j) = max(abs(exact-Vout));
end 
figure
loglog((h),(error))
title('loglog plot of error vs h')
xlabel 'h'
ylabel 'Error'
plot(log(h), log(error));
xlabel 'log(h)'
ylabel 'log(Error)'

