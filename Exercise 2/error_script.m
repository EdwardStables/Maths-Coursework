%circuit values, time values, and initial conditions
r=1000;
c=100*10^(-9);
w=20000*pi;
qi = 500*10^(-9);
ti=0.000; tf=0.008;
T = 100*10^(-6); 
f = 1/T;

%forming the equations, and an array of h values to be tested
Vin =@(t) 5*cos(2*pi*f*t);
func=@(t,q) ((Vin(t)-q/c)/r);
x=[5:1:9];
h= 10.^(-x);
error = zeros(1, size(h,2)); 
for j=1:size(h,2)
    %Get the RK2 output of the equation
    [t,y] =RK2(func,qi,[ti tf], h(j));
    Vout = y/c;
    %form the exact equation from the same time values
    exact= ((5*c*cos(w*t)+5*w*r*c^2*sin(w*t)+(qi*(1+4*pi^2)-5*c)*exp(-t/(r*c)))/((1+(w*r*c)^2)))/(100*10^-9);
    %add the error value to the error array
    error(j) = max(abs(exact-Vout));
end 
figure
loglog((h),(error))
title('Error vs h (Ralston)')
xlabel 'h'
ylabel 'Error'
figure;
plot(log(h), log(error));
xlabel 'log(h)'
ylabel 'log(Error)'
title('Loglog of Error vs h (Ralston)')

