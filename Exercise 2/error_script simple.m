%Values
R=1000; %Resistor
C=100*10^(-9); %Capacitor
w=20000*pi; %Frequency
qi = 500*10^(-9); %initial charge
h=0.000001; %step value
ti=0.0000; tf=0.0012; %initial and final time
T = 100*10^(-6); %Cosine period
f = 1/T; %Cosine frequency

%Equations
Vin =@(t) 5*cos(2*pi*f*t);
eqn=@(t,Vout) ((Vin(t)-Vout/c)/r);
%To obtain the results for the 3 different RK2 methods, change the
%commented lines in the RK2 function
[t,y] = RK2(eqn, qi, [ti tf], h);
Vout = y/c; %This section deals with charge, therefore the output needs to
%be converted to voltage

exact= ((5*C*cos(w*t)+5*w*R*C^2*sin(w*t)+... %Exact derivation is included in the report
    (qi*(1+4*pi^2)-5*C)*exp(-t/(R*C)))/((1+(w*R*C)^2)))/(100*10^-9);
error=exact-Vout; %Error calculated as the difference between 

%Output plotting
figure;
plot(t, error,'g');
title('Error')
xlabel 't/s'
ylabel 'Vout/V'



    
    
    