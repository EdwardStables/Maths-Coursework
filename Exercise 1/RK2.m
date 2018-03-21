function [xFull, yFull] = RK2(eqn, init, time, h)
%This function returns two vectors, one of the time range used for
%calculations, and the other with the output of the differential equation. 

%eqn: the differential equation to be used, in the form y' = f(x,y)
%init: the value of the output at the first time point (initial condition)
%time: [t0 tf] defining the time range to be used
%h: the step size to be used

%xFull and yFull are defined at the start for efficiency
xFull = (time(1):h:time(2))'; %define the time axis 
yFull = 0*xFull; %yFull is defined in terms of xFull, ensuring that they are of the same size
yFull(1)=init; %defining the initial condition

%-----------------------------------------------------------------------
%UNCOMMENT FOR HEUN:
%a=0.5; b=0.5; p=1; q=1;
%UNCOMMENT FOR MIDPOINT:
%a=0; b=1; p=0.5; q=0.5;
%UNCOMMENT FOR RALSTON:
a=1/4; b=3/4; p=2/3; q=2/3;
%-----------------------------------------------------------------------


for i = 1:(time(2)-time(1))/h %the loop is repeated for N iterations, where N is the number of sections
    xTemp = xFull(i); %xTemp and yTemp are defined for efficency so that xFull and yFull don't need to be constantly accessed
    yTemp = yFull(i);
    k1 = eqn(xTemp, yTemp); %This part defines k1 and k2 to produce yFull
    k2 = eqn(xTemp+p*h, yTemp+h*q*k1);
    yFull(i+1) = yTemp + h*(a*k1+b*k2);
end  


end

