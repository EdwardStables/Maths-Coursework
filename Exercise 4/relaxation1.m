
close all

h = 0.02; %step size
a0 = 0; %initial value of x
b0 = 0; %initial value of y
a = 1; %final value of x
b = 1; %final value of y

rmax = 0.0001; %residue value for which terminates

Nx = round((a - a0)/h)+1; %number of steps in x and y direction
Ny = round((b - b0)/h)+1;



x = a0:h:a; %values of x,y to correct plots
y = b0:h:b;


U = zeros(Nx,Ny); %creates initial matrix of zero

for i = 1:Nx %sets edge values of U equal to boundary conditions defined by p1,p2,p3,p4
    U(1,i) = p3(x(i)); %boundary on x axis
    U(Ny,i) = p1(x(i)); %boundary on y = b
end

for i = 1:Ny
    U(i,1) = p4(y(i)); %boundary on y axis
    U(i,Nx) = p2(y(i)); %boundary on x = a
end
    


%calculate mean value of boundaries

k = 0;

for i = 1:Nx
    k = k + U(i,1) + U(i,Nx);
end

for i = 1:Ny
    k = k + U(1,i) + U(Ny,i);
end

k = k/(2*Nx+2*Ny-4); 

%set interior points equal to mean value of boundaries
for i = 2:Nx-1
    for j = 2:Ny-1
        U(i,j) = k;
    end
end

r = zeros(Nx, Ny); %create matrix for residue

for i = 2:Nx-1
    for j = 2:Ny-1
        r(i,j)=1;
    end
end

%do relaxation method on values inside matrix
while max(max(r)) > rmax %terminates when max value of residue in grid < rmax
    for i = 2:Nx-1
        for j = 2:Ny-1
            Unew = 0.25*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1));
            r(i,j) = abs(Unew - U(i,j)); %calculate residue
            U(i,j) = Unew;  
        end
    end
end

figure
subplot(1,2,1);     %plot equipotential lines of solution
contour(x,y,U);
xlabel("x");
ylabel("y");

subplot(1,2,2);
surf(x,y,U);        %plot solution
xlabel("x");
ylabel("y");

set(gcf, 'Position', [100,100,1200,500]) %resize figure

%boundary conditions defined as functions
%uncomment lines for different conditions
function f = p1(x) %boundary u(x,b)
f = 0; %boundary conditions given in task

% f = -sin(pi*x); %other boundary conditions for testing

% f = 1-x;

% f = 0;
end

function f = p2(y) %boundary u(a,y)
f = 0; %boundary conditions given in task

% f = sin(pi*y); %other boundary conditions for testing

% f = 1-y;

% f = 0;
end

function f = p3(x) %boundary u(x,0)
%     if(x < 0.2||x > 0.8) %boundary conditions given in task
%  	  f = 0;
%     else
%       f = 1;
%     end
%     
% f = -sin(pi*x); %other boundary conditions for testing

% f = x;

f = 1;
end

function f = p4(y) %boundary u(0,y)
 f = 0; %boundary conditions given in task

% f = sin(pi*y); %other boundary conditions for testing

% f = y;

% f = 0;
end




