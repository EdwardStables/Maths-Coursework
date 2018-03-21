%same as relaxation2 but works with err1
%create rmax before calling script

h = 0.01; %step size
a0 = 0; %initial value x,y
a = 1; % final value x,y

N = round((a - a0)/h) + 1; %number of steps 

x = a0:h:a; %values of x and y to correct plot
y = a0:h:a;


%create an array with values of function 
for i = 1:N
    for j = 1:N
        u(j,i) = f(x(i),y(j));
    end
end

%relaxation part

%create NxN matrix for solution
U = zeros(N,N);
%create value to calculate mean
k = 0;

%set edges of U to boundary conditions
for i = 1:N
    U(1,i) = f(x(i),0);
    U(N,i) = f(x(i),a);
    U(i,1) = f(0,y(i));
    U(i,N) = f(a,y(i));
    
    %Sum values of edges for mean calculation
    k = U(1,i) + U(N,i) + U(i,1) + U(i,N) + k;
end

k = k/(4*N);

%set values inside matrix to mean value of boundaries
for i = 2:N-1
    for j = 2:N-1
        U(i,j) = k;
    end
end

r = zeros(N,N); %create matrix for residue

for i = 2:N-1
    for j = 2:N-1
        r(i,j)=1;
    end
end


while max(max(r)) > rmax %only terminates when max value of residue in grid < rmax
   
    for i = 2:N-1  %only for non boundary points in array
        for j = 2:N-1
            
            %calculate new value for U(i,j)
            Unew = 0.25*(U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - g(x(i),y(j))*h^2);
            
            %calculate residue
            r(i,j) = abs(Unew - U(i,j)); 
            
            %replace U(i,j) with its new value
            %if(r(i,j) > rmax)
            U(i,j) = Unew;
            %end
        end
    end
end

%calculate error
err = abs(u - U);


%function 
function a = f(x,y)
 a = x^2 + y + 2*x*y - 1;
end

%forcing term of poisson equation
function m = g(x,y)
m = 2;
end