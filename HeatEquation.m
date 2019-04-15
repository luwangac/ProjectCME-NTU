function U = HeatEquation(dt,dx)
mu = dt/(dx^2); 
cols = length([0:dx:1]);
rows =  length([0:dt:1]);
u = zeros(rows,cols); %u(t,x) 1st row is t = 0, second row is t = 1*dt, etc
%Boundary Conditions
t = 0:dt:1;
u(:,1) = 0; u(:,end) = 0;
tol = 1d-6;
x = 0:dx:1;
%Initial conditions (i.e. when t = 0, which u(1,x))
u(1,1:ceil((1/5)*(length(x)))) = -1*(x(1:find(abs(x-1/5)<tol)));
u(1,ceil((1/5)*(length(x)))+1: ceil((7/10)*(length(x)))) = x(find(abs(x-1/5)<tol)+1: find(abs(x-7/10)<tol))-2/5;
u(1,ceil((7/10)*(length(x)))+1: end) = 1-x(find(abs(x-7/10)<tol)+1:end);
for xk = 2:length(x)-1
    for ti = 1:length(t)
        u(ti+1,xk) = mu*u(ti,xk+1) + (1-2*mu)*u(ti,xk)+ mu*u(ti,xk-1);

    end
end
U = u;
end