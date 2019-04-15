
%%
%Question 1.1, Part a
clear; close all
initialCond = [0.15264,-0.02281,0.38127]; %Vc1, Vc2, iL
dt = 0.4; iterations = 17000;
RK4(iterations,dt,initialCond);
%%
%Question 1.2, Part b
clear; close all;
dt = 0.4; iterations = 17000;
initialCond = [2.532735,1.285458e-3,-3.3674];
RK4(iterations,dt,initialCond);

%%
%Question 1.3, initial conditions for which solution does not diverge...
clear; close all;
dt = 0.04; iterations = 1000;
initialCond = [2.532735,1.29,0.2227];
RK4(iterations,dt,initialCond);
%%
%Question 2: The Heat Equation
%Part a
dt = 0.01; dx = 0.1;
x = 0:dx:1;
u = HeatEquation(dt,dx);
g = figure;
subplot(2,2,1);
plot(x,u(1,:)); xlabel('x'); ylabel('u'); title(sprintf('At t = 0') );
subplot(2,2,2);
plot(x,u(3,:)); xlabel('x'); ylabel('u'); title(sprintf('At t = 0.02 '));
subplot(2,2,3);
plot(x,u(5,:)); xlabel('x'); ylabel('u'); title(sprintf('At t = 0.04'));
g.Name = sprintf(' dt = %g  ', dt);
%part b
%%
dt = 0.005; dx = 0.1;
u = HeatEquation(dt,dx);
x = 0:dx:1;
g = figure;
subplot(2,2,1);
plot(x,u(1,:)); xlabel('x'); ylabel('u'); title(sprintf('At t = 0 and dt = %g', dt));
subplot(2,2,2);
plot(x,u(6,:)); xlabel('x'); ylabel('u'); title(sprintf('At t = 0.025 and dt = %g', dt));
subplot(2,2,3);
plot(x,u(11,:)); xlabel('x'); ylabel('u'); title(sprintf('At t = 0.05 and dt = %g', dt));
