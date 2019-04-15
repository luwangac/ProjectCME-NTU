%%
iterations = 10000;
dt = 0.04; h =dt;
t = 0:dt:iterations*dt;
%%
%Runge-Kutta Method
%The Constants and functions

%THE CHAOTIC ATTRACTOR

C1 = 1/9; C2 = 1; L = 1/7; G = 0.7; m0 = -0.5; m1 = -0.8; Bp = 1;
g = @(v)m0*v+0.5*(m1-m0)*abs(v+Bp)+0.5*(m0-m1)*abs(v-Bp);

%the ODEs
Vc1dot = @(t,Vc1,Vc2,iL)(1/C1)*G*(Vc2-Vc1) + g(Vc1);
Vc2dot = @(t,Vc1,Vc2,iL)(1/C2)*G*(Vc1-Vc2) + iL;
iLdot = @(t,Vc1,Vc2,iL)(-1/L)*Vc2;
f1 = Vc1dot; f2 = Vc2dot; f3 = iLdot;
%initialize vectors
Vc1 =  zeros(length(t),1); Vc2 =  zeros(length(t),1); 
iL =  zeros(length(t),1); 
%the initial conditions
Vc1(1) = 0.15264; Vc2(1) = -0.02281; iL(1) = 0.38127;

for i = 1:iterations
    %get gradient at t0; i.e. plug into the dot equations above
    k1Vc1 = f1(t(i),Vc1(i),Vc2(i),iL(i));
    k1Vc2 = f2(t(i),Vc1(i),Vc2(i),iL(i));
    k1iL = f3(t(i),Vc1(i),Vc2(i),iL(i));
    %first mid-point value estimate using gradient of t = t0;
    Vc1m1 = Vc1(i)+0.5*h*k1Vc1;
    Vc2m1 = Vc2(i)+0.5*h*k1Vc2;
    iLm1 = iL(i)+0.5*h*k1iL;
    %use the newly-obtained function value to estimate gradient at midpoint
    k2Vc1 = f1(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    k2Vc2 = f2(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    k2iL = f3(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    %second mid-point function value estimate using new gradient
    Vc1m2 = Vc1(i)+0.5*h*k2Vc1;
    Vc2m2 = Vc2(i)+0.5*h*k2Vc2;
    iLm2 = iL(i)+0.5*h*k2iL;
    %Get a more accurate gradient estimate using the new function value
    k3Vc1 = f1(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    k3Vc2 = f2(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    k3iL = f3(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    %estimate the value of function at t = t0+h using the midpoint gradient
    Vc1he =  Vc1(i)+ h*k3Vc1;
    Vc2he = Vc2(i)+ h*k3Vc2;
    iLhe = iL(i)+ h*k3iL;
    %use the estimate of the function at  t = t0+h to estimate gradient at
    %t =t0+h
    k4Vc1 = f1(t(i)+h,Vc1he,Vc2he,iLhe);
    k4Vc2 = f2(t(i)+h,Vc1he,Vc2he,iLhe);
    k4iL = f3(t(i)+h,Vc1he,Vc2he,iLhe);
    %use a weighted sum as the overall gradient
    Vc1(i+1) = Vc1(i)+(k1Vc1/6+k2Vc1/3+k3Vc1/3+k4Vc1/6)*h;
    Vc2(i+1) = Vc2(i)+(k1Vc2/6+k2Vc2/3+k3Vc2/3+k4Vc2/6)*h;
    iL(i+1) = iL(i)+(k1iL/6+k2iL/3+k3iL/3+k4iL/6)*h;

end
subplot(2,2,1)
plot(iL,Vc1); xlabel('iL'); ylabel('Vc1');
subplot(2,2,2);  
plot(iL,Vc2); xlabel('iL'); ylabel('Vc2');
subplot(2,2,3); 
plot(Vc1,Vc2);xlabel('Vc1'); ylabel('Vc2');
subplot(2,2,4)
plot(t,Vc1);xlabel('t'); ylabel('Vc1');
%%
%Runge-Kutta Method
%The Constants and functions

% THE SADDLE-TYPE PERIODIC ORBIT

C1 = 1/9; C2 = 1; L = 1/7; G = 0.7; m0 = -0.5; m1 = -0.8; Bp = 1;
g = @(v)m0*v+0.5*(m1-m0)*abs(v+Bp)+0.5*(m0-m1)*abs(v-Bp);

%the ODEs
Vc1dot = @(t,Vc1,Vc2,iL)(1/C1)*G*(Vc2-Vc1) + g(Vc1);
Vc2dot = @(t,Vc1,Vc2,iL)(1/C2)*G*(Vc1-Vc2) + iL;
iLdot = @(t,Vc1,Vc2,iL)(-1/L)*Vc2;
f1 = Vc1dot; f2 = Vc2dot; f3 = iLdot;
%initialize vectors
Vc1 =  zeros(length(t),1); Vc2 =  zeros(length(t),1); 
iL =  zeros(length(t),1); 
%the initial conditions
Vc1(1) = 2.532735; Vc2(1) = 1.285458e-3; iL(1) = -3.367482;

for i = 1:length(t)-1
    %get gradient at t0; i.e. plug into the dot equations above
    k1Vc1 = f1(t(i),Vc1(i),Vc2(i),iL(i));
    k1Vc2 = f2(t(i),Vc1(i),Vc2(i),iL(i));
    k1iL = f3(t(i),Vc1(i),Vc2(i),iL(i));
    %first mid-point value estimate using gradient of t = t0;
    Vc1m1 = Vc1(i)+0.5*h*k1Vc1;
    Vc2m1 = Vc2(i)+0.5*h*k1Vc2;
    iLm1 = iL(i)+0.5*h*k1iL;
    %use the newly-obtained function value to estimate gradient at midpoint
    k2Vc1 = f1(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    k2Vc2 = f2(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    k2iL = f3(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    %second mid-point function value estimate using new gradient
    Vc1m2 = Vc1(i)+0.5*h*k2Vc1;
    Vc2m2 = Vc2(i)+0.5*h*k2Vc2;
    iLm2 = iL(i)+0.5*h*k2iL;
    %Get a more accurate gradient estimate using the new function value
    k3Vc1 = f1(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    k3Vc2 = f2(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    k3iL = f3(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    %estimate the value of function at t = t0+h using the midpoint gradient
    Vc1he =  Vc1(i)+ h*k3Vc1;
    Vc2he = Vc2(i)+ h*k3Vc2;
    iLhe = iL(i)+ h*k3iL;
    %use the estimate of the function at  t = t0+h to estimate gradient at
    %t =t0+h
    k4Vc1 = f1(t(i)+h,Vc1he,Vc2he,iLhe);
    k4Vc2 = f2(t(i)+h,Vc1he,Vc2he,iLhe);
    k4iL = f3(t(i)+h,Vc1he,Vc2he,iLhe);
    %use a weighted sum as the overall gradient
    Vc1(i+1) = Vc1(i)+(k1Vc1/6+k2Vc1/3+k3Vc1/3+k4Vc1/6)*h;
    Vc2(i+1) = Vc2(i)+(k1Vc2/6+k2Vc2/3+k3Vc2/3+k4Vc2/6)*h;
    iL(i+1) = iL(i)+(k1iL/6+k2iL/3+k3iL/3+k4iL/6)*h;

end
figure
subplot(2,2,1)
plot(iL,Vc1); xlabel('iL'); ylabel('Vc1');
subplot(2,2,2);  
plot(iL,Vc2); xlabel('iL'); ylabel('Vc2');
subplot(2,2,3); 
plot(Vc1,Vc2);xlabel('Vc1'); ylabel('Vc2');
subplot(2,2,4)
plot(t,Vc1);xlabel('t'); ylabel('Vc1');
%%
%Runge-Kutta Method
%The Constants and functions

% LARGE PERIODIC ATTRACTOR

C1 = 1/9; C2 = 1; L = 1/7; G = 0.7; m0 = -0.5; m1 = -0.8; Bp = 1;
g = @(v)m0*v+0.5*(m1-m0)*abs(v+Bp)+0.5*(m0-m1)*abs(v-Bp);

%the ODEs
Vc1dot = @(t,Vc1,Vc2,iL)(1/C1)*G*(Vc2-Vc1) + g(Vc1);
Vc2dot = @(t,Vc1,Vc2,iL)(1/C2)*G*(Vc1-Vc2) + iL;
iLdot = @(t,Vc1,Vc2,iL)(-1/L)*Vc2;
f1 = Vc1dot; f2 = Vc2dot; f3 = iLdot;
%initialize vectors
Vc1 =  zeros(length(t),1); Vc2 =  zeros(length(t),1); 
iL =  zeros(length(t),1); 
%the initial conditions
Vc1(1) = -3.08832; Vc2(1) = -1.0423; iL(1) = 6.93155;

for i = 1:length(t)-1
    %get gradient at t0; i.e. plug into the dot equations above
    k1Vc1 = f1(t(i),Vc1(i),Vc2(i),iL(i));
    k1Vc2 = f2(t(i),Vc1(i),Vc2(i),iL(i));
    k1iL = f3(t(i),Vc1(i),Vc2(i),iL(i));
    %first mid-point value estimate using gradient of t = t0;
    Vc1m1 = Vc1(i)+0.5*h*k1Vc1;
    Vc2m1 = Vc2(i)+0.5*h*k1Vc2;
    iLm1 = iL(i)+0.5*h*k1iL;
    %use the newly-obtained function value to estimate gradient at midpoint
    k2Vc1 = f1(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    k2Vc2 = f2(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    k2iL = f3(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    %second mid-point function value estimate using new gradient
    Vc1m2 = Vc1(i)+0.5*h*k2Vc1;
    Vc2m2 = Vc2(i)+0.5*h*k2Vc2;
    iLm2 = iL(i)+0.5*h*k2iL;
    %Get a more accurate gradient estimate using the new function value
    k3Vc1 = f1(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    k3Vc2 = f2(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    k3iL = f3(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    %estimate the value of function at t = t0+h using the midpoint gradient
    Vc1he =  Vc1(i)+ h*k3Vc1;
    Vc2he = Vc2(i)+ h*k3Vc2;
    iLhe = iL(i)+ h*k3iL;
    %use the estimate of the function at  t = t0+h to estimate gradient at
    %t =t0+h
    k4Vc1 = f1(t(i)+h,Vc1he,Vc2he,iLhe);
    k4Vc2 = f2(t(i)+h,Vc1he,Vc2he,iLhe);
    k4iL = f3(t(i)+h,Vc1he,Vc2he,iLhe);
    %use a weighted sum as the overall gradient
    Vc1(i+1) = Vc1(i)+(k1Vc1/6+k2Vc1/3+k3Vc1/3+k4Vc1/6)*h;
    Vc2(i+1) = Vc2(i)+(k1Vc2/6+k2Vc2/3+k3Vc2/3+k4Vc2/6)*h;
    iL(i+1) = iL(i)+(k1iL/6+k2iL/3+k3iL/3+k4iL/6)*h;

end
figure
subplot(2,2,1)
plot(iL,Vc1); xlabel('iL'); ylabel('Vc1');
subplot(2,2,2);  
plot(iL,Vc2); xlabel('iL'); ylabel('Vc2');
subplot(2,2,3); 
plot(Vc1,Vc2);xlabel('Vc1'); ylabel('Vc2');
subplot(2,2,4)
plot(t,Vc1);xlabel('t'); ylabel('Vc1');


%%
%Runge-Kutta Method
%The Constants and functions

% LARGE PERIODIC ATTRACTOR

C1 = 1/9; C2 = 1; L = 1/7; G = 0.7; m0 = -0.5; m1 = -0.8; Bp = 1;
g = @(v)m0*v+0.5*(m1-m0)*abs(v+Bp)+0.5*(m0-m1)*abs(v-Bp);

%the ODEs
Vc1dot = @(t,Vc1,Vc2,iL)(1/C1)*G*(Vc2-Vc1) + g(Vc1);
Vc2dot = @(t,Vc1,Vc2,iL)(1/C2)*G*(Vc1-Vc2) + iL;
iLdot = @(t,Vc1,Vc2,iL)(-1/L)*Vc2;
f1 = Vc1dot; f2 = Vc2dot; f3 = iLdot;
%initialize vectors
Vc1 =  zeros(length(t),1); Vc2 =  zeros(length(t),1); 
iL =  zeros(length(t),1); 
%the initial conditions
Vc1(1) = -3.08832; Vc2(1) = -1.0423; iL(1) = 6.93155;

for i = 1:length(t)-1
    %get gradient at t0; i.e. plug into the dot equations above
    k1Vc1 = f1(t(i),Vc1(i),Vc2(i),iL(i));
    k1Vc2 = f2(t(i),Vc1(i),Vc2(i),iL(i));
    k1iL = f3(t(i),Vc1(i),Vc2(i),iL(i));
    %first mid-point value estimate using gradient of t = t0;
    Vc1m1 = Vc1(i)+0.5*h*k1Vc1;
    Vc2m1 = Vc2(i)+0.5*h*k1Vc2;
    iLm1 = iL(i)+0.5*h*k1iL;
    %use the newly-obtained function value to estimate gradient at midpoint
    k2Vc1 = f1(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    k2Vc2 = f2(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    k2iL = f3(t(i)+0.5*h,Vc1m1,Vc2m1,iLm1);
    %second mid-point function value estimate using new gradient
    Vc1m2 = Vc1(i)+0.5*h*k2Vc1;
    Vc2m2 = Vc2(i)+0.5*h*k2Vc2;
    iLm2 = iL(i)+0.5*h*k2iL;
    %Get a more accurate gradient estimate using the new function value
    k3Vc1 = f1(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    k3Vc2 = f2(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    k3iL = f3(t(i)+0.5*h,Vc1m2,Vc2m2,iLm2);
    %estimate the value of function at t = t0+h using the midpoint gradient
    Vc1he =  Vc1(i)+ h*k3Vc1;
    Vc2he = Vc2(i)+ h*k3Vc2;
    iLhe = iL(i)+ h*k3iL;
    %use the estimate of the function at  t = t0+h to estimate gradient at
    %t =t0+h
    k4Vc1 = f1(t(i)+h,Vc1he,Vc2he,iLhe);
    k4Vc2 = f2(t(i)+h,Vc1he,Vc2he,iLhe);
    k4iL = f3(t(i)+h,Vc1he,Vc2he,iLhe);
    %use a weighted sum as the overall gradient
    Vc1(i+1) = Vc1(i)+(k1Vc1/6+k2Vc1/3+k3Vc1/3+k4Vc1/6)*h;
    Vc2(i+1) = Vc2(i)+(k1Vc2/6+k2Vc2/3+k3Vc2/3+k4Vc2/6)*h;
    iL(i+1) = iL(i)+(k1iL/6+k2iL/3+k3iL/3+k4iL/6)*h;

end
figure
subplot(2,2,1)
plot(iL,Vc1); xlabel('iL'); ylabel('Vc1');
subplot(2,2,2);  
plot(iL,Vc2); xlabel('iL'); ylabel('Vc2');
subplot(2,2,3); 
plot(Vc1,Vc2);xlabel('Vc1'); ylabel('Vc2');
subplot(2,2,4)
plot(t,Vc1);xlabel('t'); ylabel('Vc1');
