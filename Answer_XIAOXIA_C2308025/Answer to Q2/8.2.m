clear
clc

%% Initialization
h = 0.01;  % Step size
t = 0:h:30;
x = zeros(length(t),2);  
x(1,:) = [0,10];  % The original function and the initial derivative
% Call function
X = FourRK(t,h,x); 

P = 10 - (X(:,1) + X(:,2));
E = 1 - X(:,1);

%% plot
figure(1)
plot(t,X(:,1),'Linewidth',2)
hold on
plot(t,X(:,2),'Linewidth',2) 
hold on
plot(t,P,'Linewidth',2)  
hold on
plot(t,E,'Linewidth',2)  
hold off
xlabel('time(s)')
ylabel('Concerntration(¦ÌM)')
legend('ES','S','P','E')

figure(2)
plot(t(1:end-1),diff(X(:,1))/h,'Linewidth',2)  
hold on
plot(t(1:end-1),diff(X(:,2))/h','Linewidth',2) 
hold on
plot(t(1:end-1),diff(P)/h,'Linewidth',2)  
hold on
plot(t(1:end-1),diff(E)/h,'Linewidth',2) 
hold off
xlabel('time(s)')
ylabel('Concerntration change rate(¦ÌM/s)')
legend('dES/dt','dS/dt','dP/dt','dE/dt')

%% Call function

% the fourth-order Runge-Kutta method
function x = FourRK(t,h,x)
for i = 2:length(t) 
K1 = Fun(t(i-1),x(i-1,:)); % K1
K2 = Fun(t(i-1)+1/2*h , x(i-1,:)+1/2*h.*K1); % K2
K3 = Fun(t(i-1)+1/2*h , x(i-1,:)+1/2*h.*K2);  % K3
K4 = Fun(t(i-1)+h , x(i-1,:)+h.*K3);  % K4
x(i,:) = x(i-1,:) + h/6.*(K1 + 2*K2 + 2*K3 + K4);
end
end

% Target Function
function dx = Fun(t,x)
k1 = 100/60; k2 = 600/60; k3 = 150/60;
e0 = 1;

% es = x(1)  s = x(2)
dx(1) = k1*(e0 - x(1))*x(2) - (k2 + k3)*x(1);  
% The first differential function
dx(2) = -k1*(e0 - x(1))*x(2) + k2*x(1);  
% The second differential function
end