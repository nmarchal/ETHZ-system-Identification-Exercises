clc
clear all
close all

load('Data_ex9.mat')

B_PLOT_Y_PRED = true ; %plot exercice 1

global i; i=1;
%will count number of iterations

global u;
u = p2_u;
y = p2_y;

% Signals periods and lengths
M = length(p2_u);
global N; N = length(p2_u);

% Problem parameters
global NB_PARAMS; NB_PARAMS = 4;
theta = ones(NB_PARAMS,1); %(a1,a2,b1,c1)



%% prob 1

a1=theta(1);a2=theta(2);b1=theta(3);c1=theta(4); 
w = zeros(N,1);
w(1) = 0;
w(2) = -a1*w(1)-0+b1*u(1);
for i = 3:N
    w(i) = -a1*w(i-1) - a2*w(i-2) + b1*u(i-1);
end

% Create data part of regressor. Assume plant at rest
e = zeros(N,1); % computed iteratively
PhiTyu(1,:) = [0    , 0,    0,      0]; % [-w(k-1) -w(k-2) u(k-1) e(k-1)]
PhiTyu(2,:) = [-w(1), 0, u(1), e(1)];
for i = 3:N
    PhiTyu(i,:) = [-w(i-1), -w(i-2), u(i-1), e(i-1)];
end    % y^ =  [-w(k-1) -w(k-2) u(k-1) e(k-1)] * [a1,a2,b1,c1]^T;

% Minimize error and optimize theta
x0 = zeros(NB_PARAMS+N,1);      % [ [a1,a2,b1,c1], e ]

[x,fval] = fmincon(@(x)ARMAXobjective(x),x0,...
        [],[],[],[],[],[],@(x)ARMAXconstraint(x,y,PhiTyu));

% Predicted values
theta = x(1:NB_PARAMS); %(a1,a2,b1,c1)
a1 = theta(1); a2 = theta(2); b1 = theta(3); c1 = theta(4);
e = x(NB_PARAMS+1:end);

% create new w(k)
w = zeros(N,1);
w(1) = 0;
w(2) = -a1*w(1)-0+b1*u(1);
for i = 3:N
    w(i) = -a1*w(i-1) - a2*w(i-2) + b1*u(i-1);
end

%create new 
Phi = zeros(N,NB_PARAMS);
Phi(1,:) = [0    , 0,    0,      0]; % [-w(k-1) -w(k-2) u(k-1) e(k-1)]
Phi(2,:) = [-w(1), 0, u(1), e(1)];
for i = 3:N
    Phi(i,:) = [-w(i-1), -w(i-2), u(i-1), e(i-1)];
end    % y^ =  [-w(k-1)  -w(k-2)  u(k-1)  e(k-1)] * [a1,a2,b1,c1]^T;

ypred = Phi * theta;   
e = y - ypred;
    
if B_PLOT_Y_PRED
    % Plot prediction
    figure ;
    plot(e)
    hold on; grid on;
%     plot(ypred)
    plot(y)
%     legend('y', 'ypred', 'e');
    legend('error',  'Original signal');
    
    figure ;
    plot(y) ; hold on ; plot(ypred) ;
    grid on ;
    legend('signal','predicted signal') ;
    
    figure
    plot(y) ; hold on ; plot(ypred) ; plot(ypred1) ;
    grid on ;
    legend('signal','predicted signal','1st predicted signal') ;
    
end 

%% part 2

%% least square estimate
NB_PARAMS = 3;
% theta_LS = ones(NB_PARAMS,1); %(a1,a2,b1)

Phi = zeros(N,NB_PARAMS);
Phi(1,:) = [0    , 0,    0]; % [-y(k-1) -y(k-2) u(k-1)]
Phi(2,:) = [-y(1), 0, u(1)];
for i = 3:N
    Phi(i,:) = [-y(i-1), -y(i-2), u(i-1)];
end    % y^ =  [-y(i-1), -y(i-2), u(i-1)] * [a1,a2,b1]^T;

theta_LS = Phi \ y;

%% instrumntal variable
    
% Define instruments
Phi = zeros(N,NB_PARAMS);
zeta = zeros(N,NB_PARAMS);
Phi(1,:) = [0    , 0,    0];  % [-y(k-1) -y(k-2) u(k-1)]
Phi(2,:) = [-y(1), 0, u(1)];
zeta(1,:) = [0    , 0,    0]; % [-x(k-1) -x(k-2) u(k-1)]
zeta(2,:) = [-x(1), 0, u(1)];
for i = 3:N
    Phi(i,:) = [-y(i-1), -y(i-2), u(i-1)];
    zeta(i,:) = [-x(i-1), -x(i-2), u(i-1)];
end 

% Define stuff for final estimate using instruments
R = zeros(NB_PARAMS,NB_PARAMS);
f_k = zeros(NB_PARAMS,1);
for k=1:N
    R   = R   + zeta(k,:)' * Phi(k,:);
    f_k = f_k + zeta(k,:)' * y(k);
end
R   = 1/N * R;   % note that it is invertible! :D
f_k = 1/N * f_k;

% Estimate using instruments (formula given in lecture)
theta_instr = R \ f_k;


%% Functions
function [error] = ARMAXobjective(x) % x = [theta; e]
    global NB_PARAMS;
    global i;
    
    i=i+1   ;
    error = sqrt(x(NB_PARAMS+1:end)'*x(NB_PARAMS+1:end)); 
end

function [c,ceq] = ARMAXconstraint(x,y,PhiTyu)
    global N;
    global NB_PARAMS;
    global u;
    
    theta = x(1:NB_PARAMS);
    a1 = theta(1); a2 = theta(2); b1 = theta(3); c1 = theta(4);
    e = x(NB_PARAMS+1:end);
    
    % create w(k)
    w = zeros(N,1);
    w(1) = 0;
    w(2) = -a1*w(1)-0+b1*u(1);
    for i = 3:N
        w(i) = -a1*w(i-1) - a2*w(i-2) + b1*u(i-1);
    end
    
    Phi = zeros(N,NB_PARAMS);
    Phi(1,:) = [0    , 0,    0,      0]; % [-w(k-1) -w(k-2) u(k-1) e(k-1)]
    Phi(2,:) = [-w(1), 0, u(1), e(1)];
    for i = 3:N
        Phi(i,:) = [-w(i-1), -w(i-2), u(i-1), e(i-1)];
    end    % y^ =  [-w(k-1)  -w(k-2)  u(k-1)  e(k-1)] * [a1,a2,b1,c1]^T;
        
    ceq = y - Phi * theta - e;
    c = [];
end

