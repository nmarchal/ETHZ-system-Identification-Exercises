clc
clear all
close all

N = 10^3;
e_u = randn(N,1) ;

%transfer functions
A = [1 -1.5 0.7] ;
B = [0 1 0.5] ;
C = [1 -1 0.2] ;

for i = 1:50
    e_u = randn(N,1) ;
    [theta(:,i),y,u] = parameters_ls(A,B,C,N,e_u) ;
end
%% compare ruth and estimation 

%transfer functions ESTIMATED
A_est = [1 theta(1,1) theta(2,1)] ;
B_est = [0 theta(3,1) theta(4,1)] ;

BA_est = tf(B_est,A_est,-1,'Variable','z^-1') ;
CA_est = tf(C,A_est,-1,'Variable','z^-1') ;

% Create the y output
y_est = lsim(BA_est,u) + lsim(CA_est,e_u);

figure(1)
plot(y) ; hold on ; 
plot(y_est)
legend('true', 'estimation')

%% histogram
figure(2) 
subplot(2,2,1)
histogram(theta(1,:))
title('a1')
subplot(2,2,2)
histogram(theta(2,:))
title('a2')
subplot(2,2,3)
histogram(theta(3,:))
title('b1')
subplot(2,2,4)
histogram(theta(4,:))
title('b2')


%%
function [theta, y,u] = parameters_ls(A,B,C,N, e_u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

BA = tf(B,A,-1,'Variable','z^-1') ;
CA = tf(C,A,-1,'Variable','z^-1') ;
C_inv = tf(1,C,-1,'Variable','z^-1') ;
L = tf([0 1 0.2],[1 -0.1 -0.12], -1,'Variable','z^-1') ;

% construc the u signal
for i = 1:N
    if i == 1
        u(i) = 0 ;
    elseif i == 2
        u(i) = e_u(i-1) ;
    else
        u(i) = 0.1*u(i-1) + 0.12*u(i-2) + e_u(i-1) + 0.2*e_u(i-2) ;
    end
end

% Create the y output
y = lsim(BA,u) + lsim(CA,e_u);

% LS model
y_F = lsim(C_inv,y) ;
u_F = lsim(C_inv,u) ;

for k = 1:N
    if k == 1
        PHI(k,:) = [0 0 0 0] ;
    elseif k == 2
        PHI(k,:) = [-y_F(k-1) 0 u_F(k-1) 0] ;
    else
        PHI(k,:) = [-y_F(k-1) -y_F(k-2) u_F(k-1) u_F(k-2)] ;
    end
end

theta = PHI\y_F ;
end