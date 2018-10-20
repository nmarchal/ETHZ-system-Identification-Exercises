clc
clear all
close all

r = 5 ;
M = 1024 ;
L = r*M ;

G = tf ([0.1 0],[1 -1.7 0.72],1) ;
H = tf([1 -0.92],[1 -0.5],1) ;
H = H*1.5 ;


%% a
e = 0.01*randn(L,1);
u1 = 2*randn(M,1);
u = [] ;
for i = 1:r
    u = [u ; u1] ;
end
time = [0:1:L-1] ;


y1 = lsim(G,u,time) ;
y2 = lsim(H,e,time) ;
y = y1+y2 ;

%% plots
figure(1)
plot(time,y);

figure(2)
plot(time,u) ;

%% b autocorrelaction

%periodic signal
autocorrelation = [] ;
for tao = -L/2 + 1: L/2 ;
    step = 0 ;
    for k = 1:M
        ktao = k-tao ;
        while (ktao<=0) 
            ktao = ktao + M ;
        end
        step = step + u(k,1)*u(ktao,1) ;
    end
    autocorrelation = [autocorrelation ; 1/M*step] ;
end
figure(3)
timeauto = time-L/2*ones(size(time)) ;
plot(timeauto,autocorrelation) ;

%final energy signal
autocorrelation = [] ;
for tao = -L/2 + 1: L/2 ;
    step = 0 ;
    for k = 1:L
        ktao = k-tao ;
        if (ktao>0 && ktao <= L) 
        step = step + u(k,1)*u(ktao,1) ;
        end
    end
    autocorrelation = [autocorrelation ; 1/M*step] ;
end
figure(4)
plot(timeauto,autocorrelation) ;
% 
% figure(5)
% autocorrelationMatlab = xcorr(u) ;
% plot(timeauto,xcorr(u))

%% c construct ETFE
u_bar = u1 ;
y = reshape(y,M,r) ;
y_bar = 1/(r-1) * sum(y(:,2:end),2) ;

U_bar = fft(u_bar) ;
Y_bar = fft(y_bar) ;

omega = (2*pi/M)*[0:M-1]';
idx = find(omega > 0 & omega < pi);

figure(6)
loglog(omega(idx),abs(U_bar(idx)))
figure(7)
loglog(omega(idx),abs(Y_bar(idx)))

Gest = Y_bar./U_bar;
figure(8)
loglog(omega(idx),abs(Gest(idx)))
figure(9)
semilogx(omega(idx),angle(Gest(idx)))
