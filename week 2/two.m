clc
clear all
close all 

N = 1024 ;
T = 1 %sample time
e= randn([N,1]) ;

%%
periodog = periodogram(e) ;
omega = [] ;
for i = 1:(N/2+1)
    omega(i,1) = 2*pi*(i-1)/N ;
end

figure(1)
loglog(omega, periodog)

%%
P = tf(1,[1 -0.9 0.5]) ;
k = [0:1:N-1];
ts = T/N ;
kts = k*ts ;
for i = 1:(N)
    omega(i,1) = 2*pi*(i-1)/N ;
end

w = lsim(P, e, kts) ;
period = fft(w) ;
period = period.*(conj(period)) ;
period = period/N ;

figure(2)
loglog(omega, period)


%% new periodogram

period = fft(w) ;
period = period.*(conj(period)) ;

