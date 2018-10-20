% same as before without periodic signal

clc
clear all
close all

G = tf([0.1 0], conv([1 -1.7 0.72],[1 -0.98 0.9]),1) ;
H = tf([0.5 -0.9*0.5],[1 -0.25],1) ;

%generate signals
N=1024 ;
m = 4 ;
time = [0:1:N-1] ;
e = randn(N,1) ;
u = 1 + 2*randn(N,1) ;

ye = lsim(H,e,time) ;
yu = lsim(G,u,time) ;

y = ye + yu ;

figure(1)
title('output y') ; hold on ;
plot(time,y); 
figure(2)
title('input u') ; hold on ;
plot(time,u);

Y = fft(y) ;
U = fft(u) ;

figure(3)
title('DTFT of input') ; hold on ;
stem(U)

G_bar = abs(Y./U) ;

omega = exp(i*2*pi/N*[0:1:(N-1)]) ; 
Gfreq = squeeze(abs(freqresp(G,omega))) ;

figure(4)
loglog( G_bar); hold on ;
loglog( Gfreq) ;
title('unsmoothened plant estimation') ;
legend('estimated tf', 'real tf') ;
ylim([10^(-2) 10^(2)]) ;

figure(5)
loglog( abs(G_bar-Gfreq)) ; hold on ;
loglog( Gfreq) ;
title('error before smoothening') ;
legend('error', 'real tf') ;

%% c) do the average

omega_avg = exp(i*2*pi/N*[0:1:(N-1)]) ;
omega_avg = reshape(omega_avg , [4,N/4]) ; 
omega_avg = omega_avg(1,:)';
Gfreq_avg = squeeze(abs(freqresp(G,omega_avg))) ;
time_avg = [0:1:(N/m - 1)] ;
    
G_bar_average_tot = [] ;
for i = 1:m
    e_avg = randn(N/m,1) ;
    u_avg = 1 + 2*randn(N/m,1) ;
    ye_avg = lsim(H,e_avg,time_avg) ;
    yu_avg = lsim(G,u_avg,time_avg) ;
    y_avg = ye_avg + yu_avg ;
    Y_avg = fft(y_avg) ;
    U_avg = fft(u_avg) ;
    G_bar_average_tot = [G_bar_average_tot abs(Y_avg./U_avg)] ;
end

G_bar_average = zeros(N/m,1)' ;
for i = 1:m ;
    G_bar_average = G_bar_average + 1/m * G_bar_average_tot(((i-1)*N/m +1):i*N/m) ;
end

figure(6)
loglog( G_bar_average); hold on ;
loglog( Gfreq_avg) ;
title('averaged plant estimation') ;
legend('estimated tf', 'real tf') ;
ylim([10^(-2) 10^(2)]) ;

%% d) hann response (frequency)
gama = [5 10 50 100] ;
figure(7) ; hold on ;
subplot(2,1,1);
xlim([-2 2]) ;
for i = 1:size(gama,2) 
    [omega,WHann] = WfHann(gama(i),N) ;
    plot(omega, WHann) ; 
end
legend ('5', '10', '50', '100') ;

%% e) hann response (time)
gama = [5 10 50 100] ;
subplot(2,1,2);
for i = 1:size(gama,2) 
    [lags,WHann] = WtHann(gama(i),N) ;
    plot(lags, WHann) ; 
end
legend ('5', '10', '50', '100') ;

%% f) smoothed signal
figure(8)
for i = 1:size(gama,2) 
    Gest = G_bar ;
    Gs = 0*Gest; % smoothed estimate
    [omega,Wg] = WfHann(gama(i),N); % window (centered)
    zidx = find(omega==0); % shift to start at zero
    omega = [omega(zidx:N);omega(1:zidx-1)]; % frequency grid
    Wg = [Wg(zidx:N) Wg(1:zidx-1)];
    a = U.*conj(U); % variance weighting
    for wn = 1:N,
        Wnorm = 0; % reset normalisation
        for xi = 1:N,
            widx = mod(xi-wn,N)+1; % wrap window index
            Gs(wn) = Gs(wn) + ...
            Wg(widx) * Gest(xi) * a(xi);
            Wnorm = Wnorm + Wg(widx) * a(xi);
        end
    Gs(wn) = Gs(wn)/Wnorm; % weight normalisation
    end
    subplot(2,2,i)
    loglog(Gs) ; hold on ;
    loglog( Gfreq) ;
    legend(int2str(gama(i)),'real')
end

title('frequency smoothened');
