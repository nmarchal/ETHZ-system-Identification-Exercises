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
uperiod = 1+ 2*randn(N/4,1) ;
u = [uperiod ; uperiod ; uperiod ; uperiod] ;

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
G_bar = reshape(G_bar , [4,N/4]) ;
G_bar = G_bar(1,:)';

omega = exp(i*2*pi/N*[0:1:(N-1)]) ;
omega = reshape(omega , [4,N/4]) ; 
omega = omega(1,:)';
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
G_bar_average_tot = [] ;
for i = 1:m
    e_avg = randn(N/m,1) ;
    time_avg = [0:1:(N/m - 1)] ;
    ye_avg = lsim(H,e_avg,time_avg) ;
    yu_avg = lsim(G,uperiod,time_avg) ;
    y_avg = ye_avg + yu_avg ;
    Y_avg = fft(y_avg) ;
    U_avg = fft(uperiod) ;
    G_bar_average_tot = [G_bar_average_tot abs(Y_avg./U_avg)] ;
end

G_bar_average = zeros(N/m,1)' ;
for i = 1:m ;
    G_bar_average = G_bar_average + 1/m * G_bar_average_tot(((i-1)*N/m +1):i*N/m) ;
end

figure(6)
loglog( G_bar_average); hold on ;
loglog( Gfreq) ;
title('averaged plant estimation') ;
legend('estimated tf', 'real tf') ;
ylim([10^(-2) 10^(2)]) ;

%% d) hann response (frequency)
gama = [5 10 50 100] ;
figure(7) ; hold on ;
xlim([-2 2]) ;
for i = 1:size(gama,2) 
    [omega,WHann] = WfHann(gama(i),N) ;
    plot(omega, WHann) ; 
end
legend ('5', '10', '50', '100') ;

%% e) hann response (time)
gama = [5 10 50 100] ;
figure(8) ; hold on ;
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
    [omega,Wg] = WfHann(gama(i),N/m); % window (centered)
    zidx = find(omega==0); % shift to start at zero
    omega = [omega(zidx:(N/m));omega(1:zidx-1)]; % frequency grid
    Wg = [Wg(zidx:(N/m)) Wg(1:zidx-1)];
    a = U.*conj(U);
    % variance weighting
    for wn = 1:(N/m),
        Wnorm = 0; % reset normalisation
        for xi = 1:(N/m),
            widx = mod(xi-wn,N/m)+1; % wrap window index
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
        title('frequency domain window')
end


%% new this week : smoothed signal (time)
figure (9) ;
for i = 1:size(gama,2) 
    j = 0 ;
    [lags,WHann] = WtHann(gama(i),N) ;
    autocorrelation_u = 0*WHann ;
    crosscorrelation_u = 0*WHann ;
    for tao = -(N/2 -1):(N/2)
        j = j+1 ;
        autocorrelation = 0 ;
        for k = 1:size(lags,1)
            ktao = k-tao ;
            while (ktao<=0) 
                ktao = ktao + N ;
            end
            while (ktao > N)
                ktao = ktao - N ;
            end
            autocorrelation = autocorrelation + u(k)*u(ktao) ;
        end
        autocorrelation_u(j) = autocorrelation ;
    end
    denominator = fft( WHann .* autocorrelation_u) ;

    j = 0 ; %reset coutner
    for tao = -(N/2 -1):(N/2)
        j = j+1 ;
        crosscorrelation = 0 ;
        for k = 1:N
            ktao = k-tao ;
            while (ktao<=0) 
                ktao = ktao + N ;
            end
            while (ktao > N)
                ktao = ktao - N ;
            end
            crosscorrelation = crosscorrelation + y(k)*u(ktao) ;
        end
        crosscorrelation_u(j) = crosscorrelation ;
    end
    numerator = fft( WHann .* crosscorrelation_u) ; 
    
    %compute smooth transfer function
    G_smooth_time = abs(numerator./denominator) ;
    G_smooth_time = reshape(G_smooth_time , [4,N/4]) ; 
    G_smooth_time = G_smooth_time(1,:)';

    %plot results
    subplot(2,2,i)
    loglog(G_smooth_time) ; hold on ;
    loglog(Gfreq) ;
    legend(int2str(gama(i)),'real')
    title('time domain window')
end 