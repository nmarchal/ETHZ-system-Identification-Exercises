clc
clear all ;
close all ;

% second order parameter
zetaz = 0.1 ;
wz = 3.0 ;
zetap = 0.1 ;
wp = 3.5 ;
Ts = 0.02 ;

% construct the plant 
G1 = tf([1 2*zetaz*wz wz^2],[1 2*zetap*wp wp^2]) 
G2 = tf([5000],[1 250 200*50]) ;
Gs = G1*G2 ;
Gdz = c2d(Gs,Ts,'zoh') ; % discrete time

% Construc the controller 
Cdz = tf([1.25 -0.75],[1 -1],Ts) ;

%feedback loop transfer function
cl_tf = feedback (Gdz,Cdz) ; % not exactly what we have

%% 1

S = 1/(1+Cdz*Gdz) ;
% S = inv(1+Cdz*Gdz) ; 
Td = 1 - S ;

%calculate the poles
[Zs,Ps,Ks] = zpkdata(S) ;
Ps = cell2mat(Ps) ; 
[Zt,Pt,Kt] = zpkdata(Td) ;
Pt = cell2mat(Pt) ; %same poles

%check poles with unit circle
if isempty(find(abs(Ps)>1))
    fprintf('the poles are within the unit circle \n')
else
    fprintf('ERROR - The poles are not with unit circle \n')
end

%% 2

%create a prbs
Nper = 21 ; %by increasing this it is much better (around 1000 its really good)
Npt = 1023 ;
input_signal = idinput([Npt,1,Nper],'prbs',[0,1],[-0.1 0.1]) ;
%discard first period 
r = input_signal(Npt+1:end);
%noise
v = sqrt(0.1)*randn(Npt*(Nper-1),1) ;
% time
time = [0:1:(Npt*(Nper-1)-1)]'*Ts;
%frequencies
omega = exp(i*2*pi/Npt*[0:1:Npt-1]);
omega_plot = 2*pi/Npt*[0:1:Npt-1];
omega_plot_pos = omega_plot(1:floor(Npt/2)) ;
%true transfer functions
G_real = squeeze(freqresp(Gdz,omega));
G_real_plot = G_real(1:floor(Npt/2)) ;

%input output method
%y
y_tot = lsim(Gdz*Cdz/(1+Gdz*Cdz),r,time) + lsim(S,v,time);
y_tot = reshape(y_tot, Npt, (Nper-1)) ;
y = mean(y_tot,2) ;

%u
u_tot =  lsim(Cdz/(1+Gdz*Cdz),r,time);
u_tot = reshape(u_tot, Npt, (Nper-1)) ;
u = mean(u_tot,2) ;

%fft
U = fft(u) ;
Y = fft(y) ;

%G estimated
G_est_2 = Y./U ;
G_est_2_pos = G_est_2(1:floor(Npt/2)) ;

%bode plot
figure;
loglog(omega_plot_pos, abs(G_real_plot)) ;
grid on ; hold on ;
loglog(omega_plot_pos, abs(G_est_2_pos)) ;
legend('real', 'estimated')
title('plant')

%% part 3
r_tot = reshape(r, Npt, (Nper-1)) ;
error_tot = r_tot - (y_tot) ;
error = mean(error_tot,2) ;

% fft
E = fft(error) ;
R = fft(r(1:Npt)) ;

% estimate sensitivity
S_est = E./R ;
S_est_pos = S_est(1:floor(Npt/2)) ;

% true transfer function 
S_real = squeeze(freqresp(1/(1+Gdz*Cdz),omega));
S_real_pos = S_real(1:floor(Npt/2)) ;

%bode plot
figure;
loglog(omega_plot_pos, abs(S_real_pos)) ;
grid on ; hold on ;
loglog(omega_plot_pos, abs(S_est_pos)) ;
legend('real', 'estimated')
title('sensitivity')

%% part 4
w = r ;
W = R ;

% true tf
GS_real = squeeze(freqresp(Gdz/(1+Gdz*Cdz),omega));
GS_real_pos = GS_real(1:floor(Npt/2)) ;

% Calculate output
y_tot = lsim(Gdz/(1+Gdz*Cdz),w,time) + lsim(S,v,time);
y_tot = reshape(y_tot, Npt, (Nper-1)) ;
y = mean(y_tot,2) ;

% dft
Y = fft(y) ;

%estimation
GS_est = Y./W ;
GS_est_pos = GS_est(1:floor(Npt/2)) ;

%bode plot
figure;
loglog(omega_plot_pos, abs(GS_real_pos)) ;
grid on ; hold on ;
loglog(omega_plot_pos, abs(GS_est_pos)) ;
legend('real', 'estimated')
title('sensitivity')

%% 5) New estimate of G

G_new_est = S_est./GS_est ;
G_new_est_pos = G_new_est(1:floor(Npt/2)) ;

%bode plot
figure;
loglog(omega_plot_pos, abs(G_real_plot)) ;
grid on ; hold on ;
loglog(omega_plot_pos, abs(G_new_est_pos)) ;
loglog(omega_plot_pos, abs(G_est_2_pos)) ;
legend('real', 'indirect', 'ETFE')
title('plant (2nd method)')

%for m the ETFE is much better, however for the prof it is quite similar
