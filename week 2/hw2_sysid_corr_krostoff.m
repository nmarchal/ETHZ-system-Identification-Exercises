clear all
close all

% Vector array length
N=1024;

%% Matlab Task 2, Section 1

e = randn(N,1,'double');

%% Matlab Task 2, Section 2

e_p = fft(e);
for i=1:N
    per(i) = 1/N*e_p(i)*conj(e_p(i));
end


%% Matlab Task 2, Section 3

num=[1];
den=[1 -0.9 0.5];

P = tf([1], [1 -0.9 0.5], 1);

sample = [0:1:N-1];

w = lsim(P,e',sample);

vec = linspace(0,2*pi,N);

for i=1:N
    P_n(i)=abs(evalfr(P,vec(i)));
end

P_sq = P_n.*conj(P_n);


%% Matlab Task 2, Section 4

w_p = fft(w);

for i=1:N
    per_w(i) = 1/N*w_p(i)*conj(w_p(i));
end

loglog(per_w); hold on;
loglog(P_sq); hold on;
loglog(abs(per_w-P_sq));

%% Matlab Task 2, Section 5


