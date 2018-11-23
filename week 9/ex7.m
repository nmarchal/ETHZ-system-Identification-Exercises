clc
clear all
close all

N = 128 ;
a = 0.5 ;
b = 1 ;
m = 1 ;
n = 1;
batch_size = 1000 ;

NORMAL = true % distribution of w (uniform if false)

var = 0.2 ;
param_uniform = 3/5 ; %var uniform = (b-a)^2/12


%% 1 - a) b) w normal and uniform distribution
iter = 0 ;
theta = [] ;
try_size = [2^6 2^8 2^10  2^12 ] ; %
error_rss = zeros(batch_size, max(size(try_size))) ;

for N = try_size 
    clear y ; clear w ;
    iter = iter + 1 
    u = randn (N,1) ;
    for counter = 1 :batch_size 
        counter ;
        if NORMAL
            w = sqrt(var)* randn (N,1) ; 
        else
            w = param_uniform - 2*param_uniform*rand(N,1) ;
        end

        y = [] ;
        y = w(1) ;
        for i = 2:N
            y = [y (a*y(i-1)+b*u(i-1)+w(i))] ;
        end
        y = y' ;
        
        % construct phi
        for k = 0:(N-1) 
            for i = 1:(n+m)
                if i > n
                    j = i-n ;
                else
                    j = i ;
                end

                if (k-j)<0
                    phi(k+1,j) = 0 ; % negative index for y or u
                elseif i > n ;% add u
                    phi(k+1,i) = u(k+1-j) ;
                else % add y
                    phi(k+1,i) = - y(k+1-j) ;
                end
            end
        end
        %least square error
        theta(:,counter,iter) = phi\y ;
        y_predicted = phi*theta(:,counter,iter) ;
        %rss error
%         size(phi*theta) 
        error_rss(counter,iter) = sum((y - y_predicted).^2)/var ;
    end
end

%% plot histogram
close all ;

for i = 1:max(size(try_size))
    figure ;
    subplot(2,1,1)
    hist(-theta(1,:,i),linspace(0.3, 0.7, 1000)) ;
    subplot(2,1,2)
    hist(theta(2,:,i),linspace(0.8, 1.2, 1000))
end


%% 2) RSS

% plot histogram
close all ;
limit_bin = ([30 110 ; 180 320 ; 850 1150 ; 3800 4500]);

for i = 1:max(size(try_size))
    figure(i) ;
    subplot(3,1,1)
%     hist(error_rss(:,i),linspace(limit_bin(i,1),limit_bin(i,2),10))
    hist(error_rss(:,i),10)
    title('histogram')
end

%normal // FOR SOME REASON NOT REALLY WHAT i EXPECT
for i = 1:max(size(try_size))
    dist = (try_size(i)-2)+ sqrt(2*(N-2))*randn(try_size(i),1) ;
    size(dist)
    figure(i) ;
    subplot(3,1,2)
%     hist(dist,linspace(limit_bin(i,1),limit_bin(i,2),10)) ;
    hist(dist,10) ;
    title('normal distribution')
end

%chi squared
for i = 1:max(size(try_size))
    dist = chi2rnd((try_size(i)-2), try_size(i),1) ;
    figure(i) ;
    subplot(3,1,3)
%     hist(dist,linspace(limit_bin(i,1),limit_bin(i,2),10)) ;
    hist(dist,10) ;
    title('chi^2 distribution')
end