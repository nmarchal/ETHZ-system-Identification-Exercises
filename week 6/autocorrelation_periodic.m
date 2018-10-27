function [ autocorrelation_y ] = autocorrelation_periodic( u, period, tao_arg )
%autocorrelation of periodic signals
% u is the signal for which we want the autocorrelation
% period is an integer
% tao is a vector of 2 elements, with the beginning and end for the
% autocorrelation

j = 0 ;
tao_min = min(tao_arg) ;
tao_max = max(tao_arg) ;

if(tao_min == tao_max)
    error('the interval for tao is incorrect \n')
end

for tao = tao_min:tao_max
    j = j+1 ;
    autocorrelation = 0 ;
    for k = 1:period
        ktao = k-tao ;
        while (ktao<=0) 
            ktao = ktao + period ;
        end
        while (ktao > period)
            ktao = ktao - period ;
        end
        autocorrelation = autocorrelation + u(k)*u(ktao) ;
    end
    autocorrelation_y(j) = autocorrelation ;
end
autocorrelation_y = autocorrelation_y./period ;
end

