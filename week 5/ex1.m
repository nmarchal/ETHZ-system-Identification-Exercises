clc
clear all
close all

signal = [] ;
order = [5 6 7 8] ;
period_length = 2.^order -1 ;
N = 1024 ;
Band = [0 1];
Range = [-1,1];
num_per = 4 ;

signal1 = fft(idinput([period_length(1),1,num_per],'prbs',Band,Range));
signal2 = fft(idinput([period_length(2),1,num_per],'prbs',Band,Range));
signal3 = fft(idinput([period_length(3),1,num_per],'prbs',Band,Range));
signal4 = fft(idinput([period_length(4),1,num_per],'prbs',Band,Range));

periodogram1 = 1/(num_per*period_length(1)).*(abs(signal1).^2) ;
periodogram1 = reshape(periodogram1 , [num_per,period_length(1)]) ;
periodogram1 = periodogram1(1,2:end)'; %first one is behaving akwardly

periodogram2 = 1/(num_per*period_length(1)).*(abs(signal2).^2) ;
periodogram2 = reshape(periodogram2 , [num_per,period_length(2)]) ;
periodogram2 = periodogram2(1,2:end)'; %first one is behaving akwardly

periodogram3 = 1/(num_per*period_length(1)).*(abs(signal3).^2) ;
periodogram3 = reshape(periodogram3 , [num_per,period_length(3)]) ;
periodogram3 = periodogram3(1,2:end)'; %first one is behaving akwardly

periodogram4 = 1/(num_per*period_length(1)).*(abs(signal4).^2) ;
periodogram4 = reshape(periodogram4 , [num_per,period_length(4)]) ;
periodogram4 = periodogram4(1,2:end)'; %first one is behaving akwardly

%MATLAB periodogram function not working

figure ; %my periodogram function
semilogx(periodogram1) ; hold on ;
semilogx(periodogram2) ; hold on ;
semilogx(periodogram3) ; hold on ;
semilogx(periodogram4) ; hold on ;

%FIXME : it is still differentfrom the correction
