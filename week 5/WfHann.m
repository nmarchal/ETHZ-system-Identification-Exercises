function [omega,WHann] = WfHann(gamma,N)
%------------------------------------------------------------------
%
%   [omega,WHann] = WfHann(gamma,N)
%
%   Create a frequency domain Hann window with width parameter gamma
%   and data length N.  The Hann window is a raised cosine.
%
%   Roy Smith,  18 October, 2017.
%
%                6 November, 2017.  Fixed bug in N even indexing.
%
%------------------------------------------------------------------

if nargin == 0,
    disp('Syntax: [omega,W] = WfHann(gamma,N)')
    return
elseif nargin ~= 2,
    error('incorrect number of input arguments (2 expected)')
    return
end

%   basic parameter checking
if length(gamma) > 1,
    error('Width parameter, gamma, must be a scalar');
end
if round(gamma) ~= gamma,
    error('Width parameter, gamma, must be an integer');
end
if gamma < 1,
    error('Width parameter, gamma, must be positive');
end
if length(N) > 1,
    error('Calculation length, N, must be a scalar');
end
if round(N) ~= N,
    error('Calculation length, N, must be an integer');
end
if N < 1,
    error('Calculation length, N, must be positive');
end

%   The simplest approach is to define the window in the time domain and
%   then transform it to the frequency domain.

lags = [floor(-N/2+1):floor(N/2)]';
wHann = 0*lags;
idx = find(abs(lags) <= gamma);
wHann(idx) = 0.5*(1+cos(pi*lags(idx)/gamma));

%   
zidx = find(lags==0);    % index of the zero point.

wH_raw = fft([wHann(zidx:N);wHann(1:zidx-1)]);
WHann(zidx:N) = wH_raw(1:N-zidx+1);  % shift +ve freq to end
WHann(1:zidx-1) = wH_raw(N-zidx+2:N);% shift > pi freq to beginning
WHann = real(WHann);   % should actually be real
omega = 2*pi/N*lags;

return
  
%------------------------------------------------------------------
