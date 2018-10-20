function [lags,wHann] = WtHann(gamma,N)
%------------------------------------------------------------------
%
%   [lags,wHann] = WtHann(gamma,N)
%
%   Create a Hann window with width parameter gamma and data length N.
%   The Hann window is a raised cosine.
%
%   Roy Smith,  18 October, 2017.
%
%------------------------------------------------------------------

if nargin == 0,
    disp('Syntax: [lags,w] = WtHann(gamma,N)')
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

lags = [floor(-N/2+1):floor(N/2)]';
wHann = 0*lags;
idx = find(abs(lags) <= gamma);
wHann(idx) = 0.5*(1+cos(pi*lags(idx)/gamma));

return
  
%------------------------------------------------------------------
