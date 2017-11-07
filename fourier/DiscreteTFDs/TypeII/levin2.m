function [tfd, t, f] = levin2(x, fs, nfreq)
% levin2 -- Compute samples of the type II Levin distribution.
%
%  Usage
%    [tfd, t, f] = levin2(x, fs, nfreq)
%
%  Inputs
%    x     signal vector
%    fs    sampling frequency of x (optional, default is 1 sample/second)
%    nfreq number of samples to compute in frequency (optional, default
%          is twice the length of x)
%
%  Outputs
%    tfd  matrix containing the binomial distribution of signal x.  If x has
%         length N, then tfd will be nfreq by N. (optional)
%    t    vector of sampling times (optional)
%    f    vector of frequency values (optional)
%
% If no output arguments are specified, then the Levin distribution is 
% displayed using ptfd(tfd, t, f).

% Copyright (C) -- see DiscreteTFDs/Copyright

% specify defaults
x = x(:);
N = length(x);

error(nargchk(1, 3, nargin));
if (nargin < 3)
  nfreq = 2*N;
end
if (nargin < 2)
  fs = 1;
end

% make the acf for positive tau
gacf = zeros(N);
for t=1:N;
  gacf(1:t,t) = conj(x(t)) * x(t:-1:1);
end

% negative tau
gacf = [gacf ; zeros(nfreq-2*N+1,N) ; conj(flipud(gacf(2:N,:)))];

%compute the tfd
tfd = real(fft(gacf));
tfd = tfdshift(tfd)/nfreq;

t = 1/fs * (0:N-1);
f = -fs/2:fs/nfreq:fs/2;
f = f(1:nfreq);

if (nargout == 0)
  ptfd(tfd, t, f);
  clear tfd
end
