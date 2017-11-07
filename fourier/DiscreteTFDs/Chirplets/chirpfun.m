function [v, g] = chirpfun(theta, x)
% f_chirp -- objective function minimized by constr
% 
%  Usage
%    v = f_chirp(theta, x)
%
%  Inputs
%    theta  initial estimate ([time freq chirp_rate duration])
%    x      signal being fitted
%
%  Outputs
%    v    value of the objective function

% Copyright (C) -- see DiscreteTFDs/Copyright

N = length(x);
v = -abs(x' * chirplets(N, [1 theta]))^2;


N = length(x);

t = theta(1);
f = theta(2);
c = theta(3);
d = theta(4);
n = (1:N)';
y = chirplets(N,[1 theta]);
z_conj = conj(sum(x.*conj(y)));

dz_dt = -sum( x .* conj(y) .* ( (n-t)/2/d^2+1i*c*(n-t)+1i*f ) );
g(1) = 2*real(dz_dt * z_conj);

dz_df = -sum( x .* conj(y) .* (-1i*(n-t)) );
g(2) = 2*real(dz_df * z_conj);

dz_dc = -sum( x .* conj(y) .* (-1i/2*(n-t).^2) );
g(3) = 2*real(dz_dc * z_conj);

dz_dd = -sum( x .* conj(y) .* ( (n-t).^2/2/d^3 - ...
                   1/2/d) );
g(4) = 2*real(dz_dd * z_conj);