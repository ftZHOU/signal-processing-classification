function [C,Ceq,DC,DCeq] = chirpcon(theta, x)
% g_chirp -- gradient vector for constr
% 
%  Usage
%    g = g_chirp(theta, x)
%
%  Inputs
%    theta    initial estimate ([time freq chirp_rate duration])
%    x    signal being fitted
%
%  Outxuts
%    g    value of the gradient

% Copyright (C) -- see DiscreteTFDs/Copyright

% f(t,f,c,d) = | z(t,f,c,d)|^2
% d/dt f(t,f,c,d) = 2 re{ (d/dt z) z*}


C = [];
Ceq = [];
DC = [];
DCeq = [];