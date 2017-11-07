function PChirplet_Example(),
% example to demonstrate the use of Polychirplet
% by Peng zhike
% email: z.peng@sjtu.edu.cn
% % June 1, 2008
%	Copyright (c) belongs to the authors of the paper 
%)	Z.K Peng, G. Meng, Z.Q. Lang, F.L Chu, W.M Zhang, Polynomial Chirplet Transform with Application to Instantaneous Frequency Estimation,
%   IEEE Transactions on Measurement and Instrumentation 60(2011) 3222-3229
%	Y Yang, Z.K Peng, G. Meng, W.M Zhang, A Novel Time Frequency Transform for the Analysis of Signals with Time-varying Frequency, 
%   IEEE Transactions on Industrial Electronics, 59(2012) 1612-1621
%   Y Yang, Z.K Peng, G. Meng, W.M Zhang, Characterize highly oscillating frequency modulation using generalized Warblet transform,
%   Mechanical Systems and Signal Processing, 26 (2012) 128-140
% 	Yang, Yang; Zhang, Wenming; Peng, Zhike; Meng, Guang, Multicomponent Signal Analysis Based on Polynomial Chirplet Transform, 
%   IEEE Transactions on Industrial Electronics, 60(2013), pp 3948-3956 
%   The citation about the papers must be included in all publications or
%   thesises as long as this program is used by anyone. 

close all
clear all
clc
format long

SampFreq = 200;
t = 0 : 1/SampFreq : 15;

c1 = 2 * pi * 10;            % initial frequency of the chirp excitation
c2 = 2 * pi * 5/2;           % set the speed of frequency change be 1 Hz/second
c3 = 2 * pi * 1/3;
c4 = 2 * pi * -1/40;

Sig = sin(c1 * t + c2 * t.^2 / 2 + c3 * t.^3 /3 + c4 * t.^4 /4);   % get the A(t)

n1 = rand(size(Sig));
n1 = n1 - mean(n1);
n1 = 3*n1 / max(n1) * max(Sig);
std(n1)
Sig = Sig + n1/2;

R1 = [0,0];
figure(1)
[Spec, f] = Polychirplet(Sig',SampFreq,R1,1024,512);

[v, I] = max(Spec,[],1);
[p, z] = polylsqr(t,f(I),4);
z
figure(2)
plot(t, f(I), t, p)

figure(3)
[Spec, f] = Polychirplet(Sig',SampFreq,z(2:end),1024,512);

[v, I] = max(Spec,[],1);
[p, z] = polylsqr(t,f(I),4);
z
figure(4)
plot(t, f(I), t, p)

figure(5)
[Spec, f] = Polychirplet(Sig',SampFreq,z(2:end),1024,512);

[v, I] = max(Spec,[],1);
[p, z] = polylsqr(t,f(I),4);
z
figure(6)
plot(t, f(I), t, p)

figure(7)
[Spec, f] = Polychirplet(Sig',SampFreq,z(2:end),1024,512);

[v, I] = max(Spec,[],1);
[p, z] = polylsqr(t,f(I),4);
z
figure(8)
plot(t, f(I), t, p)
