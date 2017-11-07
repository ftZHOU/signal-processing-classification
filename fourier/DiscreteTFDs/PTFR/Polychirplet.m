function [Spec,f] = Polychirplet(Sig,SampFreq,Ratio,N,WinLen);

% Polynomial chirplet transform
%
%	   Sig       : the signal to be analyzed
%    SampFreq    : sampling frequency 
%      Ratio     : Coefficients of polynomial
%        N       : the number of frequency bins (default : length(Sig)).
%	WinLen       : the length of window used to locate the signal in time.
%
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
%   thesises as long as this program is used by anyone. %

if(nargin < 3),
    error('At least 3 inputs are required!');
end

SigLen = length(Sig);

if (nargin < 5),
    WinLen = SigLen / 4;
end

if (nargin < 4),
    N = SigLen;
end

if(N > 512),
    N = 512;
end

RatioNum = length(Ratio);

dt = (0:(SigLen-1))';
dt = dt / SampFreq;

df = zeros(size(dt));

for k = 1:RatioNum,
    df = df + Ratio(k) * dt.^k;
end

kernel = zeros(size(dt));

for k = 1:RatioNum,
    kernel = kernel + Ratio(k)/(k + 1) * dt.^(k + 1);
end

WinLen = ceil(WinLen / 2) * 2;
t = linspace(-1,1,WinLen)';
WinFun = exp(log(0.005) * t.^2 );
WinFun = WinFun / norm(WinFun);

Lh = (WinLen - 1)/2; 

Spec = zeros(N,SigLen) ;  
Rt = zeros(N,1);
Rdt = zeros(N,1);

wait = waitbar(0,'Please wait...');
for iLoop = 1:SigLen,

    waitbar(iLoop/SigLen,wait);
    
    tau = -min([round(N/2)-1,Lh,iLoop-1]):min([round(N/2)-1,Lh,SigLen-iLoop]);
    temp = floor(iLoop + tau);
    Rt(1:length(temp)) = kernel(temp); 
    Rdt(1:length(temp)) = dt(temp);  
    
    temp1 = floor(Lh+1+tau);
    rSig = Sig(temp);
    rSig = Hilbert(real(rSig));

    rSig = rSig .* conj(WinFun(temp1));
    Spec(1:length(rSig),iLoop) = rSig;
    
    Spec(:,iLoop) = Spec(:,iLoop) .* exp(-j * 2.0 * pi * (Rt - df(iLoop) * Rdt)); 

%    ft = abs(fft(Spec(:,iLoop)));
%    Spec(:,iLoop) = ft;
    
end;

Spec = fft(Spec); 
Spec = abs(Spec);

close(wait);

Spec = Spec(1:round(end/2),:);
[nLevel, SigLen] = size(Spec);

f = [0:nLevel-1]/nLevel * SampFreq/2;
t = (0: SigLen-1)/SampFreq;

[fmax fmin] = FreqRange(Sig);
fmax = fmax * SampFreq;
fmin = fmin * SampFreq;

clf
%=====================================================%
% Plot the result                                     %
%=====================================================%
set(gcf,'Position',[20 100 350 300]);	    
set(gcf,'Color','w');					    

mesh(t,f,Spec);  
axis([min(t) max(t) fmin fmax]);
ylabel('Freq / Hz');
xlabel('Time / Sec')
Info = 'C = ';

for i = 1:RatioNum,
    Info = [Info,num2str(Ratio(i),4), '  '];
end

if RatioNum == 0,
    Info = 'C = 0';
end

%title(['Nonlinear Chirplet[',Info,']']);
