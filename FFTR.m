function [fftR,fftAxe] = FFTR(data,timeSamp,fftLength,fftDim);

% Computes and displays (default) the one-sided FFT for REAL 'data'
% [fftR,fftAxe] = FFTR(data,timeSamp,fftLength,fftDim)
% plus the frequency axis, with harmonic amplitude preservation
% 'timeSamp' = 1 (default)
% 'fftLength' equals data length (default) or specified
% length or the next power of two (fftLength = 'n')
% 'fftDim' applies the FFT operation across the
%  dimension DIM.
%	
% Example:
%   length_Signal = 512;
%   time = (0:length_Signal-1)'/length_Signal;
%   time_Sampling = median(diff(time));
%   data = cos(64*pi*time);
% 
%   figure(1);clf
%   subplot(2,1,1)
%   plot(time,data);axis tight;grid on;xlabel('Time (s)');ylabel('Amplitude (a. u.)');
%   subplot(2,1,2);hold on
%   [fft_R,fft_Axe] = FFTR(data,time_Sampling);
%   [fft_R_Max,fft_R_Idx] = max(fft_R);
%   h1 = plot(fft_Axe,fft_R,'x');axis tight;grid on;xlabel('Frequency (Hz)');ylabel('Amplitude (a. u.)');
%   h2 = plot(fft_Axe(fft_R_Idx),fft_R(fft_R_Idx),'ro');
%   set([h1,h2],'LineWidth',6)

% Uses: fft
%
%	Author: Laurent C. Duval
%	Institution: IFP Energies nouvelles, Technology Department
%	Website: http://www.laurent-duval.eu
%	Created: 05/07/2002
%	Modified: 21/08/2016

% Avoid complex data inputs
if ~isreal(data) 
   error('Input data is complex')
end
if nargin < 4
   fftDim = 1;
end

if nargin < 2
   timeSamp = 1;
end
freqSamp = 1/timeSamp;
freqNyq = freqSamp/2;
[dataLength,dataNb] = size(data);

% Choose FFT length
% If not given, simply dataLength, else 'n' for the next
% power of two
if (nargin < 3) 
   fftLength = dataLength;
else
   if ~isnumeric(fftLength)
      switch lower(fftLength)
      case 'n'
         fftLength = 2^nextpow2(dataLength);
      otherwise
         error('Ambiguous ''fftLength'' argument!')
      end
   end
end
nbFreq = ceil((fftLength+1)/2);
dataFft = fft(data,fftLength,fftDim);
fftR = 2*abs(dataFft(1:nbFreq,:));
% Compensate for DC 2-factor
fftR(1,:) = fftR(1,:)/2;
% Compensate for last frequency 2-factor 
if ~rem(fftLength,2),
   fftR(nbFreq,:) = fftR(nbFreq,:)/2;
end
% Data length independency in spectrum peak amplitude
fftR = fftR/dataLength;
fftAxe = (0:nbFreq-1)'*2*freqNyq/fftLength;
if nargout == 0
      semilogx(fftAxe,20*log10(fftR));axis tight;grid on
      xlabel('Freq. (Hz)')
      ylabel('Magnitude (dB)')
end


