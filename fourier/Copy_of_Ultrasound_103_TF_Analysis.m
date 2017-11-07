%%% Compare a pure frequency and  a time-frequency approach for maximum location
% IFPEN: 2014/09/15
% Spectral analysis of a series of ultrasound burst signals
%
% Uses FFTR.m, available on the webpage
% Uses stft2.M, from 
%
% WORK: extend to all signals, compare estimators for classification with
% respect to additive noise

% FLAG for: pre-window and pre-filter (1) or not (0)
f_PreFilter_PreWindow = 1;

% Change the path wrt the location of files
addpath('C:\Users\duvall\Documents\User_Personal\Laurent_Personnel_Confidentiel\Site_laurent-duval-lcd-siva\_WWW_Hidden\Documents-Supelec-FI');
load('Data_2013_Ultrasound_Spectral_Analysis_EA');
% Add the path to the toolbox and subdirs
addpath(genpath('D:\User\Laurent\Matlab\toolbox\DiscreteTFDs\'));

% Set basic signal parameters
time_Sampling = 0.0000010000;
freq_Sampling = 1/time_Sampling;

% Choose one specific signal and set time axis
index_Data_Demo = 150;
data_Demo = data_Mat(:,index_Data_Demo);
data_Length = length(data_Demo);
data_Demo_Time_Axis = time_Sampling*([1:data_Length]-1)';

% Definie apodizing window parameters
time_Burst_Max = 0.6e-3;
taper_Burst_Length = 20;
taper_Fractional_Order = 2;
sample_Burst_Max = max(find(data_Demo_Time_Axis <= time_Burst_Max))-1;

% Define low-cut, high-cut filter parameters (adapted to this signal, maybe not suitable for  other signals)
f_Low = 0.8e4;f_High = 20e4;
[B_cut,A_cut] = cheby2(2,40,[f_Low f_High]*2/freq_Sampling);


if f_PreFilter_PreWindow
  % Window the signal and fill with zeros
  data_Demo_Window = data_Demo.*...
    [Window_Raised_Frac_Sine(sample_Burst_Max,taper_Burst_Length,taper_Fractional_Order);zeros(data_Length-sample_Burst_Max,1)];
  
%    data_Demo_Window = data_Demo.*...
%     [hamming(sample_Burst_Max);zeros(data_Length-sample_Burst_Max,1)];
 
  % Define low-cut, high-cut filter parameters (adapted to this signal, maybe not to other signals)
   data_Demo_Filt = filtfilt(B_cut,A_cut,data_Demo_Window);
   data_Demo_Filt = data_Demo_Window;
else
  data_Demo_Filt = data_Demo;
end

% Set the time-frequency parameters
nfreq = 8*64;decf =1;
f_Subsampling = 1;
data_Demo_Stft = data_Demo_Filt(1:f_Subsampling:end);
% Compute the stft2 time-frequency representation
[tfd, t, f] = stft2(data_Demo_Stft, freq_Sampling/f_Subsampling, nfreq, decf);
f_red = f(nfreq/2-1:end);
abs_TFD = abs(tfd(1:nfreq/2+1,:));
[abs_FFT,abs_Axis] = FFTR(data_Demo_Stft,1/(freq_Sampling/f_Subsampling));

h1 = figure(1);clf
h=suptitle(['Signal n°',num2str(index_Data_Demo)]);set(h,'Interpreter','none');
subplot(3,1,1)
plot(t,data_Demo_Stft);axis tight
xlabel('Time');ylabel('Amp.')
subplot(3,1,2)
hold on
plot(abs_Axis,20*log10(abs_FFT))
xlabel('Freq. (Hz)');ylabel('Amp. (dB)')
[abs_FFT_Value,abs_FFT_Index] = max(abs_FFT);
plot(abs_Axis(abs_FFT_Index),20*log10(abs_FFT_Value),'or');
ylabel(['Freq. (Hz) [f_{max} = ',num2str(abs_Axis(abs_FFT_Index)),' Hz]'])

axis tight; grid on

subplot(3,1,3);
% figure(2);clf;
hold on
h= imagesc(t,f_red,(flipud(abs_TFD))); axis tight
xlabel('Time');ylabel('Freq. (Hz)')
% Find the maximum of the time-frequency plane
[abs_TFD_Max,abs_TFD_Idx] = max(abs_TFD(:)); [frequency_Idx,time_Idx] = ind2sub(size(abs_TFD),abs_TFD_Idx);
h=plot(t(time_Idx),f_red((end-frequency_Idx)),'or','MarkerSize',20) ;
set(gca,'YDir','reverse')
xlabel(['Time (s) [t_{max} = ',num2str(t(time_Idx)),' s]']);
ylabel(['Freq. (Hz) [f_{max} = ',num2str(f_red(frequency_Idx)/10),' Hz]'])

