%clear all;
%close all;
clc;

% FLAG for: pre-window and pre-filter (1) or not (0)
f_PreFilter_PreWindow = 1;

% Change the path wrt the location of files
addpath('\\ntelev.metz.supelec.centralesupelec.local\Users\masclet_hug\fourier');
load('Supelec_2012_SIR_Spectral_Analysis_EA_v001.mat');
% Add the path to the toolbox and subdirs
addpath(genpath('\\ntelev.metz.supelec.centralesupelec.local\Users\masclet_hug\fourier\DiscreteTFDs'));

% Set basic signal parameters
time_Sampling = 0.0000010000; %Sampling 1 MHz
freq_Sampling = 1/time_Sampling;

%-----------------------------------------------------------------

n = size(data_Mat,2)

analyse = zeros(4,n);

for s = 1:size(analyse,2) 
    
    % Choose one specific signal and set time axis
    index_Data_Demo = s;
    data_Demo = data_Mat(:,index_Data_Demo);
    data_Length = length(data_Demo);
    data_Demo_Time_Axis = time_Sampling*([1:data_Length])';

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
    [abs_FFT_Value,abs_FFT_Index] = max(abs_FFT);

    fq_max = abs_Axis(abs_FFT_Index);
    amp_max = 20*log10(abs_FFT_Value);
    % Find the maximum of the time-frequency plane
    [abs_TFD_Max,abs_TFD_Idx] = max(abs_TFD(:));
    [frequency_Idx,time_Idx] = ind2sub(size(abs_TFD),abs_TFD_Idx);

    t_max = t(time_Idx);
    f_max = f_red(frequency_Idx)/10;
    freq_bar = freqBar(abs_FFT, abs_Axis);
    
    analyse(1,s) = t_max;
    analyse(2,s) = fq_max; %Frequence maximale 
    analyse(3,s) = freq_bar; %Frequence barycentrique
    analyse(4,s) = amp_max; %Amplitude 
end

analyse
