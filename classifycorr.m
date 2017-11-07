% start from the demo existed, plot 1448 signals with (time-frequancyMax)
% in one plan to classify 

% set the working path
addpath('/Users/zhoufengting/projetduval/Data');
load('4096samples_1448waveforms.mat');
% Add the path to the toolbox and subdirs
addpath(genpath('/Users/zhoufengting/fourier/DiscreteTFDs'));

% Set basic signal parameters
time_Sampling = 0.0000010000;
freq_Sampling = 1/time_Sampling;

% creat one matrix to register all the 
freq_time_matrix=[];
freq_matrix=[];
time_matrix=[];

% choose all signals and set time axis
for index_Data_Demo=1:1448
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
    
    % Window the signal and fill with zeros
    % data_Demo_Window = data_Demo.*...
    %[Window_Raised_Frac_Sine(sample_Burst_Max,taper_Burst_Length,taper_Fractional_Order);zeros(data_Length-sample_Burst_Max,1)];
  
    % data_Demo_Window = data_Demo.*...
    %[hamming(sample_Burst_Max);zeros(data_Length-sample_Burst_Max,1)];
    
    % data_Demo_Window = data_Demo.*...
    % [blackman(sample_Burst_Max);zeros(data_Length-sample_Burst_Max,1)];
 
    % Define low-cut, high-cut filter parameters (adapted to this signal, maybe not to other signals)
    %data_Demo_Filt = filtfilt(B_cut,A_cut,data_Demo_Window);
    data_Demo_Filt = data_Demo;
    
    % Set the time-frequency parameters
    nfreq = 8*64;decf =1;
    f_Subsampling = 1;
    data_Demo_Stft = data_Demo_Filt(1:f_Subsampling:end);
    % Compute the stft2 time-frequency representation
    [tfd, t, f] = stft2(data_Demo_Stft, freq_Sampling/f_Subsampling, nfreq, decf);
    f_red = f(nfreq/2-1:end);
    abs_TFD = abs(tfd(1:nfreq/2+1,:));
    [abs_FFT,abs_Axis] = FFTR(data_Demo_Stft,1/(freq_Sampling/f_Subsampling));
    
    % find the maximum frequency and its relevent time 
    [abs_FFT_Value,abs_FFT_Index] = max(abs_FFT);
    time_matrix=[time_matrix,abs_FFT_Index];
    freq_matrix=[freq_matrix,abs_FFT_Value];
    
end

freq_time_matrix(1,:)=time_matrix;
freq_time_matrix(2,:)=freq_matrix;

% figure
% plot(freq_time_matrix(1,:),freq_time_matrix(2,:),'p')