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

% du coup pour une tranche de 40ms, ça fait 40e-3/time_Sampling = 400
% samples
stationnarite=400;

% creat one matrix to register all the 
freq_time_matrix=[];
freq_matrix=[];
time_matrix=[];

% set order
order=1;

% choose all signals and set time axis
% for the first time deal with just first column 
for index_Data_Demo=1:1
    data_Demo = data_Mat(:,index_Data_Demo);
    data_Length = length(data_Demo);
    Time=data_Length*time_Sampling;
    
    % chaque colonne de signaux ont que 40ms du coup il y a que 10 tranches
    for ii=1:9
    dataFuite=data_Demo(stationnarite*ii:stationnarite*(ii+1));
    % Levinson durbin
    [aa, sigma2, ref, ff, mydsp] = mylevinsondurbin(dataFuite',order,freq_Sampling);
    %retrieve the maximum DSP in that signal column 
    %[maxValue_DSP, indexMaxDSP]=max(mydsp);
%     [maxValue_ff, indexMaxff]=max(ff);
%     frequency=abs(ff(indexMaxDSP));
%     freq_matrix=[freq_matrix,frequency];
%     time_matrix=[time_matrix,indexMaxff];

    end
    
   
    
end

freq_time_matrix(1,:)=time_matrix;
freq_time_matrix(2,:)=freq_matrix;


