clear all;
close all;

addpath(genpath('/Users/zhoufengting/projetduval/Data/DiscreteTFDs'));

%retrieve music signal data and plot the signal in time domaine
[musicSignal,Fs]=audioread('fluteircam.wav');
%this whole music period is 13 seconds
Tm=length(musicSignal)/Fs;


%equivalent sample number for 40ms is 40ms*Fs(=32KHz)=1280
stationnarite=1280;
%pick a random order 
order=100;

%register the results for all the clips
frePisa=[];
dspPisa=[];

%select the frequency useful(ideally 20 notes finally) and their duration
notePisa=[];
durationPisa=[];

%register last clip's frequency for reference
lastFreq=0;

%initialize the final music restoration 
music=[];

%calculate the running time as one factor of the performance
tic;

% set the recursion
rec=2;



for i=1:341
    %resample the signal music by 40ms
    sampleNb = i;
    dataMusic=musicSignal(stationnarite*(sampleNb):stationnarite*(sampleNb+1));
    
    %calculate the DSP for each clip
    [ff, mydsp] = mypisarenko(dataMusic', order, Fs, rec);
    %retrieve the maximum DSP in that clip
    [maxValue_DSP, indexMaxDSP]=max(mydsp);
    %retrieve the fondemental frenquency in this clip
    frequency=abs(ff(indexMaxDSP));
    frePisa=[frePisa,frequency];
    dspPisa=[dspPisa,maxValue_DSP];

end

for ii=1:(length(frePisa)-1)
        if((abs(frePisa(ii+1)-lastFreq)>10))
             notePisa=[notePisa,lastFreq];
             durationPisa=[durationPisa,ii];
             lastFreq= frePisa(ii+1);
        else
            lastFreq=(lastFreq+frePisa(ii+1))/2;
        end   
end

% delete the frequency too high/low;
Locate=find(notePisa<300);
notePisa(Locate)=[];
durationPisa(Locate)=[];
Locate=find(notePisa>1000);
notePisa(Locate)=[];
durationPisa(Locate)=[];
% until this step the Burg can represent the music acceptable while with lots of
% noises(frequency not exists) very short, vert unstable 

% so I add another condition of the duration to detect if that is a noise
% peak
for jj=1:(length(notePisa)-1)
    if((durationPisa(jj+1)-durationPisa(jj))<5)
       notePisa(jj+1)=[1];
       
    end
end
 Location=find(notePisa==1);
 notePisa(Location)=[];
 durationPisa(Location)=[];

tt=toc;



%all the part below is the figure or sound presentation 

%plot the signal in the time domaine
% t=[(1/Fs):(1/Fs):Tm];
% t1=[(Tm/341):(Tm/341):Tm];
% figure(1)
% plot(t,musicSignal);
% hold on;
% title('Time Domaine');
% xlabel('time/s');
% ylabel('magnitude');
% hold off;

% %plot the dsp in order of time
% t1=[(Tm/341):(Tm/341):Tm];
% figure(2)
% plot(t1,dspBurg,'p');
% hold on;
% title('DSP of Burg');
% xlabel('time/ s');
% ylabel('DSP／ W/Hz');

% %plot the frequency in order of time
% t1=[(Tm/341):(Tm/341):Tm];
% figure(3)
% plot(t1,freBurg,'p');
% hold on;
% title('DSP of Burg');
% xlabel('time/ s');
% ylabel('frequency／Hz');

% % %associate two lists of data useful and plot it
% C(2,:)=notePisa;
% C(1,:)=durationPisa;
% figure(3)
% stairs(C(1,:)*0.04,C(2,:),'LineWidth',2);
% hold on;
% grid on;
% xlabel('time/s');
% ylabel('frequency/Hz');

% %plot one window's enveloppe 
% A(1,:)=ff;
% A(2,:)=mydsp;
% B=sortrows(A',1);
% D=B';
% env=envelope(D(2,:),300,'peaks');
% env(env<0)=0;
% env(env<0.05)=0;
% figure(4)
% subplot 211
% plot(D(1,:),D(2,:),'linewidth',2);
% subplot 212
% plot(D(2,:),'linewidth',2);
% hold on;
% plot(env,'linewidth',2);
% hold off;

% try to replay the music(without adding the strength factor)
for j= 1:length(durationPisa)
    if(j==1)
        t=0:1/Fs:durationPisa(j)*0.04;
        music=[music,sin(2*pi*notePisa(j)*t)];
    else
        t=durationPisa(j-1)*0.04:1/Fs:durationPisa(j)*0.04;
        music=[music,sin(2*pi*notePisa(j)*t)];
    end

end
 sound(music,Fs);
 
 

