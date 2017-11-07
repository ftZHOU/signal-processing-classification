clear all;
close all;

addpath(genpath('/Users/zhoufengting/fourier/DiscreteTFDs'));

%retrieve music signal data and plot the signal in time domaine
[musicSignal,Fs]=audioread('fluteircam.wav');
%this whole music period is 13 seconds
Tm=length(musicSignal)/Fs;


%equivalent sample number for 40ms is 40ms*Fs(=32KHz)=1280
stationnarite=1280;
%pick a random order 
order=45;

%register the results for all the clips
freBurg=[];
dspBurg=[];

%select the frequency useful(ideally 20 notes finally) and their duration
noteBurg=[];
durationBurg=[];

%register last clip's frequency for reference
lastFreq=0;

%initialize the final music restoration 
music=[];

%calculate the running time as one factor of the performance
tic;

%calculate the mean of the whole music
px=periodogram(musicSignal);
pxMean=mean(px);



for i=1:341
    %resample the signal music by 40ms
    sampleNb = i;
    dataMusic=musicSignal(stationnarite*(sampleNb):stationnarite*(sampleNb+1));
    
    %calculate the DSP for each clip
    [aa, sigma2, ref, ff, mydsp] = myburg(dataMusic',order,Fs);
    pMean=mean(mydsp);
    %retrieve the maximum DSP in that clip
    [maxValue_DSP, indexMaxDSP]=max(mydsp);
    %retrieve the fondemental frenquency in this clip
    %try to add another condition of the peridogram
    %if(pMean>0.01*pxMean)
       frequency=abs(ff(indexMaxDSP));
       freBurg=[freBurg,frequency];
       dspBurg=[dspBurg,maxValue_DSP];
    
    %end

end

for ii=1:(length(freBurg)-1)
        if((abs(freBurg(ii+1)-lastFreq)>10))
             noteBurg=[noteBurg,lastFreq];
             durationBurg=[durationBurg,ii];
             lastFreq= freBurg(ii+1);
        else
            lastFreq=(lastFreq+freBurg(ii+1))/2;
        end   
end

% delete the frequency too high/low;
Locate=find(noteBurg<300);
noteBurg(Locate)=[];
durationBurg(Locate)=[];
Locate=find(noteBurg>1000);
noteBurg(Locate)=[];
durationBurg(Locate)=[];
% until this step the Burg can represent the music acceptable while with lots of
% noises(frequency not exists) very short, vert unstable 

% so I add another condition of the duration to detect if that is a noise
% peak
for jj=1:(length(noteBurg)-1)
    if((durationBurg(jj+1)-durationBurg(jj))<5)
       noteBurg(jj+1)=[1];
       
    end
end
 Location=find(noteBurg==1);
 noteBurg(Location)=[];
 durationBurg(Location)=[];

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

% %associate two lists of data useful and plot it
% C(2,:)=noteLeDu;
% C(1,:)=durationLeDu;
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
% for j= 1:length(durationBurg)
%     if(j==1)
%         t=0:1/Fs:durationBurg(j)*0.04;
%         music=[music,sin(2*pi*noteBurg(j)*t)];
%     else
%         t=durationBurg(j-1)*0.04:1/Fs:durationBurg(j)*0.04;
%         music=[music,sin(2*pi*noteBurg(j)*t)];
%     end
% 
% end
%   sound(music,Fs);
%  
%  
% %  record that audio 
%   filename = ('Burg.wav'); 
%   audiowrite(filename,music,Fs);

