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
order=250;

%register the results for all the clips
freLeDu=[];
dspLeDu=[];

%select the frequency useful(ideally 20 notes finally) and their duration
noteLeDu=[];
durationLeDu=[];

%register last clip's frequency for reference
lastFreq=0;

%initialize the final music restoration 
music=[];

%calculate the running time as one factor of the performance
tic;

for i=1:341
    %resample the signal music by 40ms
    sampleNb = i;
    dataMusic=musicSignal(stationnarite*(sampleNb):stationnarite*(sampleNb+1));
    
    %calculate the DSP for each clip
    [aa, sigma2, ref, ff, mydsp] = mylevinsondurbin(dataMusic',order,Fs);
    %retrieve the maximum DSP in that clip
    [maxValue_DSP, indexMaxDSP]=max(mydsp);
    %retrieve the fondemental frenquency in this clip
    frequency=abs(ff(indexMaxDSP));
    freLeDu=[freLeDu,frequency];
    dspLeDu=[dspLeDu,maxValue_DSP];

end


% 
% for ii=1:340
%     if(freLeDu(ii)>300 & freLeDu(ii)<900 & abs(freLeDu(ii+1)-freLeDu(ii))>15)
%         noteLeDu=[noteLeDu,lastFreq];
%         durationLeDu=[durationLeDu,ii];
%     elseif(freLeDu(ii)>300 & freLeDu(ii)<900 & abs(freLeDu(ii+1)-freLeDu(ii))<=15)
%         lastFreq=(lastFreq+freLeDu(ii))/2;
%     end
% end

for ii=1:(length(freLeDu)-1)
        if((abs(freLeDu(ii+1)-lastFreq)>10))
             noteLeDu=[noteLeDu,lastFreq];
             durationLeDu=[durationLeDu,ii];
             lastFreq= freLeDu(ii+1);
        else
            lastFreq=(lastFreq+freLeDu(ii+1))/2;
        end   
end

% delete the frequency too high/low;
Locate=find(noteLeDu<300);
noteLeDu(Locate)=[];
durationLeDu(Locate)=[];
Locate=find(noteLeDu>1000);
noteLeDu(Locate)=[];
durationLeDu(Locate)=[];
% 
for jj=1:(length(noteLeDu)-1)
    if((durationLeDu(jj+1)-durationLeDu(jj))<5)
       noteLeDu(jj+1)=[1];
       
    end
end
 Location=find(noteLeDu==1);
 noteLeDu(Location)=[];
 durationLeDu(Location)=[];
%  
% tt=toc;
% all the part below is the figure or sound presentation 

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
% plot(t1,dspLeDu,'p');
% hold on;
% title('DSP of LevinsonDurbin');
% xlabel('time/ s');
% ylabel('DSP／ W/Hz');

% %plot the frequency in order of time
% t1=[(Tm/341):(Tm/341):Tm];
% figure(3)
% plot(t1,freLeDu,'p');
% hold on;
% title('DSP of LevinsonDurbin');
% xlabel('time/ s');
% ylabel('frequency／Hz');

% %associate two lists of data useful and plot it
C(2,:)=noteLeDu;
C(1,:)=durationLeDu;
figure(3)
stairs(C(1,:)*0.04,C(2,:),'LineWidth',2);
hold on;
grid on;
xlabel('time/s');
ylabel('frequency/Hz');

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

%try to replay the music(without adding the strength factor)
for j= 1:length(noteLeDu)
    if(j==1)
        t=0:1/Fs:durationLeDu(j)*0.04;
        music=[music,sin(2*pi*noteLeDu(j)*t)];
    else
        t=durationLeDu(j-1)*0.04:1/Fs:durationLeDu(j)*0.04;
        music=[music,sin(2*pi*noteLeDu(j)*t)];
    end

end
 sound(music,Fs);


%  record that audio 
%  filename = ('LeDu.wav'); 
%  audiowrite(filename,music,Fs);


