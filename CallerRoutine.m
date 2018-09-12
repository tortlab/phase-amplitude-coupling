%% Constructing an example

clear
clc

data_length = 2^15;
srate=1024;
dt = 1/srate;
t=dt*(1:data_length);

nonmodulatedamplitude=2; % increase this to get less modulation (lower MI value)

Phase_Modulating_Freq=10;
Amp_Modulated_Freq=80;

lfp=(0.2*(sin(2*pi*t*Phase_Modulating_Freq)+1)+nonmodulatedamplitude*0.1).*sin(2*pi*t*Amp_Modulated_Freq)+sin(2*pi*t*Phase_Modulating_Freq);
lfp=lfp+1*randn(1,length(lfp));

%% Working with actual LFP signals (example .mat file available at github)

% % load('LFP_HG_HFO.mat')
% % lfp = lfpHG; % or lfp = lfpHFO;
% % data_length = length(lfp);
% % srate = 1000;
% % dt = 1/srate;
% % t = (1:data_length)*dt;


%% Plotting the signal

clf 
subplot(2,1,1)
plot(t,lfp)
xlim([0 1])
set(gca,'fontsize',14)
xlabel('time (s)')
ylabel('mV')

%% Define the amplitude- and phase-frequencies

PhaseFreqVector=2:2:50;
AmpFreqVector=10:5:200;

PhaseFreq_BandWidth=4;
AmpFreq_BandWidth=20;


%% Define phase bins

nbin = 18; % number of phase bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    position(j) = -pi+(j-1)*winsize; 
end

%% Filtering and Hilbert transform

'CPU filtering'
tic
Comodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);

for ii=1:length(AmpFreqVector)
    Af1 = AmpFreqVector(ii);
    Af2=Af1+AmpFreq_BandWidth;
    AmpFreq=eegfilt(lfp,srate,Af1,Af2); % filtering
    AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj);
    Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % filtering 
    PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % getting the phase time series
end
toc

%% Compute MI and comodulogram

'Comodulation loop'

counter1=0;
for ii=1:length(PhaseFreqVector)
counter1=counter1+1;

    Pf1 = PhaseFreqVector(ii);
    Pf2 = Pf1+PhaseFreq_BandWidth;
    
    counter2=0;
    for jj=1:length(AmpFreqVector)
    counter2=counter2+1;
    
        Af1 = AmpFreqVector(jj);
        Af2 = Af1+AmpFreq_BandWidth;
        [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
        Comodulogram(counter1,counter2)=MI;
    end
end
toc

%% Plot comodulogram

clf
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulogram',30,'lines','none')
set(gca,'fontsize',14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Use the routine below to look at specific pairs of frequency ranges:

Pf1 = 6;
Pf2 = 12;
Af1 = 60;
Af2 = 100;

[MI,MeanAmp] = ModIndex_v1(lfp,srate,Pf1,Pf2,Af1,Af2,position);

bar(10:20:720,[MeanAmp,MeanAmp]/sum(MeanAmp),'k')
xlim([0 720])
set(gca,'xtick',0:360:720)
xlabel('Phase (Deg)')
ylabel('Amplitude')
title(['MI = ' num2str(MI)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Or use the routine below to make a comodulogram using ModIndex_v1; this takes longer than
%% the method outlined above using ModIndex_v2 because in this routine multiple filtering of the same
%% frequency range is employed (the Amp frequencies are filtered multiple times, one
%% for each phase frequency). This routine might be the only choice though
%% for computers with low memory, because it does not create the matrices
%% AmpFreqTransformed and PhaseFreqTransformed as the routine above

tic

% define phase bins
nbin = 18; % number of phase bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    position(j) = -pi+(j-1)*winsize; 
end

% define analyzed frequency ranges
PhaseFreqVector=2:2:50;
AmpFreqVector=10:5:200;
PhaseFreq_BandWidth=4;
AmpFreq_BandWidth=20;


Comodulogram=zeros(length(PhaseFreqVector),length(AmpFreqVector));

counter1=0;
for Pf1=PhaseFreqVector
    counter1=counter1+1;
    Pf1 % just to check the progress
    Pf2=Pf1+PhaseFreq_BandWidth;
    
    counter2=0;
    for Af1=AmpFreqVector
        counter2=counter2+1;
        Af2=Af1+AmpFreq_BandWidth;
        
[MI,MeanAmp]=ModIndex_v1(lfp,srate,Pf1,Pf2,Af1,Af2,position);
Comodulogram(counter1,counter2)=MI;

    end
end

toc

%% Plot comodulogram

clf
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulogram',30,'lines','none')
set(gca,'fontsize',14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
colorbar




