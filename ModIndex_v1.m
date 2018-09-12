% Programmed by Adriano Tort, CBD, BU, 2008
% 
% Phase-amplitude cross-frequency coupling measure:
%
% See Tort et al J Neurophysiol 2010 and the Supp Info of Tort et al PNAS 2008
% 
% [MI,MeanAmp] = ModIndex_v1(lfp,srate,Pf1,Pf2,Af1,Af2,position)
%
% Inputs:
% 
% lfp is the data in a vector format
% 
% srate is the sampling frequency
% 
% Pf1 and Pf2 define the frequency range (in Hz) investigated as the
% "phase-modulating" (for example, for theta take Pf1=6 and Pf2=12)
% 
% Af1 and Af2 define the frequency range investigated as the "amplitude
% modulated" by the phase frequency (e.g., low gamma would be Af1=30 Af2=55)
% 
% position = phase bins (left boundary)
%
% Outputs:
% MI = modulation index
% MeanAmp = Amplitude distribution per phase bin (non-normalized); to
% normalize, do MeanAmp = MeanAmp/sum(MeanAmp)
 
function [MI,MeanAmp] = ModIndex_v1(lfp,srate,Pf1,Pf2,Af1,Af2,position)

% the eegfilt routine employed below is obtained from the EEGLAB toolbox 
% (Delorme and Makeig J Neurosci Methods 2004)

PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % this is just filtering 
 
Phase=angle(hilbert(PhaseFreq)); % this is getting the phase time series
 
AmpFreq=eegfilt(lfp,srate,Af1,Af2); % just filtering
 
Amp=abs(hilbert(AmpFreq)); % getting the amplitude envelope
 
% Now we search for a Phase-Amp relation between these frequencies by
% caclulating the mean amplitude of the AmpFreq in each phase bin of the
% PhaseFreq
 
% Computing the mean amplitude in each phase:

nbin=length(position);
winsize = 2*pi/nbin;

MeanAmp=zeros(1,nbin); 
for j=1:nbin   
I = find(Phase <  position(j)+winsize & Phase >=  position(j));
MeanAmp(j)=mean(Amp(I)); 
end
% the center of each bin (for plotting purposes) is position+winsize/2
 
% quantifying the amount of amp modulation by means of a
% normalized entropy index (Tort et al PNAS 2008):

MI=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);

end
