% [MI,MeanAmp]=ModIndex_v2(Phase, Amp, position)
%
% Phase-amplitude cross-frequency coupling measure:
%
% Inputs:
% Phase = phase time series
% Amp = amplitude time series
% position = phase bins (left boundary)
%
% Outputs:
% MI = modulation index (see Tort et al PNAS 2008, 2009 and J Neurophysiol 2010)
% MeanAmp = amplitude distribution over phase bins (non-normalized)
 
function [MI,MeanAmp]=ModIndex_v2(Phase, Amp, position)

nbin=length(position);  
winsize = 2*pi/nbin;
 
% now we compute the mean amplitude in each phase:
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
