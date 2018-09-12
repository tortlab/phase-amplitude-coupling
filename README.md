# phase-amplitude-coupling
Matlab routines for assessing phase-amplitude coupling by means of the modulation index (MI) and comodulogram, as described in Tort et al., J Neurophysiol 2010.

To get started, run the CallerRoutine.m routine using cell mode in Matlab. This routine will compute the comodulogram for a synthetic LFP example. The file LFP_HG_HFO.mat provides two actual LFP signals (sampled at 1000 Hz); one has prominent theta-HG coupling, the other strong theta-HFO coupling (see Scheffer-Teixeira et al., Cereb Cortex 2012). For examining the actual LFPs, uncomment lines 21-26 in CallerRoutine.m

These routines use the function eegfilt.m from the EEGLab Toolbox (Delorme & Makeig, J Neurosci Methods 2004).
