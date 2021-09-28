function y = filtSRP(x,fs)
% filtSRP   Saccade-Related Potential (SRP) filter
% Filters the columns of x with a matched filter composed of the average
% normalized SRP of 5 subjects at the Radial Electro-Oculo-Gram channel 
% ( REOG = (HEOGL+HEOGR+VEOGS+VEOGI)/4 - Pz ) at a sampling rate of 1024 Hz.
% If the signal is sampled at a rate fs other than 1024 Hz, it is resampled
% to 1024 Hz.

%       Author(s): A.S. Keren, 8-24-09
%       Leon Deouell's Cognitive Neuroscience Laboratory
%       The Interdisciplinary Center for Neural Computation (ICNC)
%       The Hebrew University of Jerusalem (HUJI)
%       alon.keren@mail.huji.ac.il

% downloaded from https://www.hcnl.org/resources on Sept 18, 2021 by O.D.

resamp = false;
if nargin > 1 && fs ~= 1024
    resamp = true;
    x = resample(x,1024,fs);
end
% Set the SRP matched filter reversed impulse response 
% (84 samples @ 1024 Hz):
SRPfilt = ...
    [         0   -0.0000   -0.0001   -0.0002   -0.0002   -0.0001    0.0001    0.0003    0.0007    0.0015    0.0028 ...
    0.0050    0.0080    0.0114    0.0151    0.0188    0.0217    0.0241    0.0267    0.0272    0.0271    0.0287 ...
    0.0329    0.0391    0.0462    0.0544    0.0605    0.0602    0.0447    0.0030   -0.0672   -0.1615   -0.2631 ...
   -0.3490   -0.3965   -0.3834   -0.3045   -0.1706   -0.0109    0.1349    0.2355    0.2789    0.2707    0.2271 ...
    0.1683    0.1100    0.0631    0.0319    0.0174    0.0142    0.0193    0.0274    0.0312    0.0303    0.0257 ...
    0.0183    0.0088   -0.0007   -0.0086   -0.0152   -0.0198   -0.0221   -0.0229   -0.0230   -0.0219   -0.0199 ...
   -0.0179   -0.0157   -0.0129   -0.0101   -0.0070   -0.0042   -0.0020   -0.0003    0.0009    0.0013    0.0013 ...
    0.0011    0.0008    0.0005    0.0002    0.0001    0.0000         0      ]';

Nt = length(SRPfilt); % Number of time-samples in the filter

% Onset of the saccadic Spike Potential relative to SRPfilt onset:
% SPonset = 28;
SPonset = 35; % changed by Olaf Dimigen*

% [**] comment OD (Sept, 2021): latency of first EOG-negative spike is 35
% note: in general, the onset of EMG muscle activity is not necessarily 
% the onset of the eye ball movment

% Convolve the filter's impulse response with the signal:
y = conv(x,SRPfilt(end:-1:1)); 

% Crop convolution edges:
y = y((Nt:end)-SPonset+1,:);   

% Resample in the original rate
if resamp
    y = resample(y,fs,1024); 
end