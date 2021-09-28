% pop_detecteogspikes - estimate likely (micro)saccade onsets based on
%                     spikes in the filtered "radial EOG
%
% Usage:
% >> pop_detecteogspikes(EEG,chans_eog,chans_parietal,filtermethod,spikethreshold,mindist,clustermode,add_reogchan,plotfig,writesac)
%
% Required inputs:
%   EEG             A continuous or epoched EEG dataset fullfilling three criteria:
%                     1. not yet corrected for ocular artifacts
%                     2. dataset contains one or more (EOG) electrode(s) near the eyes,
%                        which were recorded against the common reference
%                        (i.e., not with a bipolar EOG montage)
%                     3. still contains higher signal frequencies (~30-80 Hz),
%                        (i.e, no strong low-pass filtering was performed yet)
%
%   chans_eog       Index of all facial EOG electrodes. This should include
%                   all non-reference electrodes place on the face, nose/nasion, below
%                   the eyes or near the left or right canthus. Electrodes
%                   above the eyes (Fp1, FP2) are less effective, but
%                   might be included if no other external EOG electrodes are
%                   available. Do not include channels that were recorded
%                   with a bipolar montage. At least one channel index must
%                   be provided here.
%
%   chans_parietal  Index of channel(s) near the parietal/parieto-
%                   occipital midline. If "Pz" is available, it is suggested to
%                   just use "Pz" here (see Keren et al., 2010, NeuroImage).
%                   Otherwise, locations like P1/P2/POz/Oz might also work,
%                   alone or in combination. At least one channel index
%                   must be provided here.
%
%   filtermethod    method to filter the rEOG before detecting peaks:
%                   1: convolution with spike potential template by Keren et al.
%                      note: this requires the Signal Processing Toolbox
%                      (!!no other methods currently implemented!!)
%
%   threshfactor    Threshold for spike detection, expressed in standard deviations
%                   of the filtered rEOG signal. Recommended: ~2.0 SDs
%
%   mindist        Minimum temporal distance between successive spikes (in samples)
%                  {experimental parameter, may need fine-tuning}
%
%   clustermode    How to deal with "clusters" of detected spikes closer to
%                  each other than "mindist" samples?
%                    1: keep all spikes (do nothing)
%                    2: keep only first spike of cluster
%                    3: keep only largest spike of cluster (probably best)
%                  {experimental parameter, may need fine-tuning}
%
%  add_reogchan    [0/1]: add rEOG signal to EEG.data?
%                  0: only detect spikes events, do not add channel
%                  1: add (filtered) rEOG as extra channel to EEG.data 
%                    {!!option "1" not yet implemented!!}
%
%  plotfig         [0/1] Plot figure showing filtered rEOG and detected
%                  spikes (strongy recommended!)
%
%  writesac        [0/1]: Add detected rEOG spike-events to EEG.event?
%                    0: detect rEOG-spikes, but do not add them to EEG.event
%                       (recommended to test detection parameters first)
%                    1: add detected rEOG-spikes as new events to EEG.event
%
% Outputs:
%                  EEG with additional "reog_spike" events in EEG.event
%
% Example call:
%    EEG = detecteogspikes(EEG,[1 2 3 4],[33 44],1,2.0,25,3,0,1,0)
%
% This example call defines channels 1 to 4 as EOG electrodes, and channel 33
% as the parietal electrodes for the rEOG computation. The rEOG
% signal will then be filtered by convolution with a time-reversed template
% of a typical spike potential (subfunction filtSRP.m by A.S. Keren). Here, the
% threshold for spike detection is set at 2.0 SDs (threshfactor =2) of the
% filtered rEOG signal. We require spikes to be > 25 sampling points
% apart from each other (mindist = 25). If spikes are closer to each
% other, we'll only keep the largest spike (clustermode = 3) of the cluster.
% We do not want to add the filtered (e.g. template-convolved) rEOG as a
% new channel to EEG.data. We want to plot a figure showing the detection
% results (plotfig = 1), but not yet write the detected spikes as events to
% EEG.event (writesac = 0).
%
% #########################################################################
% This function uses the subfunction "filtSRP()" written by A.S. Keren.
% If you use this function, please make sure to cite the following paper,
% which also explains the details of rEOG-based saccade onset detection:
%
% Keren, A. S., Yuval-Greenberg, S., & Deouell, L. Y. (2010).
% Saccadic spike potentials in gamma-band EEG: characterization, detection
% and suppression. Neuroimage, 49(3), 2248-2263.
% Subfunction included with permission.
%
% #########################################################################
%
% Author: od
% Copyright (C) 2009-2021 Olaf Dimigen, HU Berlin
% olaf@dimigen.de
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, 51 Franklin Street, Boston, MA 02110-1301, USA

function [EEG, com] = pop_detecteogspikes(EEG,chans_eog,chans_parietal,filtermethod,threshfactor,mindist,clustermode,add_reogchan,plotfig,writesac)

com = '';

if nargin < 1
    help(mfilename);
    return;
end

try
    if nargin < 10
        % pop up dialogue
        [chans_eog,chans_parietal,threshfactor,mindist,clustermode,plotfig,writesac] = dlg_detecteogspikes(mfilename, EEG.chanlocs, EEG.srate);
        
        % features not yet implemented (in GUI)
        filtermethod = 1;
        add_reogchan = 0;
    end
    
    % detect eog spikes
    EEG = detecteogspikes(EEG,chans_eog,chans_parietal,filtermethod,threshfactor,mindist,clustermode,add_reogchan,plotfig,writesac);
    
catch err
    if (strcmp(err.identifier,'MATLAB:unassignedOutputs'))
        return
    else
        rethrow(err);
    end
end

% return history command string
allArgs = vararg2str({chans_eog,chans_parietal,filtermethod,threshfactor,mindist,clustermode,add_reogchan,plotfig,writesac});
com = sprintf('EEG = %s(EEG,%s)',mfilename,allArgs);
return
end