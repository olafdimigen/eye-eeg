% detecteogspikes - estimate likely (micro)saccade onsets based on spikes
%                   in the filtered "radial" EOG
%
% Usage:
% >> detecteogspikes(EEG,chans_eog,chans_parietal,filtermethod,spikethreshold,mindist,clustermode,add_reogchan,plotfig,writesac)
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
%  add_reogchan    0: only detect events, do not add channel
%                  1: add (filtered) rEOG as channel to EEG.data [not yet
%                  implemented!)
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

function EEG = detecteogspikes(EEG,chans_eog,chans_parietal,filtermethod,threshfactor,mindist,clustermode,add_reogchan,plotfig,writesac)

if isempty(chans_eog)
    error('\n%s(): you did not provide an index for EOG channel(s)',filename)
end
if isempty(chans_parietal)
    error('\n%s(): you did not provide an index for parietal channel(s)',filename)
end

% which method is used for filtering the rEOG?
switch filtermethod
    case 1
        fmethod = 'Keren et al. template';
    case 2
        fmethod = 'First derivative';
    otherwise
        error('filtermethod not recognized')
end

fprintf('\n%s(): Detecting spikes in radial EOG montage',mfilename)
fprintf('\n---------------------------------------------------')
fprintf('\nEOG channels used:                 %s',strjoin({EEG.chanlocs(chans_eog).labels}))
fprintf('\nParietal channels used:            %s',strjoin({EEG.chanlocs(chans_parietal).labels}))
fprintf('\n---------------------------------------------------')
fprintf('\nrEOG filter method:                %s',fmethod)
fprintf('\nThreshold multiplier (SDs):        %.2f',threshfactor)
fprintf('\nMin. dist. between spikes (in ms): %.2f',mindist*1000/EEG.srate)
fprintf('\nClustermode:                       %i',clustermode)
fprintf('\n---------------------------------------------------')

% number of epochs & time points
nepochs = size(EEG.data,3);
nsample = size(EEG.data,2);

% compute threshold for spike detetion across all epochs
% rather than separately for each epoch (likely more robust)
globaltresh = true;
if globaltresh
    switch filtermethod
        case 1
            greog = double(mean(EEG.data(chans_eog,:),1) - mean(EEG.data(chans_parietal,:),1));
            greog_filt = filtSRP(greog', EEG.srate);
            gthresh = threshfactor .* std(greog_filt);
        case 2
            error('Filtermethod 2 not yet implemented, please choose method 1')
        otherwise
            error('Filter method (filtermethod) unknown')
    end
end

% for epoched datasets, create concatenated rEOG of individual epochs
% just for the purpose of plotting (slow!)
if plotfig && nepochs > 1
    reogf_concat = [];
end

allpeaks = [];
for e = 1:nepochs
    
    switch filtermethod
        case 1
            
            %% compute rEOG
            % mean of all EOG channel(s) minus mean of all parietal channel(s)
            reog = double(mean(EEG.data(chans_eog,:,e),1) - mean(EEG.data(chans_parietal,:,e),1));
            
            %% filter rEOG by convolution with time-reversed SP template
            % using function filtSRP by Keren et al., 2010, NeuroImage
            %             if ~exist('resample','file') % check for SP add-on toolbox
            %                 warning('\n%s: Filtermethod 1 uses function resample() of the Signal Processing Toolbox',mfilename)
            %                 warning('It looks like you do not have the Signal Processing Toolbox installed')
            %                 warning('This will cause an error. Please try another filtermethod.')
            %             end
            reog_filt = filtSRP(reog', EEG.srate);
                        
        case 2
            error('Filtermethod 2 not yet implemented, please choose method 1')
        otherwise
            error('Filter method (filtermethod) unknown')
    end
    
    
    %% determine threshold
    if globaltresh
        thresh = gthresh; % compute across all epochs (if epoched data)
    else
        thresh = threshfactor .* std(reog_filt);
    end
    % alternative here: use median-based SD more robust against outliers
    % medianbased_sd = sqrt( median(reog_filt.^2) - (median(reog_filt) )^2 );
    % thresh = threshfactor .* medianbased_sd;
    
    %% find values exceeding threshold (spikes)
    % NOTE: we only search for positive spikes here, an alternative could
    % be to search for spikes in the absolute signal, leading to more
    % "double spikes"
    ix = find(reog_filt > thresh);
    % ix = find(abs(reog_filt) > thresh);
    
    % find adjacent above-threshold values
    if ~isempty(ix)
        peaks = findsequence2(ix);
        
        % set spike to largest peak
        % among successive above-threshold values, define peaks as sample
        % with highest rEOG value
        for p = 1:size(peaks,1)
            
            % if "peak" lasts several samples, take sample with maximum rEOG value
            if peaks(p,3) > 1
                nsmp = peaks(p,3);
                [mx,mxix] = max( reog_filt(peaks(p,1):peaks(p,2)) ); % get maximum
                %[mx,mxix] = max( abs(reog_filt(peaks(p,1):peak(p,2))) );
                peaks2(p,1) = peaks(p,1)+mxix-1;
                peaks2(p,2) = mx; % spike amplitude
            else % "peak" only lasts 1 sample
                peaks2(p,1) = peaks(p,1);
                peaks2(p,2) = reog_filt(peaks(p,1)); % spike amplitude
            end
        end
    else
        peaks2 = [];
    end
    
    if ~isempty(peaks2)
        
        %% merge clusters of detected spikes
        % we use same logic as in function "mergesacc" for eye tracker-based
        % saccade detection, following code is taken from there, help mergesacc()
        nspk = size(peaks2,1);
        keep = true(nspk,1); % vector of spikes to keep
        
        % determine temporal distance to preceding rEOG-spike
        if nspk > 1
            tempdist = [NaN; peaks2(2:end,1)-peaks2(1:end-1,1)];
        else
            tempdist = NaN;
        end
        
        % get index of spikes that form temporal clusters
        ix_close = find(tempdist < mindist);
        
        if ~isempty(ix_close)
            % get beginning, end, length of spike clusters
            clustlist      = findsequence2(ix_close);
            clustlist(:,1) = clustlist(:,1)-1; % include first saccade
            clustlist(:,3) = clustlist(:,3)+1; % update cluster length
            
            % go tru spike clusters
            for n = 1:size(clustlist,1)
                
                ixx = clustlist(n,1):clustlist(n,2);
                
                % how to treat spike clusters?
                switch clustermode
                    case 1
                        % do nothing, keep all spikes
                    case 2
                        % keep first spike of each cluster
                        keep(ixx(2):ixx(end)) = 0;
                    case 3
                        % keep largest spike of each cluster
                        keep(ixx) = 0;
                        [tmp,ix_max] = max(peaks2(ixx,2));
                        keep(ixx(ix_max)) = 1;
                end
            end
        end
        
        % remove "other" spikes of cluster
        if ismember(clustermode,2:3)
            % remove other saccades of cluster
            %fprintf('\n---------------------------------------------------');
            %fprintf('\nNumber of spikes before merging spike-clusters: %i',size(peaks2,1))
            peaks2 = peaks2(keep,:);
            %fprintf('\nNumber of spikes after merging spike-clusters:  %i',size(peaks2,1))
            %fprintf('\n---------------------------------------------------');
        end
        
        % add epoch number
        peaks2(:,3) = e;
        
        %% recompute latencies of spike events relative to entire dataset
        % (only necessary for epoched data)
        offset = (e-1)*nsample;
        peaks2(:,1) = peaks2(:,1)+offset;
        
    end
    
    allpeaks = [allpeaks; peaks2]; % slow but easy :)
    clear peaks* % important
    
    
    if plotfig && nepochs > 1
        reogf_concat = [reogf_concat; reog_filt]; % concatenated trials
    end
end

%% control figure
if plotfig
    
    %
    if nepochs == 1
        reogf_concat = reog_filt;
    end
    
    figure('name','rEOG detection'); hold on;
    title(sprintf('Results: rEOG-based saccade onset estimation: %i spikes detected',size(allpeaks,1)));
    % filtered reog
    h1 = plot(reogf_concat,'b','linewidth',0.7);
    yl = ylim;
    axis tight
    
    % threshold as green line
    h2 = plot([1 length(reogf_concat)],[thresh thresh],'g--','linewidth',1.2);
    xlabel('Time (in sampling points)')
    ylabel('Filtered radial EOG (rEOG)')
    
    if ~isempty(allpeaks)
        % detected spikes as red dots
        spikyspikes = NaN(1,length(reogf_concat));
        spikyspikes(allpeaks(:,1)) = reogf_concat(allpeaks(:,1));
        h3 = plot(spikyspikes,'r.');
        legend([h1 h2 h3 ],{'filtered rEOG','threshold','detected spike'},'location','best')
    else
        legend([h1 h2],{'filtered rEOG','threshold'},'location','best') % no spikes detected
    end
    
end
fprintf('\nNumber of detected rEOG spikes:   %i\n',size(allpeaks,1));

%% write new "rEOG saccade" events to EEG.event
if writesac && ~isempty(allpeaks)
    colnames      = {'latency','reog_spikeamplitude','epoch'};
    newevents     = [allpeaks(:,1) allpeaks(:,2) allpeaks(:,3)];
    eventname     = 'reog_spike';
    % add events
    EEG = addevents(EEG,newevents,colnames,eventname);
    
    fprintf('\n%i ''reog_spike'' events added to EEG.event',size(allpeaks,1));
    fprintf('\n---------------------------------------------------\n');
else
    fprintf('\nNo events were added to EEG.event (test mode, writesac == 0).');
    fprintf('\n---------------------------------------------------\n');
end

%% add channel with (filtered) radial EOG to EEG.data
if add_reogchan
   warning('\n%s:adding the filtered rEOG as a channel to EEG.data is not yet implemented',mfilename) 
end
    
