% detecteyemovements() - detect saccades & fixations in eye tracking data. 
%                Saccade detection is based on the algorithm by 
%                Engbert & Mergenthaler (2006). Saccades are defined as 
%                (monocular or binocular) outliers in 2D velocity space.
%                Velocity thresholds for saccade detection are determined 
%                adaptively as a multiple of the (median-based) SD of all 
%                data samples in the epoch. Fixations are defined as 
%                intervals in-between saccades. Eye movements can be added 
%                as new events to EEGLAB's event structure. For various 
%                other options, see below.
%
% Usage:
%   >> EEG = detecteyemovements(EEG,left_eye_xy,right_eye_xy,vfac,mindur,...
%            degperpixel,smooth,globalthresh,clusterdist,clustermode,
%            plotfig,writesac,writefix)
%
% Required inputs:
%   EEG          - [string] EEG struct, also containing synchronized eye 
%                  tracking data (see pop_importeyetracker)
%   left_eye_xy  - [vector of two channel indices], 
%                  specifying channel indices of X- (first value) and 
%                  Y-component (second value) of left eye. Leave empty []
%                  if the left eye was not recorded.
%   right_eye_xy - [vector of two channel indices], 
%                  specifying channel indices of X- (first value) and 
%                  Y-component (second value) of right eye. Leave empty []
%                  if the right eye was not recorded.
%   vfac         - [double] velocity factor ("lambda") to determine 
%                  the velocity threshold for saccade detection 
%                  (cf. Engbert & Mergenthaler, 2006)
%   mindur       - [integer] minimum saccade duration (in samples)
%                  (cf. Engbert & Mergenthaler, 2006)
%   degperpixel  - [double] visual angle of one screen pixel
%                  if this value is left empty [], saccade characteristics 
%                  are reported in the original data metric (pixel?) 
%                  instead of in degrees of visual angle
%   smooth       - [0/1] if set to 1, the raw data is smoothed over a 
%                  5-sample window to suppress noise
%                  noise. Recommended for high native ET sampling rates.
%   globalthresh - [0/1]. Use the same thresholds for all epochs? 
%                  0: Adaptive velocity thresholds are computed 
%                  individually for each data epoch. 
%                  1: Adaptive velocity thresholds are first computed for 
%                  each epoch, but then the mean thresholds are applied to 
%                  each epochs (i.e. same detection parameters are used for
%                  all epochs). Setting is irrelevant if the input data is 
%                  still continuous (= only one data epoch).
%   clusterdist  - [integer] value in sampling points that defines the
%                  minimum allowed fixation duration between two saccades. 
%                  If the off- and onsets of two temp. adjacent sacc. are 
%                  closer together than 'clusterdist' samples, these 
%                  saccades are regarded as a "cluster" and treated 
%                  according to the 'clustermode' setting (see below).
%                  clusterdist is irrelevant if clustermode == 1.
%   clustermode  - [1,2,3,4]. Integer between 1 and 4.
%                  1: keep all saccades, do nothing
%                  2: keep only first saccade of each cluster
%                  3: keep only largest sacc. of each cluster
%                  4: combine all movements into one (longer) saccade
%                     this new saccade is defined as the movement that 
%                     occurs between the onset of the 1st saccade in the
%                     cluster and the offset of the last sacc. in cluster
%                     WARNING: CLUSTERMODE 4 is experimental and untested!  
%   plotfig      - [0/1] Show a figure with eye movement properties?
%                  0: do not plot a figure. 
%                  1: plot a figure displaying properties of detected 
%                  saccades & fixations
%   writesac     - [0/1]: Add saccades to EEG.event?
%                  0: detect saccades, but do not store them in EEG.event. 
%                  1: add detected saccades as new events to EEG.event.
%   writefix     - [0/1]: Add fixations to EEG.event?
%                  0: detect fixations, but do not add them to EEG.event.
%                  1: add detected fixations as new events to EEG.event.
%                  Note: It is recommended to first test the parameters of
%                  saccade/fixation detection without adding events.
%                  For this, set writesac and writefix to 0.
%
% Outputs:
%   EEG         - EEG structure. If writesac or writefix were set to 1, 
%                 the EEG structure (EEG.event/EEG.urevent/EEG.epoch)
%                 will contain additional "saccade" and "fixation" events
%                 with their respective properties
%
% See also: vecvel, velthresh, microsacc_plugin, binsacc, saccpar, 
%           mergesacc, addevents
%
%
% MAJOR UPDATE 03/2018: bad ET intervals, as marked by "bad_ET" events
% in EEG.event (see function pop_rej_contin.m) are excluded from eye movement 
% detection. Bad intervals are now ignored when estimating velocity thresholds
% and saccade/fixation events overlapping with these intervals are removed
%
% An example call of the function might look like this: 
% >> EEG = detecteyemovements(EEG,[],[33 34],6,4,0.037,1,0,25,4,1,1,0)
%
% In this example, the eye position data for the right eye is stored in 
% channels 33 (horiz.) and 34 (vertical). The left eye was not recorded. 
% The velocity threshold is set to 6 times the (median-based) 
% SD of all velocity samples in the epoch. The minimum duration of
% saccades to be detected is 4 samples. In the experiment, one screen 
% pixel corresponded to 0.037 degrees of visual angle. 
% The raw data is smoothed prior to saccade detection (smooth: 1). 
% Adaptive velocity thresholds (X and Y-threshold for each eye) are 
% determined individually for each data epoch (globalthresh: 0). For saccades 
% separated by fixations of less than 25 samples, only the first saccade 
% is kept (clusterdist: 25, clustermode: 2). A figure with the 
% results is plotted. Detected saccades are stored as new events in 
% EEG.event, but fixations are not stored.
% 
% The eye movement detection is based on:
%
% Engbert, R., & Kliegl, R. (2003). Microsaccades uncover the orientation
% of covert attention. Vision Research, Vol. 43, 1035-1045
%
% ...as well as...
%
% Engbert, R., & Mergenthaler, K. (2006). Microsaccades are triggered by 
% low retinal image slip, PNAS, Vol. 103 (18), 7192-7197
%
% Author: od
% Copyright (C) 2009-2018 Olaf Dimigen, HU Berlin
% olaf.dimigen@hu-berlin.de 

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

function EEG = detecteyemovements(EEG,left_eye_xy,right_eye_xy,vfac,mindur,degperpixel,smooth,globalthresh,clusterdist,clustermode,plotfig,writesac,writefix)

allsac = [];
allfix = [];

% data of which eye is available?
ldata  = false;
if length(left_eye_xy) == 2, ldata = true; end
rdata  = false;
if length(right_eye_xy) == 2, rdata = true; end

if length(left_eye_xy) == 1 || length(right_eye_xy) == 1 
    error('%s(): For each recorded eye, horizontal (X) and vertical (Y) gaze channel must be specified.',mfilename);
end

nsample   = size(EEG.data,2);
nepochs   = size(EEG.data,3);
badepochs = zeros(nepochs,1);
nbadsmp   = 0;

% preallocate storage for sacc. detection thresholds
l_msdx = NaN(nepochs,1);
l_msdy = NaN(nepochs,1);
r_msdx = NaN(nepochs,1);
r_msdy = NaN(nepochs,1);

% warn message if back-to-back saccades are detected
clusterwarning = false;

% critical bugfix 2013-10-01, by OD
% due to bug in third-party function "smoothdata":
% function smoothdata() was removed from the toolbox
% function vecvel() was updated to incorporate different levels of smoothing:
% options:
% - smoothlevel 0: no smoothing, simple diff()
% - smoothlevel 1: 3-point window
% - smoothlevel 2: 5-point window
if smooth
    smoothlevel = 2; % 5-point smoothing
else
    smoothlevel = 0; % no smoothing
end


%% screen feedback

% detection feedback
if ~ldata && ~rdata
    error('%s(): Please correctly specify channels containing eye tracking data.\nFor each recorded eye, both horiz. (X) and vert. (Y) channel must be specified.',mfilename)
else
    % summary of detection parameters
    fprintf('\n--------------------------------------------------------------------')
    fprintf('\nDetecting saccades after Engbert & Mergenthaler (2006)\n')
    fprintf('\nVelocity threshold factor (vfac):  %.2f SD',vfac);
    fprintf('\nMinimum saccade duration (mindur): %.2f samples (%.2f ms)',mindur,mindur*1000/EEG.srate);
    if ~isempty(degperpixel) | isnan(degperpixel) % bugfix 2016-11-12 by OD: added case if degperpixel = NaN (from GUI input) 
        fprintf('\nVisual angle per screen pixel:     %f°',degperpixel);
        metric = 'deg';
    else
        fprintf('\nWARNING: No input provided for degperpixel!\nSpatial saccade properties are given in original metric (pixel?)');
        degperpixel = 1;
        metric = 'pix';
    end
    if nepochs < 2, fprintf('\n-- Using continuous data.'), else fprintf('\n-- Using epoched data.'); end
    if ldata && rdata, fprintf('\n-- Using binocular data:'); else fprintf('\n-- Using monocular data:'); end
    if ldata
        fprintf('\n\tLeft horiz.: \"%s\"',EEG.chanlocs(left_eye_xy(1)).labels);
        fprintf('\n\tLeft verti.: \"%s\"',EEG.chanlocs(left_eye_xy(2)).labels);
    end
    if rdata
        fprintf('\n\tRight horiz.: \"%s\"',EEG.chanlocs(right_eye_xy(1)).labels);
        fprintf('\n\tRight verti.: \"%s\"',EEG.chanlocs(right_eye_xy(2)).labels);
    end    
    if smooth, fprintf('\n-- Raw data is smoothed in 5-sample window.'); else fprintf('\n-- Raw data is not smoothed.'); end
    if nepochs > 1
        if globalthresh
            fprintf('\n-- Velocity thresholds computed globally across all %i epochs',nepochs);
        else
            fprintf('\n-- Velocity thresholds computed individually for each epoch');
        end
    end
    fprintf('\n-- Treatment of saccade clusters:');
    switch clustermode    
        case 1
            fprintf('\n\tAll saccades are kept');
        case 2
            fprintf('\n\tSaccades separated by fixations < %i ms are clustered.',clusterdist*(1000/EEG.srate));
            fprintf('\n\tFirst sacc. of each cluster is kept.');
        case 3
            fprintf('\n\tSaccades separated by fixations < %i ms are clustered.',clusterdist*(1000/EEG.srate));
            fprintf('\n\tLargest sacc. of each cluster is kept.');
        case 4
            fprintf('\n\tSaccades separated by fixations < %i ms are clustered.',clusterdist*(1000/EEG.srate));
            fprintf('\n\tClusters are combined into one saccade.');
        otherwise
            error('%s(): Unknown input for clustermode, should be: 1,2,3,4.',mfilename)
    end
    if plotfig, fprintf('\n-- A figure with eye movement properties is plotted.'); end        
    if writesac && writefix, fprintf('\n-- Saccades and fixations will be added to EEG.event.');
    elseif writesac,         fprintf('\n-- Saccades will be added to EEG.event.');
    elseif writefix,         fprintf('\n-- Fixations will be added to EEG.event.');
    else                     fprintf('\n-- Saccades and fixations are detected, but NOT stored anywhere.');
    end
end

%% get "bad_ET" intervals from EEG.event structure
badvector = zeros(1,size(EEG.data,2)*size(EEG.data,3));
ix_badETevent = ismember({EEG.event.type},'bad_ET');
if ~isempty(ix_badETevent)
    bad_lat     = [EEG.event(ix_badETevent).latency];
    bad_dur     = [EEG.event(ix_badETevent).duration];
    bad_ET      = [bad_lat; bad_dur]';
    bad_ET(:,3) = bad_ET(:,1)+bad_ET(:,2)-1;
    
    % create long vector (as long as EEG) indicating bad samples
    for j = 1:size(bad_ET,1)
        badvector(bad_ET(j,1):bad_ET(j,3)) = 1;
    end
end
% reshape "badvector" to 3D if data is already epoched
badvector = reshape(badvector,1,size(EEG.data,2),size(EEG.data,3));



%% pre-compute saccade velocity thresholds for all epochs
% this enables the option to use the same threshold for all epochs, not
% provided by the original Engbert & Mergenthaler implementation
% new 03/2018: exclude bad ET intervals from velocity estimation
fprintf('\nComputing adaptive velocity thresholds...')
% if any(badvector(:))
%     fprintf('\n-- Found %i \"bad_ET\" events marking bad eye-tracking intervals in EEG.event.',length(ix_badETevent))
%     fprintf('\n-- These intervals (%.2f%% of data) will be ignored when computing velocity thresholds.',(sum(badvector)/length(badvector(:))*100))
% end
for e=1:nepochs
    ix_goodET = ~badvector(:,:,e); % get index of non-bad ET samples
    if ldata
        l = EEG.data([left_eye_xy(1) left_eye_xy(2)],ix_goodET,e)';
        vl = vecvel(l,EEG.srate,smoothlevel);
        [l_msdx(e), l_msdy(e)] = velthresh(vl);
    end
    if rdata
        r = EEG.data([right_eye_xy(1) right_eye_xy(2)],ix_goodET,e)';
        vr = vecvel(r,EEG.srate,smoothlevel);
        [r_msdx(e), r_msdy(e)] = velthresh(vr);
    end
end


%% detect saccades & fixations
for e=1:nepochs
    
    %% saccades of left eye
    sac = [];
    if ldata
        l = EEG.data([left_eye_xy(1) left_eye_xy(2)],:,e)'; % don't exclude bad_ET samples here (otherwise discont. high-velocity jumps in time series)
        % bad/missing samples in eye track?
        badsmp = sum(sum(l<=0)); if badsmp > 0, badepochs(e) = 1; nbadsmp = nbadsmp + badsmp; end
        vl = vecvel(l,EEG.srate,smoothlevel); % get eye velocities
        % detect monocular saccades
        if globalthresh % use precomputed velocity thresholds (mean of all epochs)
            sacL = microsacc_plugin(l,vl,vfac,mindur,mean(l_msdx),mean(l_msdy));
        else % compute velocity thresholds from this epoch only
            sacL = microsacc_plugin(l,vl,vfac,mindur,l_msdx(e),l_msdy(e));        
        end
    end
    %% saccades of right eye
    if rdata
        r = EEG.data([right_eye_xy(1) right_eye_xy(2)],:,e)';
        % bad/missing samples in eye track?
        badsmp = sum(sum(r<=0)); if badsmp > 0, badepochs(e) = 1; nbadsmp = nbadsmp + badsmp; end        
        vr = vecvel(r,EEG.srate,smoothlevel); % get eye velocities
        % detect monocular saccades
        if globalthresh % use precomputed velocity thresholds (mean of all epochs)
            sacR = microsacc_plugin(r,vr,vfac,mindur,mean(r_msdx),mean(r_msdy));
        else 
            % compute velocity thresholds from this epoch only
            sacR = microsacc_plugin(r,vr,vfac,mindur,r_msdx(e),r_msdy(e));             
        end
    end
    
    %% binocular saccades
    if ldata && rdata
        [sac, tmp, tmp] = binsacc(sacL,sacR); clear tmp
        sac = saccpar(sac); % average saccade characteristics of both eyes
        sac = mergesacc(sac,(l+r)./2,clusterdist,clustermode); % merge nearby saccades (e.g. glissades)
    elseif ldata
        sac = sacL; clear sacL;
        sac = saccpar([sac sac]);
        sac = mergesacc(sac,l,clusterdist,clustermode);
    elseif rdata
        sac = sacR; clear sacR;
        sac = saccpar([sac sac]);
        sac = mergesacc(sac,r,clusterdist,clustermode);    
    end
    
    %% update various saccade metrics
    if ~isempty(sac)
        
        % define saccade duration as difference between saccade offset and
        % saccade onset sample. In saccpar(), monocular saccade durations
        % of both eyes are averaged, leading to uneven values (e.g.: 10.5
        % samples) different from the difference between onset and offset
        % values (which are the monocular extremes).
        % Instead: use difference between offset and onset
        sac(:,3) = sac(:,2)-sac(:,1)+1;
        
        % report saccade velocity/distance/amplitude as visual angles
        sac(:,[5 6 8]) = sac(:,[5 6 8]) .* degperpixel;
        
        % report saccade angles in degree rather than radians
        sac(:,[7 9]) = sac(:,[7 9]) * 180/pi;
        
        % add index of corresp. data epoch
        sac(:,10) = e;
        
        % store screen location for start and end of saccade
        if ldata && rdata
            gazexy = (l+r)./2; % binoc. recordings: average across eyes
        elseif ldata
            gazexy = l;
        elseif rdata
            gazexy = r;
        end
        % get position immediatly before sacc. onset and after sacc. offset
        startsmp = sac(:,1)-1; endsmp = sac(:,2)+1;
        if startsmp(1) < 1, startsmp(1,1) = 1; end
        if endsmp(end) > size(gazexy,2), endsmp(end) = size(gazexy,2); end
        sac(:,11) = gazexy(startsmp,1);
        sac(:,12) = gazexy(startsmp,2);
        sac(:,13) = gazexy(endsmp  ,1);
        sac(:,14) = gazexy(endsmp  ,2);
    end
    
    % columns of [sac]:
    % 1: saccade onset (sample)
    % 2: saccade offset (sample)
    % 3: duration (samples)
    % 4: delay between eyes (samples)
    % 5: vpeak (peak velocity)
    % 6: saccade "distance" (eucly. dist. between start and endpoint)
    % 7: saccade angle (based on saccade "distance")
    % 8: saccade "amplitude" (eucly. dist. of min/max in full saccade trajectory)
    % 9: saccade angle (based on saccade "amplitude")
    %10: index of corresponding data epoch (1 in case of contin. data)
    %11: horizontal (x) gaze position before start of saccade (pixel)
    %12: vertial (y) gaze position before start of saccade (pixel)
    %13: horizontal (x) gaze position after end of saccade (pixel)
    %14: vertical (y) gaze position after end of saccade (pixel)
    
    
    %% remove saccades that occured during "bad_ET" intervals [*]
    % Delete all saccades whos onset or offset occurs during a 
    % bad_ET interval (signal jumps to blinks are otherwise detected as 
    % saccades)
    %
    % [*] Note: It is not clear/trivial how to treat blinks in the context 
    % of EM detection, since we do not know what really happened during a
    % blink or loss of the signal. For example, is a long fixation with a 
    % blink in the middle really two fixations? Or should it be treated as 
    % one long fixation? Or should all fixations that are ended or started 
    % by blinks or are interrupted by them be completely removed?
    % The solution here is preliminary: I first remove all saccades that 
    % started or ended during a blink interval. The remaining saccades are 
    % used to define fixatinos. Finally, all fixations are removed that 
    % overlap with "bad_ET" intervals in EEG.event. So this approach 
    % removes any eye movement event that overlaps with a bad_ET interval. 
    % OD, 2018-03-08
        
    %% delete saccades starting/ending during "bad_ET" intervals
    badETsmp = find(badvector(:,:,e));
    if ~isempty(sac)
        ix_fakesac = find(ismember(sac(:,1),badETsmp) | ismember(sac(:,2),badETsmp));
        sac(ix_fakesac,:) = [];
    %         if ~isempty(ix_fakesac)
    %             fprintf('\n\n-- Removed %i saccades that occured during "bad_ET" intervals',length(ix_fakesac))
    %         end
    end
       
    %% get fixations
    nsac = size(sac,1);
    fix = [];
    if ~isempty(sac)
        if nsac > 1
            for f = 1:nsac-1
                
                fix(f,1) = sac(f,2)+1;
                fix(f,2) = sac(f+1,1)-1;

                % catch special case: if Engbert algorithms are applied 
                % without any saccade clustering [>> mergesacc()] there can
                % be back-to-back saccades with intervening "fixations" of
                % zero sample duration. Catch this by setting the duration
                % of these fixations to one sample
                if fix(f,1) > fix(f,2)
                    fix(f,1) = fix(f,2); % 1-sample fixation
                    clusterwarning = true;
                end
                
            end
        end
        % if epoch does not begin with saccade, add first fixation
        if sac(1,1) > 1
            fix = [[1 sac(1,1)-1]; fix];
        end
        % if epoch does not end with saccade, add last fixation
        if sac(end,2) < nsample
            fix = [fix;[sac(end,2)+1 nsample]];
        end
        
        % add more fixation properties
        for f = 1:size(fix,1)
            
            % fixation duration (in samples!)
            fix(f,3) = fix(f,2)-fix(f,1)+1;
            
            % mean fix. position: left eye
            if ldata
                fix(f,4) = mean( l(fix(f,1):fix(f,2),1) );
                fix(f,5) = mean( l(fix(f,1):fix(f,2),2) );
            else
                fix(f,4) = NaN;
                fix(f,5) = NaN;
            end
            % mean fix. position: right eye
            if rdata
                fix(f,6) = mean( r(fix(f,1):fix(f,2),1) );
                fix(f,7) = mean( r(fix(f,1):fix(f,2),2) );
            else
                fix(f,6) = NaN;
                fix(f,7) = NaN;
            end
            % binocular fixation position
            fix(f,8) = nanmean(fix(f,[4 6]));
            fix(f,9) = nanmean(fix(f,[5 7]));
        end
        fix(:,10) = e; % add index of corresp. data epoch
        
        % recompute latencies of eye movement events 
        % (only necessary for epoched data)
        offset = (e-1)*nsample;
        sac(:,[1 2]) = sac(:,[1 2])+offset;
        % special case: single sacc. lasts entire epoch
        if ~isempty(fix) 
            fix(:,[1 2]) = fix(:,[1 2])+offset;
        end
    end
    
    % columns of [fix]:
    % 1: fixation onset (sample)
    % 2: fixation offset (sample)
    % 3: duration (samples)
    % 4: mean fix position (L eye X)
    % 5: mean fix position (L eye Y)
    % 6: mean fix position (R eye X)
    % 7: mean fix position (R eye Y)
    % 8: mean fix position (L/R average X)
    % 9: mean fix position (L/R average Y)
    %10: index of corresponding data epoch (1 in case of contin. data)
    
    %% remove fixations overlapping with "bad_ET" intervals
    % this includes fixations interrupted by a blink (!)
    % development note: loop is slow, need to implement more efficiently
    if ~isempty(fix) & any(badvector(:))
        badfix = false(size(fix,1),1); % preallocate logical
        
        % go tru fixations, check whether "bad"
        for p = 1:size(fix,1)
            fixsmp = fix(p,1):fix(p,2);
            if any(ismember(fixsmp,badETsmp))
                badfix(p) = true;
            end
        end
        % remove fixations overlapping with "bad_ET" intervals
        fix(badfix,:) = [];
        %         if any(badfix)
        %             fprintf('\n-- Removed %i fixations that overlapped with "bad_ET" intervals',sum(badfix))
        %         end
    end
      
    % slow, but simple:
    allsac = [allsac;sac];
    allfix = [allfix;fix];
    
end % epoch loop


%% remove artificial eye movements caused by boundaries (data breaks)
% Remove all EMs whose onset is detected in temporal proximity to boundary. 
% Otherwise, data breaks will likely result in additional fake sacc./fix.
% Applies only if eye movements are detected in the continuous data
if nepochs == 1

    ix_bnd = find(cellfun(@(x) strcmp(x,'boundary'),{EEG.event.type})); % bug fix: now robust against numeric types
    
    % any data breaks?
    % (due to manual rejections or function pop_rej_eyecontin)
    if ~isempty(ix_bnd)
        
        % minimum distance from data break in milliseconds (hard-coded)
        BOUNDDIST_MS = 50;
        BOUNDDIST    = round(BOUNDDIST_MS / (1000/EEG.srate)); 
        % data break latencies
        bound_lats   = round([EEG.event(ix_bnd).latency]); 

        % mark all samples close to data break
        boundvector  = zeros(1,EEG.pnts);
        for b = 1:length(bound_lats)
            lowr = bound_lats(b)-BOUNDDIST;
            uppr = bound_lats(b)+BOUNDDIST;
            if lowr <= 0, lowr = 1; end
            if uppr > EEG.pnts, uppr = EEG.pnts; end
            boundvector(lowr:uppr) = 1;
        end
        nearboundsmp = find(boundvector);
        
        % option 1: event onset is close to boundary
        % fakesac = find(ismember(allsac(:,1),nearboundsmp));
        % fakefix = find(ismember(allfix(:,1),nearboundsmp));
        
        % option 2: event on- or offset is close to boundary
        fakesac = find(ismember(allsac(:,1),nearboundsmp) | ismember(allsac(:,2),nearboundsmp));
        fakefix = find(ismember(allfix(:,1),nearboundsmp) | ismember(allfix(:,2),nearboundsmp));

        allsac(fakesac,:) = [];
        allfix(fakefix,:) = [];
        fprintf('\n--------------------------------------------------------------------');
        fprintf('\nFound %i data breaks (boundary events) in the continuous data',length(ix_bnd));
        fprintf('\nRemoving eye movements that might be artifacts of data breaks:');
        fprintf('\nRemoved %i saccades  < %i ms away from a boundary',length(fakesac),BOUNDDIST_MS);
        fprintf('\nRemoved %i fixations < %i ms away from a boundary',length(fakefix),BOUNDDIST_MS);
    end
end

if clusterwarning && clustermode == 1
    fprintf('\n--------------------------------------------------------------------');
    fprintf('\n\n*** WARNING! ***\nDetected pairs of saccades that ended/started on adjacent data samples.');
    fprintf('\nYou should probably use saccade clustering! (clustermode: 2,3,4)');
end


%% user feedback if bad eye-tracking data (values <= 0) but no bad_ET events
if sum(badepochs)>0 & isempty(ix_badETevent)
    warning('\nI found bad or missing data (i.e. values <= 0) in the ET channels but no \"bad_ET" events!');
    fprintf('\nDetails:'); 
    fprintf('\n-- %i of %i epochs contained gaze position values <= 0',sum(badepochs),nepochs);
    fprintf('\n-- Total number of bad samples <= 0: %i',nbadsmp);
    fprintf('\n-- Did you detect or reject bad intervals with out-of-range values (''Reject data based on eyetrack'')?');
    fprintf('\n-- Did you subtract a baseline from eye channels (explaining negative values)?');
    fprintf('\n-- Note that unmarked blinks or missing data in the eye-track will...:')
    fprintf('\n    ...distort the velocity threshold for saccade detection!');
    fprintf('\n    ...be erroneously detected as additional saccades/fixations!');
end


%% user feedback: saccade & fixation detection
fprintf('\n--------------------------------------------------------------------');
fprintf('\nVelocity thresholds used:'); if nepochs > 1, fprintf(' (mean across epochs):'); end;
if ldata, fprintf('\n\tLeft eye.  Horiz.: %.2f %s/s. Vert.: %.2f %s/s',mean(l_msdx(e)*vfac*degperpixel),metric,mean(l_msdy*vfac*degperpixel),metric); end
if rdata, fprintf('\n\tRight eye. Horiz.: %.2f %s/s. Vert.: %.2f %s/s',mean(r_msdx(e)*vfac*degperpixel),metric,mean(r_msdy*vfac*degperpixel),metric); end
fprintf('\n--------------------------------------------------------------------')
if ~isempty(allsac)
    fprintf('\n%i saccades detected:',size(allsac,1));
    fprintf('\n\tMedian amplitude:  %.2f %s',median(allsac(:,6)),metric);
    fprintf('\n\tMedian duration:   %.2f ms',median(allsac(:,3))*1000/EEG.srate);
    fprintf('\n\tMedian peak veloc. %.2f %s/s',median(allsac(:,5)),metric);
end
if ~isempty(allfix)
    fprintf('\n%i fixations detected:',size(allfix,1));
    fprintf('\n\tMedian duration:   %.2f ms',median(allfix(:,3))*1000/EEG.srate);
    if ldata,fprintf('\n\tMedian fix. pos. left  eye: Horiz.: %.2f px. Vert.: %.2f px',median(allfix(:,4)),median(allfix(:,5))); end
    if rdata,fprintf('\n\tMedian fix. pos. right eye: Horiz.: %.2f px. Vert.: %.2f px',median(allfix(:,6)),median(allfix(:,7))); end    
end


%% plot figure with eye movements properties
if plotfig
    if ldata && rdata
        fprintf('\n--------------------------------------------------------------------')
        fprintf('\nPlotting eye movement properties...\nFor plotting, fixation locations are averaged across eyes.')
        ploteyemovements(allsac(:,6),allsac(:,5),allsac(:,7),allfix(:,3)*1000/EEG.srate,mean(allfix(:,[4 6]),2),mean(allfix(:,[5 7]),2),metric);
    elseif ldata
        fprintf('\n--------------------------------------------------------------------')
        fprintf('\nPlotting eye movement properties of left eye...')
        ploteyemovements(allsac(:,6),allsac(:,5),allsac(:,7),allfix(:,3)*1000/EEG.srate,allfix(:,4),allfix(:,5),metric);
    else
        fprintf('\n--------------------------------------------------------------------')
        fprintf('\nPlotting eye movement properties of right eye...')
        ploteyemovements(allsac(:,6),allsac(:,5),allsac(:,7),allfix(:,3)*1000/EEG.srate,allfix(:,6),allfix(:,7),metric);
    end
end

%% write eye movements to EEG.event & EEG.urevent
if writefix || writesac

    % check: are there already eye movements events in EEG.event?
    em_types = {'saccade','fixation','L_saccade','R_saccade','L_fixation','R_fixation'};
    if any(cellfun(@(x) any(strcmp(x,em_types)),{EEG.event.type})) % updated [v.0.337]
        fprintf('\n--------------------------------------------------------------------')
        fprintf('\n*** WARNING! ***:\nFound existing eye movement events in EEG.event!');
        fprintf('\nNew events will be added to the already existing events!');
    end
    
    % write saccades
    if writesac
        fprintf('\n--------------------------------------------------------------------')
        fprintf('\nAdding %i saccades to EEG.event...\n',size(allsac,1))
        EEG = addevents(EEG,allsac(:,[1 3 5 6 7 10 11:14]),{'latency','duration','sac_vmax','sac_amplitude','sac_angle','epoch','sac_startpos_x','sac_startpos_y','sac_endpos_x','sac_endpos_y'},'saccade');
    end
    
    % write fixations
    if writefix
        fprintf('\n--------------------------------------------------------------------')
        fprintf('\nAdding %i fixations to EEG.event...\n',size(allfix,1))
        if ldata && rdata % binocular
            EEG = addevents(EEG,allfix(:,[1 3 8 9 10]),{'latency','duration','fix_avgpos_x','fix_avgpos_y','epoch'},'fixation');
            %EEG = addevents(EEG,allfix(:,[1 3 4 5 6 7 10]),{'latency','duration','fixposition_ly','fixposition_lx','fixposition_ry','fixposition_rx','epoch'},'fixation');
        elseif ldata % left eye recorded
            EEG = addevents(EEG,allfix(:,[1 3 4 5 10]),{'latency','duration','fix_avgpos_x','fix_avgpos_y','epoch'},'fixation');        
        elseif rdata % right eye recorded
            EEG = addevents(EEG,allfix(:,[1 3 6 7 10]),{'latency','duration','fix_avgpos_x','fix_avgpos_y','epoch'},'fixation');
        end
    end
    fprintf('--------------------------------------------------------------------')
end
fprintf('\nDone.\n\n')
end