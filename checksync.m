% checksync() - estimate quality of the data synchronization between
% eye-tracker and EEG by computing the crosscorrelation between the
% (horizontal) gaze position signal and the electro-oculogram (or other
% lateral frontal EEG channels that receive strong corneoretinal artifacts). 
% Since both signals reflect the rotation of the eye ball, there should be
% virtually no time lag between both signals, that is, the cross
% correlation should peak at a lag of zero. This function computes the
% crosscorrelation between gaze and EOG, plots it, and stores the results
% in the field "EEG.misc.xcorr_EEG_ET" of the EEG structure.
%
% Usage:
%   >> EEG = pop_checksync(EEG,eye_x,heog_l,heog_r,plotfig)
%
% Required inputs:
%   EEG          - [string] EEG struct, also containing synchronized eye 
%                  tracking data (see pop_importeyetracker)
%   gaze_x        - [one channel index], 
%                  specify channel index of the X-component (horizontal)
%                  gaze position signal of the eye tracker. If you recorded
%                  from both eyes (binocular), you can also enter a *vector* 
%                  with two entries (x of left and right eye)
%   heog_l       - [channel index], index of EOG channel on LEFT side of 
%                  the head/face. If no proper EOG electrode was placed
%                  near the eye (bad idea!), use a lateral EEG channel that 
%                  is as close to the eye as possible
%   heog_r       - [channel index], index of EOG channel on RIGHT side of 
%                  the head/face. If no proper EOG electrode was placed
%                  near the eye (bad idea!), use a lateral EEG channel that 
%                  is as close to the eye as possible
%  plotfig       - show a plot with result of the crosscorrelation function
%
% Outputs:
%   EEG         - EEG structure with cross-correlation info added to the
%                 field EEG.misc.xcorr_EEG_ET
%
% See also: pop_checksync, pop_importeyetracker, synchronize
%
%
% An example call of the function might look like this: 
% >> EEG = pop_checksync(EEG,65,1,2,1)
%
% In this example, the horizontal gaze position of the left eye from the 
% eye tracker is stored in channel 65 of the synchronized EEG/ET dataset. 
% The signal of the horizontal EOG electrodes, placed on the left and right 
% canthus of each eye, is stored in channels 1 and 2, respectively. 
% The last input indicates that EYE-EEG should plot a figure showing the 
% cross-correlation function between the ET channels and the, in this case, 
% bipolar-referenced (left minus right EOG electrode) EOG. 
% If the synchronization is good, the time lag of the maximum (peak) of the
% cross-correlation should be near zero.
% 
% IMPORTANT NOTE: To get clean results, you should first remove
% bad or missing data from the eye tracking channels, e.g. due to
% eye blinks, using method pop_rej_eyecontin(). Otherwise, the cross-corr.
% function will be distorted by this missing data. In some cases, it might 
% also be necessary to remove intervals with very bad EOG recordings, but 
% in our experience, some bad EOG data is not much of a problem.
%
% Note that this method was first applied in Dimigen et al., 2011, JEP:GEN. 
% Please cite this paper when using this method. Thanks!
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


function [EEG] = checksync(EEG, chan_gaze_x, chan_heog_l, chan_heog_r, plotfig)

fprintf('\n%s(): Computing cross correlation function between eye tracker and EOG channels...',mfilename);
fprintf('\n\nNote: you need clean eye-tracking and EOG data for this function to produce sensible results.');
fprintf('\nFor example, it is usually necessary to first reject intervals with bad or missing ET data (e.g. blinks)');
fprintf('\nbefore computing the cross-correlation function with the EOG. Otherwise the results will be distorted.');
fprintf(sprintf('\nYou can do this by clicking: \"Reject data based on eye tracker\" > \"Reject bad contin. data\"'));

MAXLAG = 100; % compute cross-corr. up to max. lag of MAXLAG samples

%% warning if data were already epoched
if size(EEG.data,3)>1
    warning(sprintf('This function was designed to run on continuous (not epoched) data.\nIt may produce unreliable results if the time series are too short!'))
end

%% get possible "bad_ET" intervals from EEG.event structure
badvector = zeros(1,size(EEG.data,2)*size(EEG.data,3));
ix_badETevent = find(ismember({EEG.event.type},'bad_ET'));
if ~isempty(ix_badETevent)
    fprintf('\nFound \"bad_ET\" events in EEG.events...')
    fprintf('\nThese bad eye-tracker intervals will be ignored when computing the cross-correlation.')
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


%% get data for cross-correlation
gaze_x = mean(EEG.data(chan_gaze_x,~badvector),1); % "mean", because user can also input several channels

if chan_heog_l == chan_heog_r
    fprintf('\n\n%s(): Same channel number was provided for the left and right HEOG electrode!',mfilename);
    fprintf('\nTherefore, it is assumed that either (1) the HEOG was already recorded in a bipolar montage as one channel (e.g., L vs. R)');
    fprintf('\nor (2) that only one EOG electrode was placed near the eyes (e.g. only at left but not at right canthus).');
    fprintf('\nCross-correlation function will therefore be based on this one EOG channel.');
    heog = EEG.data(chan_heog_l,~badvector);
else
    fprintf('\n%s(): Computing one bipolar EOG channel from both EOG electrodes...',mfilename);
    % subtract right minus left EOG channels to obtain *positive* correlations with hor. ET
    heog   = EEG.data(chan_heog_r,~badvector)-EEG.data(chan_heog_l,~badvector);
end

%% compute cross correlation
% based on entire recording (may include artifacts & loss of tracking!)

% it seems like xcorr is part of the signal processing toolbox
% workaround (from https://stackoverflow.com/questions/7396814/cross-correlation-in-matlab-without-using-the-inbuilt-function)
% check whether signal processing toolbox installed
v=ver; [installedToolboxes{1:length(v)}] = deal(v.Name);

% if license('test','Signal_Toolbox') % alternative 
if all(ismember('Signal Processing Toolbox',installedToolboxes))  % not sure about name
    [xc,lags] = xcorr(gaze_x,heog,MAXLAG); % get cross-correlation
else
    corrLength = MAXLAG*2-1;
    %corrLength=length(a)+length(b)-1;
    lags = [MAXLAG-corrLength:corrLength-MAXLAG];
    xc=fftshift(ifft(fft(gaze_x,corrLength).*conj(fft(heog,corrLength))));
end

[maxvalue,ix] = max(abs(xc)); % find maximum xc
sampleDiff = lags(ix);  % find lag with maximum xc

%% user feedback
if sampleDiff < 0
    fprintf('\n\n-- Maximum cross-correlation is observed at lag of %i samples (= %.2f ms):',sampleDiff,sampleDiff*(1000/EEG.srate));
    fprintf('\n-- The eye tracker signal leads the EOG signal');
elseif sampleDiff > 0
    fprintf('\n\n-- Maximum cross-correlation is observed at lag of %i samples (= %.2f ms):',sampleDiff,sampleDiff*(1000/EEG.srate));
    fprintf('\n-- The eye tracker signal lags behind the EOG signal');
else
    fprintf('\n\n-- Maximum cross-correlation is observed at lag of %i samples (= %.2f ms):',sampleDiff,sampleDiff*(1000/EEG.srate));
    fprintf('\n-- Gaze and EOG seem perfectly aligned');
end

%% also report correlation coefficient (r) of the two full time series (at lag 0)
r = corrcoef(gaze_x,heog);
fprintf('\n-- At a lag of zero, the Pearson correlation coefficient between both signals is r = %.2f',r(1,2));

%% write values to "EEG.etc" so they are stored with the dataset
fprintf('\n\nThe cross-correlation function is stored in the field \"EEG.etc.xcorr_eyeeeg\".');
EEG.etc.xcorr_eyeeeg = [xc; lags]; % first row: xc-values, second row: lags

%% show figure with cross-correlation function
if plotfig
    fprintf('\n%s(): Plotting cross-correlation function... ',mfilename);
    
    figure('Name','synccheck: ET/EOG cross-correlation');
    %% plot cross-correlations between -MAXLAGS and +MAXLAGS
    subplot(1,2,1); hold on;
    title({'Cross-correlation Eye-tracker & EOG';'\fontsize{8}'});
    plot(lags,xc,'k.-','markersize',4) % plot xc function
    % highlight lag zero and lag with maximum absolute cross-correlation
    ylimits = ylim;
    %ylimits(2) = ylimits(2)+0.1*ylimits(2);
    h1 = plot([0 0],[ylimits(1) ylimits(2)],'color',[0.6 0.6 0.6],'linewidth',2.0,'linestyle',':'); % lag zero
    h2 = plot([lags(ix) lags(ix)],[ylimits(1) xc(ix)],'r','linewidth',1.2); % lag with max xc
    % xlabel, also explain what pos./neg. lag means
    l = legend([h1 h2],{'lag zero','max. absolute cross-corr.'},'Location','Best','box','on');
    %xlabel({'\fontsize{10}Lag (in samples)',['\fontsize{8}ET leads EOG '
    %char(8592) ' 0 ' char(8594) ' ET lags  EOG']}); % arrows do not work in R2012
    xlabel({'\fontsize{10}Lag (in samples)',['\fontsize{8}ET leads EOG <--  0  -->  ET lags  EOG']});
    ylabel('Cross-correlation value')
    ylim(ylimits)
    axis square, box on
    
    %% zoom in: plot cross-correlation between lag -5 and 5 only
    subplot(1,2,2); hold on;
    %title(sprintf('Cross-correlation ET and EOG (zoom)'),'fontweight','bold')
    title({'Cross-correlation Eye-tracker & EOG';'\fontsize{8}(zoom)'});
    xlim([-5 5]);
    plot(lags,xc,'k.-','markersize',6) % plot xc function
    % highlight lag zero and lag with maximum absolute cross-correlation
    ylimits = ylim;
    h3 = plot([0 0],              [ylimits(1) ylimits(2)],'color',[0.6 0.6 0.6],'linewidth',2.0,'linestyle',':'); % lag zero
    h4 = plot([lags(ix) lags(ix)],[ylimits(1) xc(ix)],'r','linewidth',1.2); % lag with max xc
    % legend
    %xlabel({'\fontsize{10}Lag (in samples)',['\fontsize{8}ET leads EOG ' char(8592) ' 0 ' char(8594) ' ET lags  EOG']});
    xlabel({'\fontsize{10}Lag (in samples)',['\fontsize{8}ET leads EOG <--  0  -->  ET lags  EOG']});
    ylabel('Cross-correlation value')
    ylim(ylimits)
    axis square, box on
    fprintf('Done.\n');    
end