% eyetrackerica() - reject independent components based on their 
%                   covariance with the simultaneously recorded eye track
%
% Usage:
%   >> EEG = eyetrackerica(EEG,sacstring,fixstring,sactolerance,...
%                          threshratio,flagmode,plotfig,topomode)
%
% Inputs:
%  EEG        - [string] EEG struct. EEG must already contain extra 
%               channels with synchronized eye tracking data
%               (>> pop_importeyetracker()) and both saccade and fixation
%               events (imported or by using pop_detecteyemovements)
%  sacstring  - [string] name of saccade event (in EEG.event).
%               This is 'saccade' per default, but may differ if saccades 
%               were imported from Eyelink raw data (e.g. 'L_saccade')
%  fixstring  - [string] name of fixation event (in EEG.event).
%               This is 'fixation' per default, but may differ if fixations 
%               were imported from Eyelink raw data (e.g. 'R_fixation')
%  sactol     - [vector of two integers], e.g. [5 10]. 
%               Extra temporal tolerance around saccade onset and offset.
%               Saccade intervals will range from saccade-onset minus 
%               window(1) samples until saccade-offset plus window(2) 
%               samples. Correspondingly, the tolerance will be substracted
%               from the fixation intervals.
%  threshratio - critical variance ratio threshold for component rejection.
%               If the variance of IC activity within saccades is more than
%               threshratio larger than the IC activity during fixations, the
%               component will be flagged for rejection.
%  flagmode  -  [integer: either 1,2, or 3]
%               1: do not flag any components for rejection (test mode)
%               2: flag components for rejection in EEG.reject.gcompreject, 
%               keep existing rejection flags (e.g. those flagged manually)
%               3: flag components for rejection in EEG.reject.gcompreject,
%               overwrite all existing flags in EEG.reject.gcompreject, i.e.
%               components formely flagged as "bad" may be changed to "good"
%  plotfig   -  <boolean> [0/1] - show figure to evaluate the rejection
%               threshold
%  topomode  -  <integer: 1,2,3, or 4> - determines whether after 
%               eyetrackerica a figure is shown with the 2D topographies of 
%               flagged (bad and/or good) components (using pop_selectcomp). 
%               Figure allows to toggle rejection flags manually.
%               Use one of four options:
%               1: plot topos of "bad" and "good" components (2 figures)
%               2: plot topos of "bad" components only (1 figure)
%               3: plot topos of "good" components only (1 figure)
%               4: do not plot topographies
%               Note: Option 1:3 requires electrode coordinates.
%
% Outputs:
%   EEG       - EEG structure with updated field EEG.reject.gcompreject. 
%               If rejectmode 1 or 2 was used, EEG.reject.gcompreject will 
%               contain updated flags flagging components with a variance
%               ratio larger than threshratio 
%   vartable  - [matrix of size ncomp x 3]
%               Columns of [vartable] contain the following information:
%               vartable(:,1): mean IC variance during saccade intervals
%               vartable(:,2): mean IC variance during fixation intervals
%               vartable(:,3): ratio, i.e., vartable(:,1)/vartable(:,2)
%               High ratios indicate that the component is more active
%               during saccades than during fixations, and may therefore
%               reflect corneoretinal or myogenic ocular artifact
%
% To actually remove the flagged components from your data, use EEGLAB 
% menu "Tools" > "Remove Components" or function pop_subcomp() afterwards
%
% See also: pop_eyetrackerica, geticvariance
%
% Example: A call of the function might look like this:
%
% >> [EEG vartable] = eyetrackerica(EEG,'saccade','fixation',[10 5],1.1,3,1,1)
%
% Explanation: compute the variance of the activity time course of each IC
% within saccade events (of type 'saccade' in EEG.event) and fixation
% events (of type 'fixation' in EEG.event). For this analysis, saccade 
% intervals are defined as lasting from 10 samples before fixation onset 
% until 5 samples after saccade offset [sactol]. Conversely, variance for
% fixation intervals is computed from 5 samples after fixation onset until 
% 10 samples before fixation offset. If the variance ratio 
% [var(sac)/var(fix)] for a component is larger than 1.1, set this 
% component to "reject". Overwrite existing rejection flags already present
% in EEG.reject.gcompreject. Show an extra figure to evaluate the setting
% for the threshold (maxcrit). Plot 2D topographies of components flagged
% as good and bad.
%
% The saccade/fixation "variance ratio" criterion was proposed by:
%
% Plöchl, M., Ossandon, J.P., & König, P. (2012). Combining EEG and 
% eye tracking: identification, characterization, and correction of eye 
% movement artifacts in electroencephalographic data. Frontiers in Human 
% Neuroscience, doi: 10.3389/fnhum.2012.00278.
%
% Author: od
% Copyright (C) 2009-2013 Olaf Dimigen & Ulrich Reinacher, HU Berlin
% olaf.dimigen@hu-berlin.de / ulrich.reinacher.1@hu-berlin.de

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

function [EEG vartable] = eyetrackerica(EEG,sacstring,fixstring,sactolerance,threshratio,flagmode,plotfig,topomode)

%% check whether ICA was computed
if isempty(EEG.icaweights)
    fprintf('%s(): Error: no ICA decomposition. Use menu "Tools > Run ICA" first.',mfilename);
    return
end

%% check whether there are saccade & fixation events of specified type
hasSac = any(ismember({EEG.event.type},sacstring));
hasFix = any(ismember({EEG.event.type},fixstring));
if ~hasSac && ~hasFix
    error('%s(): Saccade string \"%s\" and fixation \"%s\" were both not found in EEG.event. Did you detect eye movements?',mfilename,sacstring,fixstring);
elseif ~hasSac
    error('%s(): Saccade string \"%s\" was not found in EEG.event. Did you already detect saccades?',mfilename,sacstring);
elseif ~hasFix
    error('%s(): Fixation string \"%s\" was not found in EEG.event. Did you already detect fixations?',mfilename,fixstring);
end

%% get variance ratio for all ICs
fprintf('\n\n%s(): Running eye-tracker informed component selection...\n',mfilename)
[varsac, varfix] = geticvariance(EEG,sacstring,fixstring,sactolerance);
varratio = varsac./varfix;
vartable = [varsac varfix varratio];


%% feedback: MATLAB console
fprintf('\n\nIC   var(sac)   var(fix)   ratio')
fprintf('\n--------------------------------')
for n = 1:length(varratio)
    fprintf('\n%i\t%.2f\t%.2f\t%.2f',n,varsac(n),varfix(n),varratio(n));
    if varratio(n) > threshratio
        fprintf(' *')
    end
    if ~rem(n,10)
        fprintf('\n--------------------------------')
    end
end

% remember existing rejection flags in EEG.reject.gcompreject
if isfield(EEG,'reject')
    old_badcomps = EEG.reject.gcompreject;
else
    old_badcomps = [];
end

%% set rejection flag for components that exceed variance ratio threshold
new_badcomps = varratio > threshratio;
new_badcomps = new_badcomps';

fprintf('\n------------------------------\n')
fprintf('\n%s(): %i of %i components exceed the variance ratio of %.3f',mfilename,sum(new_badcomps),length(new_badcomps),threshratio);

%% how to handle super-threshold components?
switch flagmode
    case 1
        % test mode, do not flags any components as "bad"
        fprintf('\n%s(): test mode. No components were flagged for rejection.',mfilename);
        
    case 2
        % add flags if there are "bad" components that are not yet flagged
        EEG.reject.gcompreject = new_badcomps | old_badcomps;
        fprintf('\n%s(): %i previously \"good\" components were flagged for rejection.',mfilename,sum(new_badcomps)-sum(old_badcomps));
        fprintf('\n%s(): Now a total of %i components is flagged for rejection.',mfilename,sum(EEG.reject.gcompreject));
        
    case 3 % overwrite existing flags with current flags
        EEG.reject.gcompreject = new_badcomps;
        fprintf('\n%s(): %i of %i components are now flagged for rejection.',mfilename,sum(new_badcomps),length(old_badcomps));
        if sum(old_badcomps)>0
            fprintf('\n%s(): Existing flags in EEG.reject.gcompreject were overwritten.',mfilename);
        end
end

if ismember(flagmode,2:3)
    fprintf('\n%s(): Use \"Tools > Remove components\" or pop_subcomp() to remove flagged components',mfilename);
end


%% show figure to evaluate variance ratio threshold
if plotfig
    fprintf('\n%s(): Plotting figure with eyetrackerica results...',mfilename);   
    figure('Name','eyetrackerica: Results');
    ncomp = length(varratio);
    
    %% ICA variance ratios vs. threshold (threshold)
    subplot(1,2,1); hold on;
    title('IC variance ratios','fontweight','bold')
    fill([0.5 ncomp+0.5 ncomp+0.5 0.5],[0 0 1 1],[.85 .85 .85],'EdgeColor','none'); % illustrate ratio of "1"
    bar(1:ncomp,varratio,'k');
    dummy = zeros(ncomp,1);
    ixbadcomp = find(varratio > threshratio);
    dummy(ixbadcomp) = varratio(ixbadcomp);
    % highlight components exceeding threshold
    bar(1:ncomp,dummy,'r');
    %h1 = plot([0 ncomp+1],[1 1],'Color',[.4 .4 .4]); % line: ratio of "1"
    h2 = plot([0 ncomp],[threshratio threshratio],'r:','linewidth',1.2); % line: threshold
    xlabel('Independent component #')
    ylabel('Variance ratio [sac/fix]')
    xlim([0.5 ncomp+0.5])
    l = legend(h2,sprintf('Threshold of %.3f',threshratio));
    set(l,'box','off','Fontsize',8);
    %axis square; 
    box on;
    
    %% plot threshold (threshratio) vs. # of bad components
    ratios = 0:0.1:3;
    for r = 1:length(ratios)
        n_rejected(r) = sum(varratio > ratios(r));
    end
    subplot(1,2,2); hold on;
    title('Effect of alternative thresholds','fontweight','bold')
    plot(ratios,n_rejected,'k.-','linewidth',1.0);
    ylim([0 ncomp+1]);
    plot([threshratio threshratio],[0 length(ixbadcomp)],'r:','linewidth',1.2)
    plot([0 threshratio],[length(ixbadcomp) length(ixbadcomp)],'r:','linewidth',1.2)
    plot(threshratio,length(ixbadcomp),'ro','linewidth',1)
    set(gca,'xdir','reverse');
    ylabel('No. of rejected components')
    xlabel('Threshold setting')
    %axis square;
    box on;
end


%% plot topos of flagged components
if flagmode == 1 && ismember(topomode,1:3)
    fprintf('\n%s(): Test mode. Rejection status of components was not changed. Not plotting any topographies.',mfilename);
else
    
    %% plot bad+good/bad/good/nothing?
    plotbad = false; plotgood = false;
    switch topomode
        case 1
            plotbad  = true; plotgood = true;
        case 2
            plotbad  = true;
        case 3
            plotgood = true;
        case 4
            fprintf('\nDone.\n')
            return
        otherwise
            fprintf('\n%s(): Input for topomode not recognized. Not plotting any topographies.\n',mfilename)
            return
    end
    
    %% test for channel locations
    if isempty(EEG.chanlocs) || ~isfield(EEG.chanlocs,'theta') || all(cellfun('isempty',{EEG.chanlocs.theta}))
        fprintf('\n%s(): Cannot plot topographies without channel location information.\n',mfilename);
        EEG = eeg_checkset(EEG,'chanloc');
    else
        % plot bad components
        if plotbad
            bad  =  find(EEG.reject.gcompreject);
            if ~isempty(bad)
                fprintf('\n%s(): Plotting topos of \"bad\" components flagged for rejection.\n',mfilename);
                pop_selectcomps(EEG,bad);
            else
            end
        end
        % plot good components
        if plotgood
            good =  find(~EEG.reject.gcompreject);
            if ~isempty(good)
                fprintf('\n%s(): Plotting topos of \"good\" components NOT flagged for rejection.\n',mfilename);
                pop_selectcomps(EEG,good);
            else
            end
        end
    end
end
fprintf('\nDone.\n')
end