% geticvariance() - compute variance of independent component activations
%                   during saccade and fixation intervals 
%
% Usage:
%  >> [varsac,varfix] = geticvariance(EEG,sacstring,fixstring,tolerance)
%
% Inputs:
%  sacstring - [string] name of saccade event in EEG.event.type.
%              this is usually 'saccade', but may differ (e.g. 'L_saccade')
%              if saccades were imported from Eyelink raw data
%  fixstring - [string] name of fixation event in EEG.event.type.
%              this is usually 'fixation', but may differ if fixations were
%              imported from Eyelink raw data
%  tolerance - [vector with two integers] defines the tolerance around
%              saccade onset and saccade offset (in samples).
%              IC variance during saccade intervals will be computed from
%              saccade-onset minus tolerance(1) samples until saccade-offset
%              plus tolerance(2) samples. Conversely, tolerances are
%              subtracted from the fixation intervals. Two positive
%              integers are expected.
%
% Outputs:
%  varsac    - [vector of length ncomp] mean variance during saccade intervals
%
%  varfix    - [vector of length ncomp] mean variance during fixation intervals
%
%              output(:,1): mean IC variance during saccade intervals
%              output(:,2): mean IC variance during fixation intervals
%              output(:,3): ratio sac/fix variance, i.e.: output(:,1)./output(:,2)
%              Ratios above 1 indicate that the component is more active
%              during saccades than during fixations and thus likely
%              to reflect corneoretinal or myogenic ocular artifact
%
% See also: eyetrackerica, pop_eyetrackerica, geticvariance
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

function [varsac,varfix] = geticvariance(EEG,sacstring,fixstring,tolerance)

if length(tolerance) ~= 2
    error('%s(): tolerance needs to be a vector with two entries.')
end
tolerance = abs(tolerance);


%% get activation time courses of ICs
if isempty(EEG.icaact)
    icaact = eeg_getica(EEG);
else
    icaact = EEG.icaact;
end

if size(EEG.data,3) == 1 % data is continouos
    continuous = true;
else
    icaact = icaact(:,:); % convert to 2D
    continuous = false;
end
epochlength = size(EEG.data,2);


%% get latencies (onset sample) of saccade & fixation events
ix_sac = find(ismember({EEG.event.type},sacstring));
ix_fix = find(ismember({EEG.event.type},fixstring));

if ~isempty(ix_sac) && ~isempty(ix_fix)
    
    % preallocate
    varsac = NaN(size(icaact,1),length(ix_sac));
    varfix = NaN(size(icaact,1),length(ix_fix));
    
    % latencies and durations of EMs (in samples)
    saclat = [EEG.event(ix_sac).latency];
    sacdur = [EEG.event(ix_sac).duration];
    fixlat = [EEG.event(ix_fix).latency];
    fixdur = [EEG.event(ix_fix).duration];
    
    % go tru saccades
    fprintf('\n-- Calculating IC variance during %i saccades...',length(ix_sac));
    for ns = 1:length(ix_sac)
        % define samples belonging to saccade (with tolerance)
        lowr = saclat(ns) - tolerance(1);
        uppr = saclat(ns) + round(sacdur(ns))-1 + tolerance(2);
       
        if lowr <= 0, lowr = 1; end
        if uppr > size(icaact,2), uppr = size(icaact,2); end
        
        if ~continuous
            same_epoch = false;
            % epoched data: make sure data is from same epoch
            epoch_lowr = floor((lowr-1)/epochlength)+1; % epoch @ start of window
            epoch_uppr = floor((uppr-1)/epochlength)+1; % epoch @ end of window           
            if epoch_lowr == epoch_uppr
                same_epoch = true;
            end
        end
        
        if continuous || same_epoch
            % get icaact variance within this saccade
            varsac(:,ns) = var(icaact(:,lowr:uppr),0,2);
        end
    end
    fprintf(' computed for %i of %i saccade intervals', sum(~isnan(varsac(1,:))),ns);
    
    % go tru fixations
    fprintf('\n-- Calculating IC variance during %i fixations...',length(ix_fix));
    for nf = 1:length(ix_fix)
        % define samples belonging to fixation (with tolerance)
        lowr = fixlat(nf) + tolerance(2);
        uppr = fixlat(nf) + round(fixdur(nf))-1 - tolerance(1);
        
        % at least 2 samples left of fixation interval to compute variance?
        if lowr < uppr-1
            if ~continuous
                % data is epoched: make sure data is from same epoch
                same_epoch = false;
                epoch_lowr = floor((lowr-1)/epochlength)+1; % epoch @ start of window
                epoch_uppr = floor((uppr-1)/epochlength)+1; % epoch @ end of window
                
                if epoch_lowr == epoch_uppr
                    same_epoch = true;
                end
            end
            if continuous || same_epoch
                % get icaact variance within this saccade
                varfix(:,nf) = var(icaact(:,lowr:uppr),0,2);
            end
        else
            % catch extreme scenario: fixation duration is only 1 sample or
            % the tolerance around saccade is set so high that the 
            % remaining fixation interval is not even 2 samples long
            
            %fprintf('\nWarning: Remaining samples to compute fixation variance was less than 2 samples long. Skipped.')
        end
    end
    fprintf(' computed for %i of %i fixation intervals', sum(~isnan(varfix(1,:))),nf);
       
        
    if sum(~isnan(varsac(1,:))) ~= ns
        fprintf('\nNote: Variance was not computed for all saccade intervals, because intervals exceeded the epoch boundaries.');
    end
    if sum(~isnan(varfix(1,:))) ~= nf
        fprintf('\nNote: Variance was not computed for all fixation intervals, because intervals exceeded the epoch boundaries.\n');
    end
        
    % get mean sac and fix variances and the variance ratio (sac/fix)
    varsac    = nanmean(varsac,2);
    varfix    = nanmean(varfix,2);
else
    % error message: no saccade or fixation events of this name found in EEG.event.type
    if isempty(ix_sac)
        if isempty(ix_fix)
            error('%s(): There are neither saccade events nor fixation events called \"%s\" or \"%s\" in EEG.event.type!', mfilename,sacstring,fixstring)
        else
            error('%s(): Fixation events found, but found no saccade events called \"%s\" in EEG.event.type!', mfilename, sacstring)
        end
    else
        error('%s(): Saccade events found, but found no fixation events called \"%s\" in EEG.event.type!', mfilename, fixstring)
    end
end