% addevents() - efficiently add new events to a EEGLAB dataset
%               adds new events to EEG.event, EEG.urvent (and EEG.epoch)
%               can be applied to continuous and epoched datasets
%
% Inputs:
%   newevents  - [matrix with n rows and c columns].
%                n = number of events to add
%                c = number of event properties (length(fieldnames))
%                The matrix must have lenght(fieldnames) columns. Each colum
%                contain the values of the corresponding fieldname.
%                Values for the following fields are mandatory:
%                - latency
%                - duration
%                - epoch (only in case of epoched datasets)
%   fieldnames - [cell of strings] names of subfields of
%                EEG.event/EEG.urevent. The following mandatory fieldnames
%                must be included in the list:
%                'latency'
%                'duration'
%                'epoch' (only in case of epoched datasets)
%   eventtype  - [string] type of event to be added. Note: Only one type of
%                event (EEG.event.type, e.g. "saccade") can be added with
%                each call of addevents(). However, many instances of this
%                particular event (e.g. 1000 saccades at once) can be added
%                with each call.
% Outputs:
%
%   EEG        - EEG datset with updated EEG.event & EEG.urevent structures
%
% See also: pop_detecteyemovements, detecteyemovements
%
% An example call might look like this:
% EEG = addevents(EEG,[1000 1 456 2; 2000 1 789 2;],{'latency','duration','someProperty','epoch'},'myEventType');
%
% This adds two new events of type "myEventType" to EEG.event and
% EEG.urevent at latencies 1000 and 2000. Events have a duration of 1 and
% both latencies fall into data epoch 2. The two new events will have the
% additional property 'EEG.event.someProperty' set to 456 and 789,
% respectively.
%
% Authors: ur & od
% Copyright (C) 2009-2017 Olaf Dimigen & Ulrich Reinacher, HU Berlin
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

function [EEG] = addevents(EEG,newevents,inputfldnames,eventtype)


% bugfix on 2015-03-27: ANT data importer does not fill the EEG.urevent structure
% leading to a crash of this function
if isempty(EEG.urevent)
    fprintf('\nWarning: EEG.urevent is empty. Recreating urevents with eeg_checkset before adding new events...\n')
    EEG = eeg_checkset(EEG,'makeur')
end

n_existingevents   = length(EEG.event);
n_existingurevents = length(EEG.urevent);
n_newevents        = size(newevents,1);

if isempty(EEG.event)
    existingfldnames = ''; % bugfix, Nov. 2017, OD, if EEG.event = [];
else
    % universal row vectors for name cells
    existingfldnames = fieldnames(EEG.event)';
end

% for compatibility with MATLAB 2010a and lower
[dim1 dim2] = size(inputfldnames);
if dim2==1 % iscolumn(inputfldnames)
    inputfldnames = inputfldnames';
    dim2 = dim1;
end
if dim2 ~= size(newevents,2)
    error('The number of columns with event data does not match the number of fieldnames.')
end
clear dim1 dim2
inputfldnames = [inputfldnames {'type'}];

%% get fields only present in EEG.event
[exfldonly,ib] = setdiff(existingfldnames, inputfldnames);

%% get fields only present in newly imported events
[newfldonly,ia] = setdiff(inputfldnames, existingfldnames);

%% make storage cell for EEG.events with ALL fields
storeOldEvent = cell(length(existingfldnames) + length(newfldonly), n_existingevents);

%% put EEG.event into storeOldEvent
storeOldEvent(1:length(existingfldnames),:) = struct2cell(EEG.event);

%% since new event properties are always numerical
storeOldEvent(length(existingfldnames)+1:end,:) = {0};

%% get type (numeric or string) of event fields in EEG.event
existing_isnumeric = all(cellfun(@isnumeric,storeOldEvent(ib,:)),2)';

%% create storage cell for new events with ALL fields
dummy2New = cell(n_newevents,length(exfldonly));
dummy2New(:, existing_isnumeric) = {0};
dummy2New(:,~existing_isnumeric) = {'none'};
storeNewEvent = [num2cell(newevents) repmat({eventtype},n_newevents,1) dummy2New];

%% build struct from cell
oldStruct = cell2struct(storeOldEvent,[existingfldnames newfldonly],1)';
newStruct = cell2struct(storeNewEvent,[inputfldnames exfldonly],2)';

%% put together old and new events
EEG.event        = [oldStruct newStruct];
existingfldnames = fieldnames(EEG.event)';


%% set pointers to urevents
for evnt = 1:size(newevents,1)
    newevent = n_existingevents + evnt;
    EEG.event(newevent).urevent = n_existingurevents + evnt;
end

%% add new urevents
% get list of non-mandatory subfields in EEG.event (e.g. 'code')
extrafields = setdiff(existingfldnames,[fieldnames(EEG.urevent);{'epoch'; 'urevent'}]);

% test whether subfields contain numeric or string info
extrafields_isnumeric2 = false(length(extrafields),1);
for fld = 1:length(extrafields)
    fldvalues2 = { EEG.event.(extrafields{fld}) };
    try
        %eval(['numeric = ~cellfun(''isclass'', fldvalues, ''char'');']);
        numeric2 = cellfun(@isnumeric, fldvalues2);
    catch %addeventsMyCellCatch
        numeric2 = mycellfun('isclass',fldvalues2,'double');
    end
    extrafields_isnumeric2(fld) = all(numeric2);
end
%%
for ue=1:n_existingurevents
    for fld = 1:length(extrafields)
        if extrafields_isnumeric2(fld)
            EEG.urevent(ue).(extrafields{fld}) = 0; % better: NaN?
        else
            EEG.urevent(ue).(extrafields{fld}) = 'none';
        end
    end
end
%%
dummy = EEG.event(n_existingevents+1:end);
dummy = rmfield(dummy,setxor(fieldnames(EEG.urevent),existingfldnames));
newurevent  = [EEG.urevent, dummy];
EEG.urevent = newurevent;

%% calculate ur-latencies of new events

%% scenario A: data is still continuous
if isempty(EEG.epoch)
    % true continuous or no consequence if EEG.xmin = 0
    latout = eeg_urlatency( EEG.event, [EEG.event(n_existingevents+1:end).latency] );
    if EEG.xmin ~= 0
        disp('deal with EEG.xmin')
        disp('pseudo-epoch')
        
        % development notes:
        % problems with special case: only one epoch, inconsistent format
        % in EEGLAB, i.e. "half continuous, half epoched" data
        % may happen if continuous data was cut into a single epoch
        % can there be epochation without the timelocking event belonging to
        % the epoch, e.g. epoch (+2ms:300ms) after stimulus :> messes with xmin
    end
    % go tru new events
    for ee = 1:n_newevents
        ee2 = n_existingevents + ee;
        EEG.urevent(EEG.event(ee2).urevent).latency = latout(ee);
    end
    
    %% scenario B: data already epoched
else
    offset = zeros(size(EEG.epoch));
    
    % retrieve infos for first existing event in epochs
    for ep = 1:length(EEG.epoch)
        ancor = EEG.epoch(ep).event(1); % index in EEG.event of 1st event in epoch
        % possible since epoch struct still mirrors old EEG.event. So there
        % is no chance for an event without correct urevent.latency.
        
        % get properties of 1st event in epoch...
        lat   = EEG.event(ancor).latency;  % latency in epoched (3D) dataset
        urev  = EEG.event(ancor).urevent;  % index in urevent
        urlat = EEG.urevent(urev).latency; % latency in org. cont. (2D) dataset
        % how is the latency of this 1st event in epoch after vs. before epochation?
        % = difference urlat-lat
        offset(ep) = urlat - lat;
    end
    
    % go tru new events
    for ee = 1:n_newevents
        ee2 = n_existingevents+ee;
        % BUGFIX: assign correct latencies for "new urevents" based on epoched dataset (Oct 18, 2017 by OD)
        EEG.urevent(EEG.event(ee2).urevent).latency = EEG.event(ee2).latency + offset(EEG.event(ee2).epoch);
        % The original, buggy line:
        %EEG.urevent(EEG.event(ee2).urevent).latency = EEG.event(ee).latency + offset(EEG.event(ee).epoch); % original code as before Oct-2017
    end
end

%% sort all of the new/merged events by latency
% do not shift this to anyplace earlier in the function.
% Calculation of urlatencies gets complicated otherwise
EEG = pop_editeventvals(EEG,'sort',{'latency' 0 }); % resort events
EEG = eeg_checkset(EEG,'eventconsistency');         % updates EEG.epochs (!)
end

%% helper function "mycellfun"
% The following function has been copied from function "eeg_checkset()"
% from the EEGLAB toolbox. The author is Arnaud Delorme.
% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
function res = mycellfun(com, vals, classtype)
res = zeros(1, length(vals));
switch com
    case 'isempty',
        for index = 1:length(vals), res(index) = isempty(vals{index}); end;
    case 'isclass'
        if strcmp(classtype, 'double')
            for index = 1:length(vals), res(index) = isnumeric(vals{index}); end;
        else
            error('unknown cellfun command')
        end;
    otherwise
        error('unknown cellfun command')
end
end