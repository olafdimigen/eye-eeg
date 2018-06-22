% parseeyelink() - parse Eyelink(TM) eye tracking data from an
% ASCII-converted text-file and save the content as a MATLAB structure
%
% Usage:
%   >>  et = parseeyelink(inputFile, outputMatFile, triggerKeyword);
%
% Inputs:
%   inputFile      - text (ascii) file, converted from raw ET recording
%   outputMatFile    - MATLAB (.mat) file to save "et" structure in
%
% Optional Inputs:
%   triggerKeyword - keyword contained in special messages used for
%                    synchronization. Text messages can be used as an
%                    alternative method for synchronization with the
%                    EEG/MEG (instead of shared trigger pulses). However,
%                    all messages used for synchronization must contain a
%                    special keyword string (e.g. "MYKEYWORD") which is
%                    followed by the respective integer number (e.g. "123")
%                    which is also send as a trigger pulseto the parallel
%                    electrophysiological recording.
%
% For more information on the optional input, please read the tutorial of
% the toolbox.
%
%   An example use of the function might look like this:
%   >> ET = parseeyelink('/myrawdata/SAMPLES.txt','/myfiles/SAMPLES.mat')
%   Optional with keyword:
%   >> ET = parseeyelink('/myrawdata/SAMPLES.txt','/myfiles/SAMPLES.mat','MYKEYWORD')
%
% Outputs:
%   et     - eyetracker structure
%
% See also:
%   pop_parsesmi, pop_importeyetracker
%
% author: ur
% Copyright (C) 2009-2018 Olaf Dimigen & Ulrich Reinacher, HU Berlin
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

function et = parseeyelink(inputFile,outputMatFile,triggerKeyword)

if nargin < 2
    help(mfilename);
    return;
end;

%% feedback
fprintf('\n%s(): Parsing Eyelink raw data. This may take a moment...',mfilename)
fprintf('\nNote: missing data samples (e.g. due to eye blinks) will be stored as zeros in parsed data.\n')
fprintf('\n-- reading txt...')

fid = fopen(inputFile,'r');
if fid == -1
    fprintf('\nThe source file ''%s'' was not found.\n',inputFile)
    return
    % throw specific error as caller?
end

B = fread(fid,'*char')';
fclose(fid);
fprintf('\n\tDone.')

%% file specifications and comments
fprintf('\n-- getting comments and column headers...')

et.comments = regexp(B,'(^[*#;/]{2}[^\r\n]*\r?\n)','match','lineanchors');

%% build column header
triggersAsColumn = false; % flag for INPUT column in data
triggersAsMessages = false;  % flag for INPUT events

edfColHeader = regexp(B,'^SAMPLES\s(GAZE|HREF|RAW)\s([^\r\n]*)','tokens','once','lineanchors');

if ~isempty(edfColHeader)
    % are all possible keywords in the matchlist?
    edfColType = regexp(edfColHeader{2},{'LEFT','RIGHT','VEL','RES','INPUT','HTARGET'}, 'once');
    
    % possible to do this with regexprrep?
    header = cell(1,17);
    header{1} = 'TIME';
    if edfColType{1}
        header{2} = sprintf('L_%s_X',edfColHeader{1});
        header{3} = sprintf('L_%s_Y',edfColHeader{1});
        header{4} = 'L_AREA';
        if edfColType{3}
            header{8} = 'L_VEL_X';
            header{9} = 'L_VEL_Y';
        end
    end
    if edfColType{2}
        header{5} = sprintf('R_%s_X',edfColHeader{1});
        header{6} = sprintf('R_%s_Y',edfColHeader{1});
        header{7} = 'R_AREA';
        if edfColType{3}
            header{10} = 'R_VEL_X';
            header{11} = 'R_VEL_Y';
        end
    end
    if edfColType{4}
        header{12} = 'RES_X';
        header{13} = 'RES_Y';
    end
    if edfColType{5}
        header{14} = 'INPUT';
        triggersAsColumn = true;
    end
    
    if edfColType{6}
        header{15} = 'targetpos_x';
        header{16} = 'targetpos_y';
        header{17} = 'targetdistance';
    end
    
    et.colheader = header(~cellfun(@isempty,header));
else
    et.colheader = {};
end
clear edfCol* header
fprintf('\n\tDone.')

%% get data
fprintf('\n-- getting data...')
dataLines =  regexp(B,'^\d[^\r\n]*','match','lineanchors'); % bug 2015-06-27a

% bug 2015-06-27a
% bugreport june, 27, 2015: new Eyelink lines reporting calibration model
% current parser cannot distinguish them from regular data lines starting with a timestamp
% e.g.:
% MSG	11797285 !CAL Cal coeff:(X=a+bx+cy+dxx+eyy,Y=f+gx+goaly+ixx+jyy)
% 32407  491.81  902.09  3.5305  12.203
% 11033  10.253  294.18 -0.2066  2.6616

% replace default empty value token '.' by '0'. (used for missing data)
dataLines = regexprep(dataLines,'(\s)\.(\s)','$10$2');

% remove corneal reflex mode quality data. Also removes superfluous white
% space. Need to keep one whitespace at the end, since sscanf will
% otherwise merge numbers over linebreaks
dataLines = regexprep(dataLines,'[\D]{2,}',' ');
dataChar = char(dataLines);
dataChar(:,end+1)= char(32);

% test once for number of columns
nColData = length(sscanf(dataChar(1,:),'%f'));

% handle weird behaviour in monoc recording that introduces single line of
% binoc zeroes at end of a block
nColDatAll = arrayfun(@(x) length(sscanf(dataChar(x,:),'%f')), 1:length(dataChar));
dataChar(nColDatAll~=nColData,:)=[];

% sscans string array by !column! for float-numbers, arranges in rows of
% nColumns, till end of string (inf). transpose back to column order.
et.data = sscanf(dataChar','%f',[nColData,inf])';

clear dataLines dataChar
fprintf('\n\tDone.')

%% check for empty columns in data array
nColh = length(et.colheader);

if nColh ~= nColData
    warning('\n%s(): The number of column names in the header (%i) differs from the number of actual data columns (%i)',mfilename,nColData,nColh)
end

%% get messages and events
fprintf('\n-- getting messages...')

% % !! Caused problems under OSX: eternal matching, HOW?
% % necessary since some messages get split up over multiple lines (beginning
% % with whitespace characters)
% B = regexprep(B,'\r?\n[\t\f ]+','\t');
% use following line to get
et.messages = regexp(B,'^[A-Z \t][^\r\n]*\r?\n','match','lineanchors');
%et.messages = regexp(B,'^[A-Z][^\r\n]*\r?\n','match','lineanchors');
if regexp(B,'^INPUT','once','lineanchors')
    triggersAsMessages = true;
end
fprintf('\n\tDone.')

%% online detection: get recording type for events
eventDataType = regexp(B,'^EVENTS\s(GAZE|HREF|RAW)','tokens','once','lineanchors');
clear B*
%% online detection: get saccades detected by Eyelink
saccadePattern = '^ESACC\s+([RL])\s+([^\r\n]*)';
tokensSaccade = regexp(et.messages,saccadePattern,'tokens','once')';

tokensSaccade = vertcat(tokensSaccade{:});
if ~isempty(tokensSaccade)
    % rare case of empty data in eyelink events
    tokensSaccade(:,2) =  regexprep(tokensSaccade(:,2),'(\s)\.(\s)','$10$2');
    
    saccades.eye = [tokensSaccade{:,1}]';
    saccades.data = cell2mat(cellfun(@str2num,tokensSaccade(:,2),'UniformOutput',false));
    
    saccades.colheader = {'latency','endtime','duration','sac_startpos_x','sac_startpos_y','sac_endpos_x','sac_endpos_y','sac_amplitude','sac_vmax','resolution_x','resolution_y'};
    saccades.colheader = saccades.colheader(1:size(saccades.data,2));
    et.eyeevent.saccades = saccades;
end

%% online detection: get fixations detected by Eyelink
fixationPattern = '^EFIX\s+([RL])\s+([^\r\n]*)';
tokensFixation = regexp(et.messages,fixationPattern,'tokens','once')';

tokensFixation = vertcat(tokensFixation{:});
if ~isempty(tokensFixation)
    % rare case of empty data in eyelink events
    tokensFixation(:,2) = regexprep(tokensFixation(:,2),'(\s)\.(\s)','$10$2');
    
    fixations.eye = [tokensFixation{:,1}]';
    fixations.data = cell2mat(cellfun(@str2num,tokensFixation(:,2),'UniformOutput',false));
    switch char(eventDataType)
        case {'GAZE','RAW'}
            fixations.colheader = {'latency','endtime','duration','fix_avgpos_x','fix_avgpos_y','fix_avgpupilsize','resolution_x','resolution_y'};
        case 'HREF' % "exotic" configuration
            fixations.colheader = {'latency','endtime','duration','href_x','href_y','avgpupilsize','displaypos_x','displaypos_y','avgpupilsize','resolution_x','resolution_y'};
    end
    fixations.colheader = fixations.colheader(1:size(fixations.data,2));
    et.eyeevent.fixations = fixations;
end

%% online detection: get eye blinks detected by Eyelink
blinkPattern = '^EBLINK\s+([RL])\s+([^\r\n]*)';
tokensBlink = regexp(et.messages,blinkPattern,'tokens','once')';

tokensBlink = vertcat(tokensBlink{:});
if ~isempty(tokensBlink)
    blinks.eye = [tokensBlink{:,1}]';
    blinks.data = str2num(char(tokensBlink{:,2}));
    blinks.colheader = {'latency','endtime','duration'};
    et.eyeevent.blinks = blinks;
end

%% clean up
clear tokens* saccades fixations blinks *Pattern


%% build table of events (parallel port inputs or messages with keyword)

% 1. check for existence of messages containing special keyword
if exist('triggerKeyword','var')
    
    if isempty(triggerKeyword)
        warning('\n\n%s(): You are using an empty string as the keyword.',mfilename)
    end
    
    % handle leading whitespaces as well as missing trailing whitespaces
    % (handle "MYKEYWORD 123" as well as "MYKEYWORD123")
    test  = regexp(et.messages,['MSG\s(\d+)\s+',triggerKeyword,'\s?(\d+)'],'tokens')';
    test2 = [test{:}]';
    test3 = cellfun(@str2double,test2,'UniformOutput',false);
    et.event = cell2mat(test3);
    if isempty(et.event)
        warning('\n\n%s(): The keyword string you have provided (''%s'') was not found in any messages of the ET file (%s). Please check your messages and make sure they contain the specified keyword.',...
            mfilename, triggerKeyword, inputFile)
    end
    
    % 2. check for parallel port inputs (stored in messages called "INPUT")
elseif triggersAsMessages
    test  = regexp(et.messages,'INPUT\s+(\d+)\s+(.*)','tokens')';
    test2 = [test{:}]';
    test3 = cellfun(@str2double,test2,'UniformOutput',false);
    et.event = cleantriggers(cell2mat(test3));
    
    % 3. check for parallel port inputs (stored in dedicated raw data column)
elseif triggersAsColumn
    et.event = cleantriggers( [et.data(:,1),et.data(:,strcmp('INPUT',et.colheader))]);
    
    
else % no events found: user feedback
    warning('\n%s(): Did not find trigger pulses or special messages that could be used as events for synchronization.',mfilename)
    if ~exist('triggerKeyword','var')
        fprintf('\nIf you have sent messages with a special keyword, please provide the keyword when calling this function. See help.')
    end
end

%% Extract timestamps of *other* messages (not eyeevent, not keyword-messages)
% New feature, Jan-2018, OD
% Identify messages that begin with 'MSG' (e.g. user-send messages) so that 
% their timestamps can be later converted into "eeg time" in function 
% synchronize(). Extracting this information should be helpful for user 
% who have coded their experimental conditions (etc.) in the form of 
% (Eyelink) messages

% go tru messages in et.messages
remove_this_msg = zeros(size(et.messages,2),1);
contains_keyword = [];
for m = 1:size(et.messages,2)
    % does message line start with 'MSG'? --> include
    contains_msg = strcmp(et.messages{m}(1:3),'MSG');
    % does message contain the user-defined keyword (sync message) --> exclude
    if exist('triggerKeyword','var')
        contains_keyword = ~isempty(strfind(et.messages{m},triggerKeyword));
    end
    % kick out
    if ~contains_msg | contains_keyword
        remove_this_msg(m) = true;
    end
end
% now delete messages not starting with 'MSG' (e.g. eyeevents) or containing
% keyword since they are already stored elsewhere in the output. Also, this
% removes a few strange line-broken message lines that contain only numbers
% or only alphabetic chars (related to the calibration model)
tmp = et.messages;
tmp(:,find(remove_this_msg)') = []; % remove msg marked above

% now extract numeric timestamp from these remaining msg
timestamps = regexp(tmp,'^-?\D+(\d+)','tokens','once','lineanchors');
timestamps = [timestamps{:}]';
timestamps = cellfun(@str2double,timestamps,'UniformOutput',false);
timestamps = cell2mat(timestamps);

% store these messages in structure 'et.othermessages'
for j = 1:size(tmp,2)
    et.othermessages(j).timestamp = timestamps(j); % numeric timestamps
    et.othermessages(j).msg       = tmp{j}; % the msg string
end
clear tmp test* edfColHeader

%% save
fprintf('\n-- saving structure as mat...')
save(outputMatFile,'-struct','et')
fprintf('\n\tDone.')

fprintf('\n\nFinished parsing eye track.\nSaved .mat file to %s.\n',outputMatFile)