% parsesmi() - parse Sensomotoric Instruments (TM) eye tracking
% data from an ASCII-converted text-file and save the content as a MATLAB
% structure
%
% Usage:
%   >>  et = parsesmi(inputFile, outputMatFile, triggerKeyword);
%
% Inputs:
%   inputFile      - path to text file, converted from raw ET recording
%   outputMatFile  - MATLAB (.mat) file to save "et" structure in
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
%   >> ET = parsesmi('/myrawdata/SAMPLES.txt','/myfiles/SAMPLES.mat')
%   Optional with keyword:
%   >> ET = parsesmi('/myrawdata/SAMPLES.txt','/myfiles/SAMPLES.mat','MYKEYWORD')
%
% Outputs:
%   et     - eyetracker structure
%
% See also:
%   pop_parsesmi,  pop_importeyetracker
%
% Author: ur
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

function et = parsesmi(inputFile,outputMatFile,triggerKeyword)

if nargin < 2
    help(mfilename);
    return;
end;

%% feedback
fprintf('\n\n%s(): Parsing Sensomotoric Instruments (IView X/RED) raw data. This may take a moment...',mfilename)
fprintf('\n-- reading txt...')

fid = fopen(inputFile,'r');
if fid == -1
   error('\nThe source file ''%s'' was not found.\n',inputFile)
end
B = fread(fid,'*char')';
fclose(fid);
fprintf('\n\tDone.')

%% file specifications and comments
fprintf('\n-- getting comment lines and column headers...')
et.comments = regexp(B,'(##[^\r\n]*\r\n)','match');

%% build column header
et.colheader = regexp(B,'Time[^\r\n]*\r\n','match');

% clear 'Type' from header, as this information is not kept
et.colheader = strrep(et.colheader,sprintf('Type\t'),'');
et.colheader = regexp(et.colheader{:},'[^\t\r\n]*','match');
triggerColIndex = find(strcmp('Trigger',et.colheader));
fprintf('\n\tDone.')

%% get data
fprintf('\n-- getting data...')
dataLines = regexp(B,'(^\d*\tSMP[^\r\n]*)','match','lineanchors')';
% remove 'Type' information => keep only numerical data
dataLines = strrep(dataLines,sprintf('SMP\t'),'');

% fast conversion from cell of characters to array with blank-pads
dataChar = char(dataLines);

% sscanf treats numbers that follow directly over linebreaks without
% whitespace-characters as one. Therefore you occasionally need an
% additional space at the end of the array.
dataChar(:,end+1) = char(32);

% test once for number of columns
nColumns = length(sscanf(dataChar(1,:),'%f'));

%sscans string array by !column! for float-numbers, arranges in rows of
%nColumns, till end of string (inf). transpose back to column order.
et.data = sscanf(dataChar','%f',[nColumns,inf])';

clear dataLines dataChar
fprintf('\n\tDone.')


%% check for empty columns in data array
nColh = length(et.colheader);
nDat  = size(et.data,2);
if nColh ~= nDat
    column_mismatch = true;
else
    column_mismatch = false;
end

%% get messages
fprintf('\n-- getting messages...')
et.messages = regexp(B,'\d*\tMSG[^\r\n]*\r\n','match');
fprintf('\n\tDone.')

%% build table with synchronization events for eye tracker

% 1: check for messages containing special KEYWORD
if exist('triggerKeyword','var')
    
    if isempty(triggerKeyword)
        warning('\n\n%s(): You are using an empty string as the keyword (before the integer event number).',mfilename)
    end
    
    % take care of leading whitespaces as well as missing trailing whitespaces
    % handle "MYKEYWORD 123" as well as "MYKEYWORD123"
    
    test1 = regexp(et.messages,['(\d+)\s+MSG\s+\d+\s+#\sMessage:\s',triggerKeyword,'\s?(\d+)'],'tokens')';  
    test2 = [test1{:}]';    
    test3 = cellfun(@str2double,test2,'UniformOutput',false);
    et.event = cell2mat(test3);
    if isempty(et.event)
        warning('\n\n%s(): The keyword string you have provided (''%s'') was not found in any messages of the ET file (%s). Please check your messages and make sure they contain the specified keyword.',...
            mfilename, triggerKeyword, inputFile)
    end
    
% 2: check for TTL pulses (in separate data column)
elseif ~isempty(triggerColIndex)
    
    % parallel port pulses often last several samples, take rising edge of
    % pulse as trigger latency. Discard samples where zero was send to port
    et.event = cleantriggers( [et.data(:,1),et.data(:,triggerColIndex)]);
    
% no events found: user feedback
else
    warning('\n%s(): Did not find trigger pulses or special messages that could be used as events for synchronization.',mfilename)
    if ~exist('triggerKeyword','var')
        fprintf('\nIf you have sent messages with a special keyword, please provide the keyword when calling this function. See help.')
    end
end

clear test* triggerColIndex


%% Extract timestamps of *other* messages (not eyeevent, not keyword-messages)
% New feature, Jan-2018, OD
% Identify messages that begin with 'MSG' (e.g. user-send messages) so that 
% their timestamps can be later converted into "eeg time" in function 
% synchronize(). Extracting this information should be helpful for user 
% who have coded their experimental conditions (etc.) in the form of 
% (SMI) messages
if isfield(et,'messages') & ~isempty(et.messages)
    
    % go tru messages in et.messages
    remove_this_msg = zeros(1,size(et.messages,2));
    contains_keyword = [];
    for m = 1:size(et.messages,2)
        % does message line start with 'MSG'? --> include
        contains_msg = ~isempty(strfind(et.messages{m},'MSG'));
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
    tmp(:,find(remove_this_msg)) = []; % remove msg marked above
    
    if ~isempty(tmp)   
        % now extract numeric timestamp from these remaining msg
        %timestamps = regexp(tmp,'^-?\D+(\d+)','tokens','once','lineanchors');
        timestamps = regexp(tmp,'^(\d+)','tokens','once','lineanchors');
        timestamps = [timestamps{:}]';
        timestamps = cellfun(@str2double,timestamps,'UniformOutput',false);
        timestamps = cell2mat(timestamps);
        
        % store these messages in structure 'et.othermessages'
        for j = 1:size(tmp,2)
            et.othermessages(j).timestamp = timestamps(j); % numeric timestamps
            et.othermessages(j).msg       = tmp{j}; % the msg string
        end
        clear tmp test* edfColHeader
    end
else
    fprintf('\n\n%s(): Did not find any other messages (MSG) lines in the data\n.',mfilename)
end

%% save
fprintf('\n-- saving structure as mat...')
save(outputMatFile,'-struct','et')
fprintf('\n\tDone.')

if column_mismatch
    % recent IDF converter versions add extra column names in header
    fprintf('\n\nNote: Number of data column names for numerical data (%i) differed from number of actual data columns (%i).',nColh,nDat);
    fprintf('\nThis is not necessarily a problem, since some IDF converter versions add extra column names (e.g., \"AUX1\").');
end

fprintf('\n\nFinished parsing eye track.\nSaved .mat file to %s.\n',outputMatFile)