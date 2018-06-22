%  eegplugin_eye_eeg() - Preprocess, synchronize, and analyze simultaneously
%    recorded eye tracking and EEG data
%
%  Description:
%   EYE-EEG is an extension for EEGLAB written to facilitate joint analyses
%   of eye tracking and electroencephalographic (EEG) data. Among other
%   things, it can be used for fixation control, objective eye artifact
%   rejection, pupillometry, saccade and fixation detection, control of
%   microsaccades, eye tracker-supported ICA component selection, basic
%   oculomotor research, or computing fixation-related potentials (FRPs).
%
%   EYE-EEG was developed by Olaf Dimigen and Ulrich Reinacher at 
%   Humboldt University at Berlin, Germany, 2009-2018
%
%   >> web: http://www2.hu-berlin.de/eyetracking-eeg
%
%   Copyright (C): Olaf Dimigen with Ulrich Reinacher, 2009-2018
%   User feedback welcome: email: olaf.dimigen@hu-berlin.de
%
%   Project made possible by a grant from Deutsche Forschungsgemeinschaft
%   to Research Group 868, http://mbd.uni-potsdam.de
%
%   If you use functions of this toolbox, please cite this reference:
%   ------------------------------------------------------------------
%   Dimigen, O., Sommer, W., Hohlfeld, A., Jacobs, A., & Kliegl, R. (2011).
%   Coregistration of eye movements and EEG in natural reading: Analyses
%   & Review, J Exp Psychol: General, Vol. 140 (4), 552-572.
%
%   If you use the saccade/fixation detection, please additionally cite:
%
%   Engbert, R., & Mergenthaler, K. (2006). Microsaccades are triggered
%   by low retinal image slip, PNAS, Vol. 103 (18), 7192-7197
%
%   If you select independent components based on the sac/fix variance
%   ratio criterion, make sure to also cite:
%
%   Plöchl, M., Ossandon, J.P., & König, P. (2012). Combining EEG and
%   eye tracking: identification, characterization, and correction of
%   eye movement artifacts in electroencephalographic data. Frontiers in
%   Human Neuroscience, doi: 10.3389/fnhum.2012.00278
%
%   Please clarify that you have used the implementation of these methods
%   in EYE-EEG, as they may differ from the original versions in some
%   regards
%
%   Installation:
%   Simply copy the folder called "eye_eegX.XX" containing the toolbox's m-files
%   as a subdirectory into the EEGLAB plugin directory and restart EEGLAB.
%
%   GUI functions
%   ------------------------------------------------------------------
%   <a href="matlab:helpwin pop_applytochannels">pop_applytochannels</a>    - apply EEGLAB functions to selected channels
%   <a href="matlab:helpwin pop_detecteyemovements">pop_detecteyemovements</a> - detect saccades and fixations
%   <a href="matlab:helpwin pop_ploteyemovements">pop_ploteyemovements</a>   - plot properties of saccades and fixations
%   <a href="matlab:helpwin pop_importeyetracker">pop_importeyetracker</a>   - import & synchronized eye tracking data
%   <a href="matlab:helpwin pop_parseeyelink">pop_parseeyelink</a>       - preprocess eye tracking data from SR Research
%   <a href="matlab:helpwin pop_parsesmi">pop_parsesmi</a>           - preprocess eye tracking data from Sensomot. Instr.
%   <a href="matlab:helpwin pop_parsetobii">pop_parsetobii</a>           - preprocess eye tracking data from Tobii
%   <a href="matlab:helpwin pop_rej_eyecontin">pop_rej_eyecontin</a>      - reject cont. data with out-of-range eye track
%   <a href="matlab:helpwin pop_rej_eyeepochs">pop_rej_eyeepochs</a>      - reject epochs with out-of-range eye track
%   <a href="matlab:helpwin pop_checksync">pop_checksync</a>      - assess quality of ET/EEG synchronization via cross-correlation
%   <a href="matlab:helpwin pop_ploteventrate">pop_ploteventrate</a>      - plot rate of saccades or fixations relative to the epoch time-locking event
%   <a href="matlab:helpwin pop_overweightevents">pop_overweightevents</a>   - create optimized dataset for training the ICA
%   <a href="matlab:helpwin pop_eyetrackerica">pop_eyetrackerica</a>      - identify ICs that covary with eye tracker

%   Non-GUI functions
%   ----------------------------------------------------
%   <a href="matlab:helpwin addevents">addevents</a>              - add many events/urevents to cont. or epoched data
%   <a href="matlab:helpwin applytochannels">applytochannels</a>        - see pop function
%   <a href="matlab:helpwin detecteyemovements">detecteyemovements</a>     - see pop function
%   <a href="matlab:helpwin eyetrackerica">eyetrackerica</a>          - see pop function
%   <a href="matlab:helpwin geticvariance">geticvariance</a>          - compute IC variance inside/outside saccades
%   <a href="matlab:helpwin parseeyelink">parseeyelink</a>           - see pop function
%   <a href="matlab:helpwin parsesmi">parsesmi</a>               - see pop function
%   <a href="matlab:helpwin parsetobii">parsetobii</a>               - see pop function
%   <a href="matlab:helpwin rej_eyecontin">rej_eyecontin</a>          - see pop function
%   <a href="matlab:helpwin hist2d">hist2d</a>                 - 2-dimensional histogram
%   <a href="matlab:helpwin cleantrigger">cleantrigger</a>           - helper function
%   <a href="matlab:helpwin findmatchingevents">findmatchingevents</a>     - helper function
%   <a href="matlab:helpwin findsequence2">findsequence2</a>          - helper function
%   <a href="matlab:helpwin ploteyemovements">ploteyemovements</a>       - plot properties of saccades & fixations
%   <a href="matlab:helpwin synchronize">synchronize</a>            - synchronize ET/EEG data
%   <a href="matlab:helpwin checksync">checksync</a>                - assess quality of ET/EEG synchronization via cross-correlation
%   <a href="matlab:helpwin overweightevents">overweightevents</a>   - create optimized dataset for training the ICA

%   Non-GUI functions for saccade & fixation detection
%   ----------------------------------------------------
%   <a href="matlab:helpwin vecvel">vecvel</a>                 - compute eye velocities (*)
%   <a href="matlab:helpwin velthresh">velthresh</a>              - compute relative velocity-thresholds (*)
%   <a href="matlab:helpwin microsacc_plugin">microsacc_plugin</a>       - detect (micro)saccades (*)
%   <a href="matlab:helpwin binsacc">binsacc</a>                - identify binocular saccades (*)
%   <a href="matlab:helpwin saccpar">saccpar</a>                - compute saccade parameters (*)
%   <a href="matlab:helpwin mergesacc">mergesacc</a>              - handle temporal saccade clusters
%
%   (*) Copyright: Ralf Engbert et al., University of Potsdam
%   http://mbd.psych.uni-potsdam.de/~ralfweb, included with permission
%   Notes:
%   vecvel() was modified to include different data smoothing options
%   velthresh() outsources code that was originally part of microsacc().
%   microsacc() was therefore renamed to microsacc_plugin().
%
%   Third-party helper functions
%   ----------------------------------------------------
%   <a href="matlab:helpwin parseq">parseq</a>                 - Copyright: R. Kliegl & J. Laubrock, University of Potsdam
%   <a href="matlab:helpwin searchclosest">searchclosest</a>          - Copyright: Dr. Murtaza Khan, drkhanmurtaza@gmail.com
%
%   Functions that control the GUI dialogues
%   ----------------------------------------------------
%   <a href="matlab:helpwin dlg_applytochannels">dlg_applytochannels</a>
%   <a href="matlab:helpwin dlg_calcvisangle">dlg_calcvisangle</a>
%   <a href="matlab:helpwin dlg_detecteyemovements">dlg_detecteyenovements</a>
%   <a href="matlab:helpwin dlg_eyetrackerica">dlg_eyetrackerica</a>
%   <a href="matlab:helpwin dlg_loadeyetracker">dlg_loadeyetracker</a>
%   <a href="matlab:helpwin dlg_parser">dlg_parser</a>
%   <a href="matlab:helpwin dlg_ploteyemovements">dlg_ploteyemovements</a>
%   <a href="matlab:helpwin dlg_rej_eyecontin">dlg_rej_eyecontin</a>
%   <a href="matlab:helpwin dlg_rej_eyeepochs">dlg_rej_eyeepochs</a>
%   <a href="matlab:helpwin dlg_syncsettings">dlg_syncsettings</a>
%   <a href="matlab:helpwin dlg_ploteventrate">dlg_ploteventrate</a>
%   <a href="matlab:helpwin dlg_checksync">dlg_checksync</a>
%   ----------------------------------------------------
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% SMI (TM)     is a trademark of SensoMotoric Instruments GmbH, Germany
% Eyelink (TM) is a trademark of SR Research Ltd., Mississauga, Canada


function vers = eegplugin_eye_eeg(fig,try_strings,catch_strings)

vers = 'eye_eeg_v0.85';

% add subfolder for dialogues and other helpers
addpath(fullfile(fileparts(which(mfilename)),'internal'));

%% callbacks for menu items

% PARSE EYE TRACKER RAW DATA
cb_parsesmi         = [try_strings.no_check, '[ET LASTCOM]  = pop_parsesmi; clear ET;', catch_strings.add_to_hist];
cb_parseeyelink     = [try_strings.no_check, '[ET LASTCOM]  = pop_parseeyelink; clear ET;', catch_strings.add_to_hist];
cb_parsetobii       = [try_strings.no_check, '[ET LASTCOM]  = pop_parsetobii; clear ET;', catch_strings.add_to_hist];

% IMPORT & SYNCHRONIZE EYE TRACK
cb_importeyetracker = [try_strings.check_cont, '[EEG LASTCOM] = pop_importeyetracker(EEG);', catch_strings.new_and_hist];
cb_checksync        = [try_strings.check_data, '[EEG LASTCOM] = pop_checksync(EEG);', catch_strings.add_to_hist];

% PREPROCESS EYE TRACK
cb_rej_eyeepochs    = [try_strings.check_epoch, '[EEG LASTCOM] = pop_rej_eyeepochs(EEG);', catch_strings.new_and_hist];
cb_rej_eyecontin    = [try_strings.check_cont,  '[EEG LASTCOM] = pop_rej_eyecontin(EEG);', catch_strings.new_and_hist];

% APPLY EEGLAB FUNCTION TO SELECTED CHANNELS ONLY
check4chanlocs      = '[EEG LASTCOM] = eeg_checkset(EEG,''chanloc''); if ~isempty(LASTCOM), if LASTCOM(1) == -1, LASTCOM = ''''; return; end; end;';

cb_envelope_f1      = [try_strings.check_data, '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_eegfilt(EEG);'');', catch_strings.new_and_hist];
cb_envelope_f2      = [try_strings.check_data, '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_ERPLAB_butter1(EEG);'');', catch_strings.new_and_hist];
cb_envelope_f3      = [try_strings.check_data, '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_ERPLAB_polydetrend(EEG);'');', catch_strings.new_and_hist];
cb_envelope_f4      = [try_strings.check_data, '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_iirfilt(EEG);'');', catch_strings.new_and_hist];
cb_envelope_f5      = [try_strings.check_data, '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_eegfiltnew(EEG);'');', catch_strings.new_and_hist];
cb_envelope_f6      = [try_strings.check_data, '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_firma(EEG);'');', catch_strings.new_and_hist];
cb_envelope_f7      = [try_strings.check_data, '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_firpm(EEG);'');', catch_strings.new_and_hist];
cb_envelope_f8      = [try_strings.check_data, '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_firws(EEG);'');', catch_strings.new_and_hist];
cb_envelope_2       = [try_strings.check_data, '[EEG tmp LASTCOM] = pop_epoch(EEG); clear tmp;' catch_strings.new_and_hist];
cb_envelope_3       = [try_strings.check_data, '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_rmbase(EEG);'');' catch_strings.new_and_hist];
cb_envelope_4       = [try_strings.check_data, check4chanlocs '[EEG LASTCOM] = pop_applytochannels(EEG,[],''pop_timtopo(EEG);'');' catch_strings.add_to_hist];

% DETECT SACCADES & FIXATIONS
cb_detectem_1       = [try_strings.no_check,  '[EEG LASTCOM] = pop_detecteyemovements(EEG);', catch_strings.new_and_hist];

% PLOT SACCADES & FIXATIONS
cb_plotem_1         = [try_strings.no_check,  '[EEG LASTCOM] = pop_ploteyemovements(EEG);', catch_strings.add_to_hist];

% PLOT SACCADES & FIXATIONS
cb_plotrate         = [try_strings.no_check,  '[LASTCOM]     = pop_ploteventrate(EEG);', catch_strings.add_to_hist];

% CREATE OPTIMIZED ICA TRAINING DATA (OPTICAT)
% cb_overweight       = [try_strings.no_check,  '[EEG_overweighted LASTCOM] = pop_overweightevents(EEG);', catch_strings.add_to_hist];
cb_overweight       = [try_strings.no_check,  'EEG_overweighted = pop_eegfiltnew(EEG); [EEG_overweighted LASTCOM] = pop_overweightevents(EEG); [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_overweighted, size(ALLEEG,2),''setname'',''Overweighted training data'',''gui'',''off''); eeglab redraw;', catch_strings.add_to_hist];


% EYETRACKER-SUPPORTED ICA
cb_etICA            = [try_strings.check_ica, '[EEG vartable LASTCOM] = pop_eyetrackerica(EEG);', catch_strings.add_to_hist];


%% create GUI menus
W_MAIN = findobj('tag','EEGLAB');

% "EYETRACKER"
% attempt to make "Eyetracker" the last menu item before "Help"
% (note: menu position may change if more plugins are added, e.g. "ERPLAB")
eyetracker_m = uimenu( W_MAIN, 'Label', 'Eyetracker', 'tag','eyetracker');

try
    nitems = get(eyetracker_m,'position');
    set(eyetracker_m,'position',nitems-1);
catch
    set(eyetracker_m,'position',nitems);
end

% set "Eyetracker" menu to default color for plugin menus
icadefs; % get PLUGINMENUCOLOR
set(eyetracker_m, 'foregroundcolor', PLUGINMENUCOLOR);


% "PARSE EYETRACKER RAW DATA"
parseEyetracker_m = uimenu( eyetracker_m, 'Label', 'Parse eyetracker raw data', 'Separator', 'on');
uimenu( parseEyetracker_m, 'Label', 'text file from Eyelink', 'CallBack', cb_parseeyelink);
uimenu( parseEyetracker_m, 'Label', 'text file from SMI', 'CallBack', cb_parsesmi);
uimenu( parseEyetracker_m, 'Label', 'text file from Tobii', 'CallBack', cb_parsetobii);

% "IMPORT & SYNCHRONIZE ET"
uimenu( eyetracker_m, 'Label', 'Import & synchronize ET', 'CallBack', cb_importeyetracker, 'Separator', 'off', 'userdata', 'epoch:off; ');
uimenu( eyetracker_m, 'Label', 'Evaluate synchronization (cross-corr.)', 'CallBack', cb_checksync, 'Separator', 'off', 'userdata', 'startup:off;'); % epoch:off;'  

% APPLY FUNCTION TO SELECTED CHANNELS
helpers_m = uimenu( eyetracker_m, 'Label', 'Apply function to selected channels', 'Separator', 'on',  'userdata', 'startup:off;');
filters_m = uimenu( helpers_m, 'Label', 'Filter the data (selected channels)', 'Separator', 'on',  'userdata', 'startup:off;');
uimenu( filters_m, 'Label', 'Basic FIR filter (legacy)', 'CallBack', cb_envelope_f1, 'Separator', 'off', 'userdata', 'startup:off;');

% test whether filter functions from bundled plugins are available
% ...from ERPLAB
if exist('pop_ERPLAB_butter1','file')
    uimenu( filters_m, 'Label', 'ERPLAB Butterworth filter'   , 'CallBack', cb_envelope_f2, 'Separator', 'on',  'userdata', 'startup:off;  ');
end
if exist('pop_ERPLAB_polydetrend','file')
    uimenu( filters_m, 'Label', 'ERPLAB Polynomial Detrending', 'CallBack', cb_envelope_f3, 'Separator', 'off', 'userdata', 'startup:off; ');
end

% ...from other plugins (not yet added to EEGLAB path)
eegplugpath = fullfile(fileparts(which('eeglab')),'/plugins');
try
    placeHolder = '';
    % placeHolder = ' '; % use this if error occurs with exotic MacOS/Unix versions
    
    % folders for IIR and FIR filter plugins
    dir_iir = dir([eegplugpath placeHolder '/iirfilt*']);
    dir_fir = dir([eegplugpath placeHolder '/firfilt*']);
    
    % iirfilt plugin installed? if so, add menu entry
    if ~isempty(dir_iir)
        filtName = 'iirfilt.m';
        if exist(fullfile(eegplugpath,dir_iir(1).name,filtName),'file')
            uimenu( filters_m, 'Label', 'Short non-linear IIR filter', 'CallBack', cb_envelope_f4, 'Separator', 'on', 'userdata', 'startup:off; ');
        end
    else
        %warning('iirfilt plugin not found in plugin folder)')
        % #bugfix OD: removed this warning in v0.44, since iirfilt.m is not "shipped" with EEGLAB per default anymore
    end
    
    % firfilt plugin installed? if so, add menu entries
    if ~isempty(dir_fir)
        filtName      = {  'pop_eegfiltnew.m' 'pop_firws.m' 'pop_firpm.m' 'pop_firma.m'};
        filtMenuLabel = { 'Basic FIR filter (new)', 'Windowed sinc FIR filter' 'Parks-McClellan (equiripple) FIR filter' 'Moving average FIR filter'};
        filtMenuCb    = { cb_envelope_f5 cb_envelope_f8 cb_envelope_f7 cb_envelope_f6 };
        for loop = 1:length(filtName)
            if exist(fullfile(eegplugpath,dir_fir(1).name,filtName{loop}),'file')
                uimenu( filters_m, 'Label', filtMenuLabel{loop}, 'CallBack', filtMenuCb{loop}, 'Separator', 'off', 'Position' , loop);
            end
        end
    else
        warning('firfilt plugin not found in plugin folder')
    end
    
catch dir_error
    fprintf('\n%s(): Could not add filter functions to plugin menu. Skipping this step.\n',mfilename)
end

uimenu( helpers_m, 'Label', 'Extract epochs (without baseline removal)', 'CallBack', cb_envelope_2, 'Separator', 'off', 'userdata', 'startup:off;');
uimenu( helpers_m, 'Label', 'Remove baseline (selected channels)'      , 'CallBack', cb_envelope_3, 'Separator', 'off', 'userdata', 'startup:off;');
uimenu( helpers_m, 'Label', 'Channel ERPs - with scalp maps'           , 'CallBack', cb_envelope_4, 'Separator', 'off', 'userdata', 'startup:off;');

% "REJECT DATA BASED ON ET"
rejectEyetrack_m = uimenu( eyetracker_m, 'Label', 'Reject data based on eye track'      , 'Separator', 'on',  'userdata', 'startup:off;  ');
uimenu( rejectEyetrack_m, 'Label', 'Reject bad cont. data', 'CallBack', cb_rej_eyecontin, 'Separator', 'off', 'userdata', 'startup:off; epoch:off;');
uimenu( rejectEyetrack_m, 'Label', 'Reject bad epochs',     'CallBack', cb_rej_eyeepochs, 'Separator', 'off', 'userdata', 'startup:off; continuous:off;');

% "DETECT SACCADES & FIXATIONS"
uimenu( eyetracker_m, 'Label', 'Detect saccades & fixations', 'CallBack', cb_detectem_1, 'Separator', 'off', 'userdata', 'startup:off;  ');

% "PLOT SACCADES & FIXATIONS"
uimenu( eyetracker_m, 'Label', 'Plot eye movement properties', 'CallBack', cb_plotem_1, 'Separator', 'off', 'userdata', 'startup:off;  ');

% "PLOT SACCADES & FIXATIONS"
uimenu( eyetracker_m, 'Label', 'Plot eye movement rate', 'CallBack', cb_plotrate, 'Separator', 'off', 'userdata', 'startup:off;  ');

% "REJECT COMPONENTS WITH ET"
uimenu( eyetracker_m, 'Label', 'Create overweighted ICA training data (OPTICAT)', 'CallBack', cb_overweight, 'Separator', 'off', 'userdata', 'startup:off;  ');

% "REJECT COMPONENTS WITH ET"
uimenu( eyetracker_m, 'Label', 'Reject components with eyetracker', 'CallBack', cb_etICA, 'Separator', 'off', 'userdata', 'startup:off;  ');

% "ABOUT THE TOOLBOX"
uimenu( eyetracker_m, 'Label', 'About this toolbox', 'CallBack', 'pophelp(''eegplugin_eye_eeg.m'');', 'Separator', 'on', 'userdata', 'startup:on;');