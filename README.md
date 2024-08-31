# EYE-EEG Toolbox

Version 1.0 of the EYE-EEG extension for EEGLAB

## What is the EYE-EEG toolbox?

This extension allows you to do synchronize and process your combined eye-tracking and EEG.

The EYE-EEG toolbox is an extension for the open-source MATLAB toolbox EEGLAB developed to facilitate integrated analyses of electrophysiological and oculomotor data. The toolbox parses, imports, and synchronizes simultaneously recorded eye tracking data and adds it as extra channels to the EEG.

Saccades and fixations can be imported from the eye tracking raw data or detected with an adaptive velocity-based algorithm. Eye movements are then added as new time-locking events to EEGLAB's event structure, allowing easy saccade- and fixation-related EEG analysis (e.g., fixation-related potentials, FRPs). Alternatively, EEG data can be aligned to stimulus onsets and analyzed according to oculomotor behavior (e.g. pupil size, microsaccades) in a given trial. Saccade-related ICA components can be objectively identified based on their covariance with the electrically independent eye tracker.

EYE-EEG closely integrates into EEGLAB and adds a top-level menu called "Eyetracker" to EEGLAB. All functions can be accessed via this menu and are saved in EEGLAB's command history. Alternatively, functions can be called from the command line, providing advanced users with the option to use them in custom scripts. Using EEGLAB's export functions, integrated datasets may also be exported to other free toolboxes like Fieldtrip or Brainstorm.

## EYE-EEG Team
EYE-EEG was written by [Olaf Dimigen](https://olaf.dimigen.de) with major contributions by Ulrich Reinacher (until 2012). The toolbox has now been used in numerous [scientific publications](https://eyetracking-eeg.org/papersusing.html).

## Reference paper 
Dimigen, O., Sommer, W., Hohlfeld, A., Jacobs, A., & Kliegl, R. (2011). Coregistration of eye movements and EEG in natural reading: Analyses & Review. Journal of Experimenta Psychology: General, 140 (4), 552-572 

## More information, tutorials, and example data to play with:
For more information and a [tutorial](https://eyetracking-eeg.org/tutorial.html) on this software, visit [https://www.eyetracking-eeg.org](https://www.eyetracking-eeg.org).
