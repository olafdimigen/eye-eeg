function eventTable = cleantriggers(timeAndTrigger)
% small helper function called by parseeyelink() and parsesmi()
%
% Take only rising flank of trigger = first sample of trigger 
% pulses that last several samples
% Discard samples where parallel port input is zero
% timeAndTrigger = [timeStamps, triggerColumn]

% identify sequences of unchanging values with parseq()
triggerIndex = parseq(timeAndTrigger(:,2));

% keep only first timestamp/trigger value for each sequence
eventTable = [timeAndTrigger(triggerIndex(:,1),1), timeAndTrigger(triggerIndex(:,1),2)];

% clear 'zero' values
eventTable(eventTable(:,2)==0,:)=[];