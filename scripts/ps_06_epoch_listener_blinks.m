% ps_06_epoch_listener_blinks(...
%     COND_DUR, SRATE,...
%     PATHIN_LISTENER_BLINKS, PATHOUT_LISTENER_BLINK_EPOCHS)
% 
% This function epochs listener blinks.
%
% The function saves a structure array with a structure for every story,
% each of which contains a structure array for every condition. In each of
% these sub-structures are stored the time-stamps of the peaks of listener
% blinks, and the ICA signal that was used to find them.
%
% Inputs:
% COND_DUR                       - The duration of each condition in 
%                                  seconds.
% SRATE                          - The EEG sample rate in hertz.
% PATHIN_LISTENER_BLINKS         - The path to where the listener blink
%                                  files are stored.
% PATHOUT_LISTENER_BLINK_EPOCHS  - The path to where the listener blink
%                                  epochs should be stored.



% -------------------------------------------------------------------------
%% ***************** ps_06_epoch_listener_blinks **************************
% -------------------------------------------------------------------------


function ps_06_epoch_listener_blinks(...
    COND_DUR, SRATE,...
    PATHIN_LISTENER_BLINKS, PATHOUT_LISTENER_BLINK_EPOCHS)

%% Get the listener blink filenames (2 per subject): --------------

listenerBlinkFiles = dir(fullfile(PATHIN_LISTENER_BLINKS,'*blinks.mat'));


%% Loop over subjects and extract epochs (merge epochs across sessions): --

wb = waitbar(0, 'Extracting listener blink epochs...');

for f = 1:length(listenerBlinkFiles)

    %% Prepare the necessary variables to begin epoching:

    % Load the listener blink data for one session:
    load(fullfile(PATHIN_LISTENER_BLINKS, listenerBlinkFiles(f).name), ...
        'Blinks', 'BlinkFits', 'eegTimes', 'EegEvents');

    % Get time-stamps of the ends of conditions (these will be
    % half-a-sample after the end of every condition):
    EegEvents = EegEvents( strcmp({EegEvents.type},'boundary') );
    if isfield(EegEvents,'duration') % This deals with problematic datasets
        EegEvents = EegEvents( isnan([EegEvents.duration]) );
    end
    conditionEnds = [[EegEvents.latency], length(eegTimes)+0.5];

    % Get the time-stamps of the blinks and the signal used to find them:
    listenerBlinkLats = [BlinkFits.maxFrame];
    icaSignal         = Blinks.signalData.signal;


    %% Initialize a loop to keep extracting epochs until data runs out: 

    s = 1;          % s -> story count
    e = 0;          % e -> epoch (condition) count
    epoTill       = 0;        
    StoriesBlinks = struct();
    isDataLeft    = true;
    isNewStory    = false;

    while isDataLeft

        %% Set helper variables to the right indices:

        if isNewStory
            % Prepare to epoch right after the end of the previous story:
            epoTill = floor(conditionEnds(s));

            s = s + 1;
            e = 0;
        end

        e       = e + 1;
        epoFrom = epoTill + 1;
        epoTill = epoTill + COND_DUR*SRATE;


        %% Epoch the blink time-stamps, ICA signal, and time vector:

        % Adjust the epoched blink time-stamps so that they are relative to
        % the start of the epoch:
        blinkLatsEpoch = ...
            listenerBlinkLats( ...
                listenerBlinkLats >= epoFrom  &  listenerBlinkLats <= epoTill ) ...
                    - epoFrom + 1;

        % ICA signal:
        icaSignalEpoch = icaSignal(epoFrom:epoTill);

        % Adjust the epoched times to they start at 0:
        timeEpoch      = eegTimes(epoFrom:epoTill) - eegTimes(epoFrom);


        %% Store the epochs:

        StoriesBlinks(s).Conditions(e).blinkLats = blinkLatsEpoch;
        StoriesBlinks(s).Conditions(e).icaSignal = icaSignalEpoch;
        StoriesBlinks(s).Conditions(e).times     = timeEpoch;
        

        %% Decide how to proceed:

        % Check where the end of the next epoch will be and determine
        % whether there's data left or it's a new story:
        nextEpoTill = epoTill + COND_DUR*SRATE;
        isDataLeft  = length(eegTimes) >= nextEpoTill;
        isNewStory  = conditionEnds(s) < nextEpoTill;

    end % End of the epoching loop

    %% Organize the data:

    switch rem(f,2)
        case 1
            FirstSession = StoriesBlinks;
            continue;
        case 0
            SecondSession = StoriesBlinks;
    end

    % Merge the epochs from the two conditions together, including an empty
    % structure for the unused story:
    AllStoriesListenerBlinks = [FirstSession, ...
                                struct('Conditions',[]), ...
                                SecondSession];
    
    %% Save the listener blink epochs:

    subId = extract(listenerBlinkFiles(f).name, digitsPattern);
    save(fullfile(PATHOUT_LISTENER_BLINK_EPOCHS, ['Sub_', subId{:}, '_blink_epochs.mat']), ...
         'AllStoriesListenerBlinks');
    clear AllStoriesListenerBlinks;

    progress = f/length(listenerBlinkFiles);
    waitbar(progress, wb);

end % End of the loop iterating through blink files

delete(wb);
disp('All of the listener blink epochs have been extracted successfully.');

end % End of the function