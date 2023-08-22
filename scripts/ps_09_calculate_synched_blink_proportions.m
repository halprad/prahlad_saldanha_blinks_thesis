function ps_09_calculate_synched_blink_proportions(...
    FRAME_RATE, SRATE, CONDITION_TYPES, STIM_KEY_TBL, ...
    PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_POST_BLINK_EPOCHS, ...
    PATHOUT_PROPORTIONS)
%% Prepare to calculate proportions: --------------------------------------

% Convert for convenience later:
stimKeyCells = table2cell(STIM_KEY_TBL);

% Get the names of all the listener blink epoch files:
BlinkEpochFiles = dir( ...
                  fullfile(PATHIN_LISTENER_BLINK_EPOCHS, '*blink_epochs.mat'));

% Load the post-blink epochs:
load(fullfile(PATHIN_POST_BLINK_EPOCHS,'post_speaker_blink_epochs.mat'), 'AllStoriesPostBlinks');

% Create a variable to help convert listener blink time-stamps to be
% comparable to the time-stamps of the post-speaker-blinks:
rateConverter = FRAME_RATE / SRATE; 


%% Loop over subjects and calculate proportions: --------------------------

proportionsCells = {};

for f = 1:length(BlinkEpochFiles)

    % Get the identifier of one subject and load their blink epochs:
    subId = extractBefore(BlinkEpochFiles(f).name, '_blink_epochs');
    load(fullfile(PATHIN_LISTENER_BLINK_EPOCHS, BlinkEpochFiles(f).name), ...
        'AllStoriesListenerBlinks');

    %% Loop over condition types and calculate proportions for each type:

    for t = 1:length(CONDITION_TYPES)

        % Find the row and column indices to the conditions of interest:
        relevantCondsBool = find(strcmp(CONDITION_TYPES{t}, stimKeyCells));
        [conditionsIdx, storiesIdx] = ...
            ind2sub(size(stimKeyCells), relevantCondsBool);

        % Count the number of synched-blinks and the total number of
        % post-blinks for each condition:
        [nSynchedBlinksVec, nPostBlinksTotalVec] = local_count_synched_blinks( ...
            AllStoriesListenerBlinks, AllStoriesPostBlinks, ...
            storiesIdx, conditionsIdx, rateConverter);

        % Calculate the proportion of blink-pauses for this condition type
        % by lumping data from all conditions of this type together:
        proportion = sum(nSynchedBlinksVec) / sum(nPostBlinksTotalVec);

        % If the proportion cannot be used in a beta regression model,
        % adjust it so that is can:
        switch proportion
            case 0
                proportion = 0.001;
                warning('Detected a proportion of 0. Adjusting it to 0.001');
            case 1
                proportion = 0.999;
                warning('Detected a proportion of 1. Adjusting it to 0.999');
        end

        % Store the information:
        proportionsCells = [proportionsCells; ...
                           {subId, CONDITION_TYPES{t}, proportion}];

    end % End of loop over condition types
end % End of loop over subjects

%% Save and notify: -------------------------------------------------------

synchedBlinkProportionsTbl = cell2table(proportionsCells, ...
                 "VariableNames", {'Subject', 'Condition', 'Proportion'});

save(fullfile(PATHOUT_PROPORTIONS,'synched_blink_proportions.mat'), ...
     "synchedBlinkProportionsTbl");

writetable(synchedBlinkProportionsTbl, ...
    fullfile(PATHOUT_PROPORTIONS,'synched_blink_proportions.csv'));

disp(['All proportions of the listener blinks in post-speaker-' ...
      'blink-intervals were calculated and saved successfully']);

end % End of main function



% -------------------------------------------------------------------------
%% *********************** Local function *********************************
% -------------------------------------------------------------------------


%% local_count_synched_blinks() ---------------------------------------------

% This function counts the number post-blinks that contain at least one
% listener blink for all the conditions of one type (e.g. audio-visual) by
% using the story and condition indices. The function outputs two vectors:
% one contains the number of synched-blinks in each of the conditions, the
% other contains the total number of pauses in each of the conditions.

function [nSynchedBlinksVec, nPostBlinksTotalVec] = ...
    local_count_synched_blinks( ...
        AllStoriesListenerBlinks, AllStoriesPostBlinks, ...
            storiesIdx, conditionsIdx, rateConverter)

% Preallocate space for all counts of this type of condition:
nSynchedBlinksVec   = nan(1, length(conditionsIdx));
nPostBlinksTotalVec = nan(1, length(conditionsIdx));

for c = 1:length(conditionsIdx) 

    % Determine the story and condition indices:
    stryId = storiesIdx(c);
    condId = conditionsIdx(c);

    % Extract the blink time-stamps of one condition:
    listenerBlinkLats = AllStoriesListenerBlinks(stryId).Conditions(condId).blinkLats;

    % Change the time-stamps to match the video frame-rate:
    toRoundUp = round(listenerBlinkLats * rateConverter) < listenerBlinkLats * rateConverter;
    listenerBlinkLats( toRoundUp) = ceil (listenerBlinkLats (toRoundUp) * rateConverter);
    listenerBlinkLats(~toRoundUp) = floor(listenerBlinkLats(~toRoundUp) * rateConverter);

    % Extract all the post-blinks of one condition:
    postBlinkLatsCells = ...
        AllStoriesPostBlinks(stryId).Conditions(condId).postBlinkLatsCells;

    % Check whether each post-blink contains at least one blink:
    synchedBlinkBool = cellfun(@(onePostBlink) ...
        0 < length(intersect(onePostBlink, listenerBlinkLats)), ...
        postBlinkLatsCells );

    % Count the number of post-blinks with listener blinks:
    nSynchedBlinksVec(c)   = sum(synchedBlinkBool);

    % Count the total number of post-blinks:
    nPostBlinksTotalVec(c) = length(synchedBlinkBool);

end
end