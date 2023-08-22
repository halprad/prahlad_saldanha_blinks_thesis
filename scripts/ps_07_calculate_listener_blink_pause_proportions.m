% ps_07_calculate_listener_blink_pause_proportions(...
%     CONDITION_TYPES, STIM_KEY_TBL, ...
%     PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_PAUSE_EPOCHS, PATHOUT_PROPORTIONS)
%
% This function calculates the proportions of speech-pauses that contain at
% least one listener blink for all the types of conditions. It saves a
% table with the proportions as a .csv file which can be imported into R,
% and as a .mat file which can be imported into MATLAB.
%
% Inputs:
% CONDITION_TYPES               - A cell array of character vectors where
%                                 each vector is the abbreviation of a type
%                                 of condition presented during the
%                                 experiment. (These must match with those
%                                 in STIM_KEY_TBL.)
% STIM_KEY_TBL                  - The stimulus key table to help extract
%                                 proportions of the select types of
%                                 conditions.
% PATHIN_LISTENER_BLINK_EPOCHS  - The path to where the listener blink
%                                 epochs are stored.
% PATHIN_PAUSE_EPOCHS           - The path to where the speech-pause epochs
%                                 are stored.
% PATHOUT_PROPORTIONS           - The path to where the tables of the
%                                 listener blink-pause proportions should
%                                 be stored.



% -------------------------------------------------------------------------
%% ********** ps_07_calculate_listener_blink_pause_proportions ************
% -------------------------------------------------------------------------

function ps_07_calculate_listener_blink_pause_proportions(...
    CONDITION_TYPES, STIM_KEY_TBL, ...
    PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_PAUSE_EPOCHS, PATHOUT_PROPORTIONS)

%% Prepare to calculate proportions: --------------------------------------

% Convert for convenience later:
stimKeyCells = table2cell(STIM_KEY_TBL);

% Get the names of all the listener blink epoch files:
BlinkEpochFiles = dir( ...
                  fullfile(PATHIN_LISTENER_BLINK_EPOCHS, '*blink_epochs.mat'));

% Load the pause epochs:
load(fullfile(PATHIN_PAUSE_EPOCHS,'pause_epochs.mat'), 'AllStoriesPauses');


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

        % Count the number of blink-pauses and the total number of pauses
        % for each condition:
        [nBlinkPausesVec, nPausesTotalVec] = local_count_blink_pauses( ...
            AllStoriesListenerBlinks, AllStoriesPauses, ...
            storiesIdx, conditionsIdx);

        % Calculate the proportion of blink-pauses for this condition type
        % by lumping data from all conditions of this type together:
        proportion = sum(nBlinkPausesVec) / sum(nPausesTotalVec);

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

listenerBlinkPauseProportionsTbl = cell2table(proportionsCells, ...
                 "VariableNames", {'Subject', 'Condition', 'Proportion'});

save(fullfile(PATHOUT_PROPORTIONS,'listener_blink_pause_proportions.mat'), ...
     "listenerBlinkPauseProportionsTbl");

writetable(listenerBlinkPauseProportionsTbl, ...
    fullfile(PATHOUT_PROPORTIONS,'listener_blink_pause_proportions.csv'));

disp(['All proportions of the listener blinks in pauses were calculated ', ...
     'and saved successfully']);

end % End of main function



% -------------------------------------------------------------------------
%% *********************** Local function *********************************
% -------------------------------------------------------------------------

%% local_count_blink_pauses() ---------------------------------------------

% This function counts the number pauses that contain at least one listener
% blink for all the conditions of one type (e.g. audio-visual) by using the
% story and condition indices. The function outputs two vectors: one
% contains the number of blink-pauses in each of the conditions, the other
% contains the total number of pauses in each of the conditions.

function [nBlinkPausesVec, nPausesTotalVec] = local_count_blink_pauses( ...
    AllStoriesListenerBlinks, AllStoriesPauses, storiesIdx, conditionsIdx)

% Preallocate space for all pause counts of this type of condition:
nBlinkPausesVec = nan(1, length(conditionsIdx));
nPausesTotalVec = nan(1, length(conditionsIdx));

for c = 1:length(conditionsIdx) 

    % Determine the story and condition indices:
    stryId = storiesIdx(c);
    condId = conditionsIdx(c);

    % Extract the blink latencies of one condition:
    listenerBlinkLats = AllStoriesListenerBlinks(stryId).Conditions(condId).blinkLats;

    % Extract all the pauses of one condition:
    pausesLatsCells = AllStoriesPauses(stryId).Conditions(condId).pauseLatsCells;

    % Check whether each pause contains at least one blink:
    blinkPausesBool = cellfun(@(onePause) ...
        0 < length(intersect(onePause, listenerBlinkLats)), ...
        pausesLatsCells );

    % Count the number of pauses with blinks:
    nBlinkPausesVec(c) = sum(blinkPausesBool);

    % Count the total number of pauses:
    nPausesTotalVec(c) = length(blinkPausesBool);

end % End of loop
end % End of local fn