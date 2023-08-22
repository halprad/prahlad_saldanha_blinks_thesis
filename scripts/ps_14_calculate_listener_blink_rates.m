function ps_14_calculate_listener_blink_rates(...
    CONDITION_TYPES, STIM_KEY_TBL, ...
    PATHIN_LISTENER_BLINK_EPOCHS, PATHOUT_RATES)

%% Prepare to calculate blinks rates: -------------------------------------

% Convert for convenience later:
stimKeyCells = table2cell(STIM_KEY_TBL);

% Get the names of all the listener blink epoch files:
BlinkEpochFiles = dir( ...
                  fullfile(PATHIN_LISTENER_BLINK_EPOCHS, '*blink_epochs.mat'));

%% Loop over subjects and calculate blink rates: --------------------------

countsCells = {};

for f = 1:length(BlinkEpochFiles)

    % Get the identifier of one subject and load their blink epochs:
    subId = extractBefore(BlinkEpochFiles(f).name, '_blink_epochs');
    load(fullfile(PATHIN_LISTENER_BLINK_EPOCHS, BlinkEpochFiles(f).name), ...
        'AllStoriesListenerBlinks');

     %% Loop over condition types and calculate rates for each type:

    for t = 1:length(CONDITION_TYPES)

        % Find the row and column indices to the conditions of interest:
        relevantCondsBool = find(strcmp(CONDITION_TYPES{t}, stimKeyCells));
        [conditionsIdx, storiesIdx] = ...
            ind2sub(size(stimKeyCells), relevantCondsBool);

        % Count the number of blinks for this type of condition:
        nBlinksVec = local_count_blinks( ...
            AllStoriesListenerBlinks, storiesIdx, conditionsIdx);

        % Calculate the average blink rate in minutes for this type of
        % condition (the condition was 30 seconds long):
        rate = mean(nBlinksVec) * 2;


        % Store the information:
        countsCells = [countsCells; ...
            {subId, CONDITION_TYPES{t}, rate}];

    end % End of loop over conditions
end % End of loop over subjects 

%% Save and notify: -------------------------------------------------------

listenerBlinkRatesTbl = cell2table(countsCells, ...
                 "VariableNames", {'Subject', 'Condition', 'Proportion'});

save(fullfile(PATHOUT_RATES,'listener_blink_rates.mat'), ...
     "listenerBlinkRatesTbl");

writetable(listenerBlinkRatesTbl, ...
    fullfile(PATHOUT_RATES,'listener_blink_rates.csv'));

disp(['All proportions of the listener blink rates were calculated ', ...
      'and saved successfully']);

end % End of main function


%% ******************** Local function ************************************

function nBlinksVec = local_count_blinks( ...
    AllStoriesListenerBlinks, storiesIdx, conditionsIdx)

% Preallocate space for all blink counts of this type of condition:
nBlinksVec = nan(1, length(conditionsIdx));

for c = 1:length(conditionsIdx) 

    % Determine the story and condition indices:
    stryId = storiesIdx(c);
    condId = conditionsIdx(c);

    % Extract the blink latencies of one condition:
    listenerBlinkLats = AllStoriesListenerBlinks(stryId).Conditions(condId).blinkLats;


    % Count the number of blinks:
    nBlinksVec(c) = length(listenerBlinkLats);

end 
end