% ps_08_shuffle_pauses_calculate_listener_blink_pause_proportions( ...
%     CONDITION_TYPE, STIM_KEY_TBL, N_ITERATIONS, ...
%     PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_PAUSE_EPOCHS, ...
%     PATHOUT_SHUFFLED_PROPORTIONS)
%
% This function shuffles the pauses within each epoch and calculates the
% propotion of shuffled pauses of one type of condition that contain
% listener blinks. The function saves a vector of proportions for every
% participant.
%
% Inputs:
% CONDITION_TYPE                - A character vector specifying the type of
%                                 condition to calculate the proportions
%                                 for. This must match an abbreviation in
%                                 the stimulus key.
% STIM_KEY_TBL                  - The stimulus key table to help extract
%                                 proportions of the selected type of
%                                 condition.
% N_ITERATIONS                  - The number of times the shuffling should
%                                 be repeated.
% PATHIN_LISTENER_BLINK_EPOCHS  - The path to where the listener blink
%                                 epochs are stored. 
% PATHIN_PAUSE_EPOCHS           - The path to where the speech-pause epochs
%                                 are stored.
% PATHOUT_SHUFFLED_PROPORTIONS  - The path to where vectors containing the
%                                 listener blink-pause proportions after
%                                 the shuffling should be stored.



% -------------------------------------------------------------------------
%% ** ps_08_shuffle_pauses_calculate_listener_blink_pause_proportions *****
% -------------------------------------------------------------------------


function ps_08_shuffle_pauses_calculate_listener_blink_pause_proportions( ...
    CONDITION_TYPE, STIM_KEY_TBL, N_ITERATIONS, ...
    PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_PAUSE_EPOCHS, ...
    PATHOUT_SHUFFLED_PROPORTIONS)

%% Prepare to calculate proportions: --------------------------------------

% Get the names of all the listener blink epoch files:
BlinkEpochFiles = dir( ...
                  fullfile(PATHIN_LISTENER_BLINK_EPOCHS, '*blink_epochs.mat'));

% Load the pause epochs:
load(fullfile(PATHIN_PAUSE_EPOCHS,'pause_epochs.mat'), 'AllStoriesPauses');

% Convert to facilitate the extraction of indices:
stimKeyCells = table2cell(STIM_KEY_TBL);

% Find the row and column indices to the conditions of interest:
relevantCondBool = find(strcmp(CONDITION_TYPE, stimKeyCells));
[conditionsIdx, storiesIdx] = ind2sub(size(stimKeyCells), relevantCondBool);

%% Loop over subjects and calculate proportions after shuffling: ----------

% Start stopwatch and open a progressbar:
timer = tic;
wb    = waitbar(0,'Shuffling pauses and computing listener blink-pause proportions...');

for f = 1:length(BlinkEpochFiles)

    %% Load the blinks and preallocate:

    % Get the identifier of one subject and load their blink epochs:
    subId = extractBefore(BlinkEpochFiles(f).name, '_blink_epochs');
    load(fullfile(PATHIN_LISTENER_BLINK_EPOCHS, BlinkEpochFiles(f).name), ...
          'AllStoriesListenerBlinks');

    % Preallocate:
    proportionsAfterShuffling = nan(1, N_ITERATIONS);

    % Notify:
    disp(['Beginning to shuffle the pauses for ', subId, ' ...']);
    disp(['Elapsed time: ', num2str(toc(timer))]);


    %% Loop to repeatedly shuffle the pauses:

    parfor i = 1:N_ITERATIONS

        % Shuffle the pauses, then count the number of blink-pauses and the
        % total number pauses for each condition:
        [nBlinkPausesVec, nPausesTotalVec] = ...
            local_shuffle_pauses_count_blink_pauses( ...
            AllStoriesListenerBlinks, AllStoriesPauses, ...
            storiesIdx, conditionsIdx);

        % Calculate the proportion of blink-pauses for this type of
        % condition by lumping data from all conditions of this type
        % together:
        proportion = sum(nBlinkPausesVec) / sum(nPausesTotalVec);

        % Store:
        proportionsAfterShuffling(i) = proportion;

    end % End of the shuffling loop


    %% Save the proportions after the shuffling:

    save( fullfile(PATHOUT_SHUFFLED_PROPORTIONS, ...
          [subId,'_shuffled_listener_blink_pause_proportions_', ...
                                         CONDITION_TYPE,'.mat']), ...
          "proportionsAfterShuffling");

    % Notify:
    disp(['Completed shuffling and calculating proportions for ', subId, ...
          ' for the ', CONDITION_TYPE, ' conditions.']);
    disp(['Elapsed time: ', num2str(toc(timer))]);
    waitbar(f/length(BlinkEpochFiles), wb);

end % End of the loop over subjects

%% Celebrate:

delete(wb);
disp('All of the shuffling and saving of listener blink-pauses has finished');
load('gong.mat'); 
sound(y, Fs);

end

% -------------------------------------------------------------------------
%% ************************ Local functions *******************************
% -------------------------------------------------------------------------

% Order of local functions:
% 1) local_shuffle_pauses_count_blink_pauses()
% 2) local_shuffle_pauses()

%% local_shuffle_pauses_count_blink_pauses() ---------------------------------------------

% This function shuffles the pauses. It then counts the number pauses that
% contain at least one listener blink for all the conditions of one type
% (e.g. audio-visual) by using the story and condition indices. The
% function outputs two vectors: one contains the number of blink-pauses in
% each of the conditions, the other contains the total number of pauses in
% each of the conditions.

function [nBlinkPausesVec, nPausesTotalVec] = ...
    local_shuffle_pauses_count_blink_pauses( ...
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

    % Shuffle the pause intervals and durations:
    shuffledPauseCells = local_shuffle_pauses(pausesLatsCells);

    % Check whether each pause contains at least one blink after the
    % shuffling:
    blinkPausesBool = cellfun(@(onePause) ...
        0 < length(intersect(onePause, listenerBlinkLats)), ...
        shuffledPauseCells );

    % Count the number of pauses with blinks:
    nBlinkPausesVec(c) = sum(blinkPausesBool);

    % Count the total number of pauses:
    nPausesTotalVec(c) = length(blinkPausesBool);

end
end


%% local_shuffle_pauses() -------------------------------------------------

% This function shuffles pauses around by shuffling the order of the pauses
% as well as the order of the inter-pause intervals. The interval between
% the start of the condition and the first pause, and the interval between
% the end of the condition and the last pause aren't changed.

function shuffledPauseCells = local_shuffle_pauses(pausesLatsCells)

% Determine the gap in samples between the onset of each pause and the
% previous pause or the start of the condition:
allGaps       = [];
lastPauseEnd  = 0;
for p = 1:length(pausesLatsCells)
    allGaps(p)    = pausesLatsCells{p}(1) - lastPauseEnd;
    lastPauseEnd  = pausesLatsCells{p}(end);
end

% Reset the sample latencies of each of the pauses to help move them
% around:
resetPauseCells = cellfun(@(onePauseCell) ...
                  onePauseCell(1,:) - onePauseCell(1), ...
                  pausesLatsCells, 'UniformOutput', false);

% Shuffle the order of the pauses so that the order of their durations is
% randomized:
shuffledPauseCells = resetPauseCells( randperm(length(resetPauseCells)) );

% Shuffle the order of gaps, but leave the first and last gaps in place as
% there may not be gaps at the edges of the condition:
randIdx       = randperm( length(allGaps) );
randIdx       = [1, randIdx(randIdx ~= 1)]; % So that the first gap is not moved
shuffledGaps  = allGaps(randIdx);

% Adjust the latencies of each of the pauses:
lastPauseEnd = 0;
for p = 1:length(shuffledPauseCells)
    shuffledPauseCells{p} = shuffledPauseCells{p} + shuffledGaps(p) + lastPauseEnd;
    lastPauseEnd          = shuffledPauseCells{p}(end);
end

end