function ps_10_shuffle_post_blinks_calculate_synched_blink_proportions( ...
    FRAME_RATE, SRATE, CONDITION_TYPE, STIM_KEY_TBL, N_ITERATIONS,...
    PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_POST_BLINK_EPOCHS, ...
    PATHOUT_SHUFFLED_PROPORTIONS)


%% Prepare to calculate proportions: --------------------------------------

% Get the names of all the listener blink epoch files:
BlinkEpochFiles = dir( ...
                  fullfile(PATHIN_LISTENER_BLINK_EPOCHS, '*blink_epochs.mat'));

% Load the post-blink interval epochs:
load(fullfile(PATHIN_POST_BLINK_EPOCHS,'post_speaker_blink_epochs.mat'), ...
    'AllStoriesPostBlinks');

% Convert to facilitate the extraction of indices:
stimKeyCells = table2cell(STIM_KEY_TBL);

% Create a variable to help convert listener blink time-stamps to be
% comparable to the time-stamps of the post-speaker-blinks:
rateConverter = FRAME_RATE / SRATE; 

% Find the row and column indices to the conditions of interest:
relevantCondBool = find(strcmp(CONDITION_TYPE, stimKeyCells));
[conditionsIdx, storiesIdx] = ind2sub(size(stimKeyCells), relevantCondBool);

%% Loop over subjects and calculate proportions after shuffling: ----------

% Start stopwatch and open a progressbar:
timer = tic;
wb    = waitbar(0,'Shuffling post-blink intervals and computing synched blink proportions...');

for f = 1:length(BlinkEpochFiles)

    %% Load the blinks and preallocate:

    % Get the identifier of one subject and load their blink epochs:
    subId = extractBefore(BlinkEpochFiles(f).name, '_blink_epochs');
    load(fullfile(PATHIN_LISTENER_BLINK_EPOCHS, BlinkEpochFiles(f).name), ...
          'AllStoriesListenerBlinks');

    % Preallocate:
    proportionsAfterShuffling = nan(1, N_ITERATIONS);

    % Notify:
    disp(['Beginning to shuffle the post-blink intervals for ', subId, ' ...']);
    disp(['Elapsed time: ', num2str(toc(timer))]);


    %% Loop to repeatedly shuffle the post-blink intervals:

    parfor i = 1:N_ITERATIONS

        % Shuffle the post-blink intervals, then count the number of
        % synched-blinks and the total number post-blink intervals for each
        % condition:
        [nSynchedBlinksVec, nPostBlinksTotalVec] = ...
            local_shuffle_post_blinks_count_synched_blinks( ...
            AllStoriesListenerBlinks, AllStoriesPostBlinks, ...
            storiesIdx, conditionsIdx, rateConverter);

        % Calculate the proportion of synched-blinks for this type of
        % condition by lumping data from all conditions of this type
        % together:
        proportion = sum(nSynchedBlinksVec) / sum(nPostBlinksTotalVec);

        % Store:
        proportionsAfterShuffling(i) = proportion;

    end % End of the shuffling loop


    %% Save the proportions after the shuffling:

    save( fullfile(PATHOUT_SHUFFLED_PROPORTIONS, ...
          [subId,'_shuffled_synched_blink_proportions_', ...
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
disp('All of the shuffling and saving of synched blinks has finished');
load('chirp.mat'); 
sound(y, Fs);

end

% -------------------------------------------------------------------------
%% ************************ Local functions *******************************
% -------------------------------------------------------------------------

% Order of local functions:
% 1) local_shuffle_post_blinks_count_synched_blinks()
% 2) local_shuffle_post_blinks()

%% local_shuffle_post_blinks_count_synched_blinks() ---------------------------------------------

% This function shuffles the post-blink intervals. It then counts the
% number intervals that contain at least one listener blink for all the
% conditions of one type (e.g. audio-visual) by using the story and
% condition indices. The function outputs two vectors: one contains the
% number of synched blinks in each of the conditions, the other contains
% the total number of post-blink intervals in each of the conditions.

function [nSynchedBlinksVec, nPostBlinksTotalVec] = ...
    local_shuffle_post_blinks_count_synched_blinks( ...
    AllStoriesListenerBlinks, AllStoriesPostBlinks, ...
    storiesIdx, conditionsIdx, rateConverter)

% Preallocate space for all counts of this type of condition:
nSynchedBlinksVec   = nan(1, length(conditionsIdx));
nPostBlinksTotalVec = nan(1, length(conditionsIdx));

for c = 1:length(conditionsIdx) 

    % Determine the story and condition indices:
    stryId = storiesIdx(c);
    condId = conditionsIdx(c);

    % Extract the blink latencies of one condition:
    listenerBlinkLats = AllStoriesListenerBlinks(stryId).Conditions(condId).blinkLats;

    % Change the time-stamps to match the video frame-rate:
    toRoundUp = round(listenerBlinkLats * rateConverter) < listenerBlinkLats * rateConverter;
    listenerBlinkLats( toRoundUp) = ceil (listenerBlinkLats (toRoundUp) * rateConverter);
    listenerBlinkLats(~toRoundUp) = floor(listenerBlinkLats(~toRoundUp) * rateConverter);

    % Extract all the post-blinks intervals of one condition:
    postBlinkLatsCells = ...
        AllStoriesPostBlinks(stryId).Conditions(condId).postBlinkLatsCells;

    % Shuffle the post-blink intervals:
    shuffledPostBlinkCells = local_shuffle_post_blinks(postBlinkLatsCells);

    % Check whether each post-blink interval contains at least one blink
    % after the shuffling:
    synchedBlinkBool = cellfun(@(onePostBlink) ...
        0 < length(intersect(onePostBlink, listenerBlinkLats)), ...
        shuffledPostBlinkCells );

    % Count the number of post-blink intervals with blinks:
    nSynchedBlinksVec(c)   = sum(synchedBlinkBool);

    % Count the total number of post-blink intervals:
    nPostBlinksTotalVec(c) = length(synchedBlinkBool);

end
end


%% local_shuffle_post_blinks() -------------------------------------------------

% This function shuffles post-blink intervals around by shuffling the order
% of the inter-post-blink intervals. The interval between the start of the
% condition and the first pause, and the interval between the end of the
% condition and the last pause aren't changed.

function shuffledPostBlinkCells = local_shuffle_post_blinks(postBlinkLatsCells)

% Determine the gap in samples between the onset of each post-blink
% interval and the previous post-blink interval or the start of the
% condition:
allGaps       = [];
lastPostBlinkEnd  = 0;
for b = 1:length(postBlinkLatsCells)
    allGaps(b)        = postBlinkLatsCells{b}(1) - lastPostBlinkEnd;
    lastPostBlinkEnd  = postBlinkLatsCells{b}(end);
end

% Reset the sample latencies of each of the post-blink intervals to help
% move them around:
resetPostBlinkCells = cellfun(@(onePostBlinkCell) ...
                      onePostBlinkCell(1,:) - onePostBlinkCell(1), ...
                      postBlinkLatsCells, 'UniformOutput', false);

% Shuffle the order of gaps, but leave the first and last gaps in place as
% there may not be gaps at the edges of the condition:
randIdx       = randperm( length(allGaps) );
randIdx       = [1, randIdx(randIdx ~= 1)]; % So that the first gap is not moved
shuffledGaps  = allGaps(randIdx);

% Adjust the time-stamps of each of the post-blink intervals:
lastPostBlinkEnd = 0;
for b = 1:length(resetPostBlinkCells)
    shuffledPostBlinkCells{b} = resetPostBlinkCells{b} + shuffledGaps(b) + lastPostBlinkEnd;
    lastPostBlinkEnd          = shuffledPostBlinkCells{b}(end);
end

end