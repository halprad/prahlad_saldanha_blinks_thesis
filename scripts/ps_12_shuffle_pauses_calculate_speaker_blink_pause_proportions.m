function ps_12_shuffle_pauses_calculate_speaker_blink_pause_proportions( ...
    SPEAKERS, STIM_KEY_TBL, N_ITERATIONS, FRAME_RATE, SRATE, COND_DUR, ...
    MIN_GAP, MIN_PAUSE, MAX_PAUSE, ...
    PATHIN_SPEAKER_BLINKS, PATHIN_PAUSES, ...
    PATHOUT_SHUFFLED_PROPORTIONS)

%% Make preparations: -----------------------------------------------------

% Get the names of all the listener blink epoch files:
SpeakerBlinkFiles = dir(fullfile(PATHIN_SPEAKER_BLINKS, '*labels.mat'));

% Load the pauses:
PauseFiles = dir(fullfile(PATHIN_PAUSES,'*pauses.mat'));

% Determine how many conditions there were for each story:
numConditions    = arrayfun(@(storyNo) ...
                   local_check_num_conditions(STIM_KEY_TBL, storyNo), ...
                   1 : width(STIM_KEY_TBL));

% Determine which stories each speaker spoke for:

storyNames     = STIM_KEY_TBL.Properties.VariableNames;

speakerStoryNosCells = cellfun(@(oneSpeaker) ...
    local_match_speakers_2_story_nos(oneSpeaker, storyNames), ...
    SPEAKERS, 'UniformOutput', false);

% Determine to which story number each pause file (and consequently each
% blink file) corresponds:
pauseFileStoryNo = arrayfun(@(pauseFile) ...
                   local_match_pause_file_2_story_no(pauseFile, storyNames), ...
                   PauseFiles);

% Create a variable to help convert pause time-stamps to be comparable to
% the time-stamps of the speaker-blinks:
rateConverter = FRAME_RATE / SRATE;  


%% Loop over speakers, shuffle the pauses, and calculate proportions: -----

wb = waitbar(0, 'Shuffling pauses and calculating proportions of speaker-blink pauses.');
timer = tic;
for sp = 1:length(SPEAKERS)

    % Preallocate
    nBlinkPauses = nan(length(speakerStoryNosCells{sp}), N_ITERATIONS);
    nPausesTotal = nan(length(speakerStoryNosCells{sp}), N_ITERATIONS);

    for st = 1:length(speakerStoryNosCells{sp})
        
        % Get the index to the appropriate story and file:
        storyNo     = speakerStoryNosCells{sp}(st);
        currFileIdx = storyNo == pauseFileStoryNo;

        % Load the speaker blinks:
        load(fullfile(PATHIN_SPEAKER_BLINKS, SpeakerBlinkFiles(currFileIdx).name), ...
            'labeledFrames');

        % Load the pauses:
        load(fullfile(PATHIN_PAUSES, PauseFiles(currFileIdx).name), 'isPause');

        % Trim them so that they match the lengths of the stories shown
        % during the experiment:
        speakerBlinkLats = labeledFrames(labeledFrames <= ...
                           FRAME_RATE*COND_DUR*numConditions(storyNo));
        isPause          = isPause(1 : SRATE*COND_DUR*numConditions(storyNo));

        % Exclude short pauses sliced between conditions, and extract a
        % cell array of the pause time-stamps and the corrected pause
        % time-series:
        pauseLatsCells = local_extract_pause_timestamps( ...
                                                     isPause, SRATE, MIN_GAP, ...
                                                     MIN_PAUSE, MAX_PAUSE);

        
        parfor i = 1:N_ITERATIONS

            % Shuffle the pause intervals and durations:
            shuffledPauseCells = local_shuffle_pauses(pauseLatsCells);

            % Convert the pause time-stamps to match the video frame rate:
            shuffledPauseCellsFrate = cellfun(@(onePauseLats) ...
                local_convert_srate_2_frate(onePauseLats, rateConverter), ...
                shuffledPauseCells, 'UniformOutput', false);

            % Check whether each pause contains at least one blink after the
            % shuffling:
            blinkPausesBool = cellfun(@(onePause) ...
                0 < length(intersect(onePause, speakerBlinkLats)), ...
                shuffledPauseCellsFrate );

            % Count the number of pauses with blinks:
            nBlinkPauses(st,i) = sum(blinkPausesBool);

            % Count the total number of pauses:
            nPausesTotal(st,i) = length(blinkPausesBool);

        end % End of shuffling loop
    end % End of loop over stories for a given speaker

    %% Save the proportions after the shuffling:

    % Calculate the proportions:
    proportionsAfterShuffling = sum(nBlinkPauses, 1) ./ sum(nPausesTotal, 1);

    save( fullfile(PATHOUT_SHUFFLED_PROPORTIONS, ...
          [SPEAKERS{sp},'_shuffled_speaker_blink_pause_proportions.mat']), ...
          "proportionsAfterShuffling");

    % Notify:
    disp(['Completed shuffling and calculating proportions for ', SPEAKERS{sp}]);
    disp(['Elapsed time: ', num2str(toc(timer))]);
    waitbar(sp/length(SPEAKERS), wb);

end % End of the loop over speakers

%% Celebrate:

delete(wb);
disp('All of the shuffling and saving of speaker blink-pauses has finished');
load('train.mat'); 
sound(y, Fs);

end % End of the main function

% -------------------------------------------------------------------------
%% ************************ local functions *******************************
% -------------------------------------------------------------------------


%% local_check_num_conditions() -------------------------------------------

% This function checks how many conditions were presented for a given
% story.

function numConditions = local_check_num_conditions(STIM_KEY_TBL, storyNo)
isNonCondition    = strcmp(STIM_KEY_TBL{:,storyNo}, '.');
firstNonCondition = find(isNonCondition, 1, 'first');
numConditions     = firstNonCondition - 1;
if isempty(numConditions)
    numConditions = length(STIM_KEY_TBL{:, storyNo});
end
end

%% local_match_speakers_2_story_nos()

function speakerStoryNos = local_match_speakers_2_story_nos( ...
                                        oneSpeaker, storyNames)

isSpeakersStory = contains(storyNames, oneSpeaker);
speakerStoryNos = find(isSpeakersStory);

end


%% local_match_pause_file_2_story_no() ------------------------------------

% This function matches the index of a pause file to the story-video which
% has the same name. As a result, the story number (in the experiment) the
% pause file corresponds to is determined.

function pauseFileStoryNo = local_match_pause_file_2_story_no( ...
                                            pauseFile, storyNames)

[~, pauseFileName]    = fileparts(pauseFile.name);
pauseFileName         = extractBefore(pauseFileName,'_pauses');
pauseFileStoryNo      = find(strcmp(pauseFileName, storyNames));
end


%% local_extract_pause_timestamps() ---------------------------------------

% This function extracts a cell array containing the time-stamps of the
% pauses in samples. The time-stamps are relative to the start of each
% epoch. Each cell corresponds to one pause. The function makes sure that
% no two consecutive pauses are too close to each other; it eliminates
% close pauses. Also, overly long pauses and short pauses sliced between
% conditions get eliminated by the function. The function also outputs a
% time-series of the corrected pauses.

function pauseLatsCells = local_extract_pause_timestamps( ...
                                    isPause, SRATE, MIN_GAP, ...
                                            MIN_PAUSE, MAX_PAUSE)

pauseLats    = find(isPause);

%% Determine where the ends of the pauses are:

% A time-stamp more than a sample ahead of the previous one is part of a
% separate pause, hence:
isPauseEnd   = 1 < (pauseLats(2:end) - pauseLats(1:end-1));

% The last stamp will always be the end of a pause.
isPauseEnd   = [isPauseEnd, true]; 

pauseEndLats   = find(isPauseEnd);

%% Segregate the pause time-stamps into different cells:

pauseLatsCells = {};
epoTill        = 0;
for p = 1:length(pauseEndLats)
    epoFrom    = epoTill + 1;
    epoTill    = pauseEndLats(p);
    pauseLatsCells{p} = pauseLats(epoFrom:epoTill);
end

%% Get rid of close pauses:

currIdx = 2;
numRemaining = length(pauseLatsCells);
while currIdx <= numRemaining - 1

    % Determine if the next pause is too close and eliminate accordingly:
    currentEndSecs = pauseLatsCells{currIdx}    (end) / SRATE;
    nextStartSecs  = pauseLatsCells{currIdx + 1}  (1) / SRATE;
    isNext2Close   = MIN_GAP > nextStartSecs - currentEndSecs;

    if isNext2Close
        pauseLatsCells(currIdx + 1) = [];
    end

    % Determine if the prior pause is too close and eliminate accordingly:
    priorEndSecs     = pauseLatsCells{currIdx - 1}(end) / SRATE;
    currentStartSecs = pauseLatsCells{currIdx}      (1) / SRATE;
    isPrior2Close    = MIN_GAP > currentStartSecs - priorEndSecs;

    if isPrior2Close
        pauseLatsCells(currIdx - 1) = [];
    end

    % Prepare to check the next pause:
    currIdx =  currIdx + 1;
    numRemaining = length(pauseLatsCells);
end

%% Get rid of very long pauses and very short pauses:

pauseDurationsSecs   = cellfun(@length, pauseLatsCells) / SRATE;  
pauseLatsCells       = pauseLatsCells(pauseDurationsSecs <= MAX_PAUSE & ...
                                      pauseDurationsSecs >= MIN_PAUSE );

end % End of local function

%% local_convert_srate_2_frate() ------------------------------------------

function onePauseLatsFrate = local_convert_srate_2_frate( ...
                             onePauseLats, rateConverter)

onePauseLatsFrate = nan(size(onePauseLats));

% Change the time-stamps to match the video frame-rate:
toRoundUp = round(onePauseLats * rateConverter) < onePauseLats * rateConverter;
onePauseLatsFrate( toRoundUp) = ceil (onePauseLats (toRoundUp) * rateConverter);
onePauseLatsFrate(~toRoundUp) = floor(onePauseLats(~toRoundUp) * rateConverter);

% Get rid of duplicate time-stamps:
onePauseLatsFrate = unique(onePauseLatsFrate);

end


%% local_shuffle_pauses() -------------------------------------------------

% This function shuffles pauses around by shuffling the order of the pauses
% as well as the order of the inter-pause intervals. The interval between
% the start of the story and the first pause, and the interval between
% the end of the story and the last pause aren't changed.

function shuffledPauseCells = local_shuffle_pauses(pauseLatsCells)

% Determine the gap in samples between the onset of each pause and the
% previous pause or the start of the condition:
allGaps       = [];
lastPauseEnd  = 0;
for p = 1:length(pauseLatsCells)
    allGaps(p)    = pauseLatsCells{p}(1) - lastPauseEnd;
    lastPauseEnd  = pauseLatsCells{p}(end);
end

% Reset the sample latencies of each of the pauses to help move them
% around:
resetPauseCells = cellfun(@(onePauseCell) ...
                  onePauseCell(1,:) - onePauseCell(1), ...
                  pauseLatsCells, 'UniformOutput', false);

% Shuffle the order of the pauses so that the order of their durations is
% randomized:
shuffledPauseCells = resetPauseCells( randperm(length(resetPauseCells)) );

% Shuffle the order of gaps, but leave the first and last gaps in place as
% there may not be gaps at the edges of the story:
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