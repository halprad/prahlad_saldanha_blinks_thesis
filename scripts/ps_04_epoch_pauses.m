% ps_04_epoch_pauses( ...
%     SRATE, COND_DUR, MIN_GAP, MIN_PAUSE, MAX_PAUSE, STIM_KEY_TBL, ...
%     PATHIN_PAUSES, PATHOUT_PAUSE_EPOCHS)
%
% This function takes continus pause vectors for each story and epochs them
% into conditions.
%
% Pauses too close together get eliminated. Also, pause which are too long,
% and short pauses sliced between conditions get eliminated.
%
% The function saves a structure array with a structure for every story,
% each of which contains a structure array for every condition. In each of
% these sub-structures, a time-series of whether each sample was a pause is
% saved, along with a cell array of the time-stamps of different pauses
% segregated into different cells.
%
% Inputs:
% SRATE                 - The EEG sample rate in hertz.
% COND_DUR              - The duration of each condition in seconds.
% MIN_GAP               - The minimum acceptable gap between pauses in
%                         seconds.
% MIN_PAUSE             - The minimum acceptable duration of a pause in
%                         seconds.
% MAX_PAUSE             - The maximum acceptable duration of a pause in
%                         seconds.
% STIM_KEY_TBL          - The stimulus key table to help organize the data.
% PATHIN_PAUSES         - The path to where the pause files are stored.
% PATHOUT_PAUSE_EPOCHS  - The path to where the pause epochs should be
%                         stored.



% -------------------------------------------------------------------------
%% ********************** ps_04_epoch_pauses ******************************
% -------------------------------------------------------------------------

function ps_04_epoch_pauses( ...
    SRATE, COND_DUR, MIN_GAP, MIN_PAUSE, MAX_PAUSE, STIM_KEY_TBL, ...
    PATHIN_PAUSES, PATHOUT_PAUSE_EPOCHS)


% Determine the interval between samples:
sampInterval   = 1000/SRATE; % In ms


%% Get the names of the pause files and useful info about them: -----------

PauseFiles     = dir(fullfile(PATHIN_PAUSES,'*pauses.mat'));

storyNames     = STIM_KEY_TBL.Properties.VariableNames;

% Determine to which story number each pause file corresponds:
pauseFileStoryNo = arrayfun(@(pauseFile) ...
                   local_match_pause_file_2_story_no(pauseFile, storyNames), ...
                   PauseFiles);

% Determine how many conditions there were for each story:
numConditions    = arrayfun(@(storyNo) ...
                   local_check_num_conditions(STIM_KEY_TBL, storyNo), ...
                   1 : width(STIM_KEY_TBL));


%% Loop over all stories and epoch the pauses by condition: ---------------

AllStoriesPauses = struct();

for s = 1:length(storyNames)
    if numConditions(s) <= 0
        % There was one story which wasn't used and had 0 conditions, so:
        AllStoriesPauses(s).Conditions = [];
    else
        currFileIdx = pauseFileStoryNo == s;
        load(fullfile(PATHIN_PAUSES, PauseFiles(currFileIdx).name), 'isPause');
        epoTill = 0;
        for e = 1:numConditions(s)

            % Epoch one condition:
            epoFrom       = epoTill + 1;
            epoTill       = epoTill + COND_DUR*SRATE;
            isPauseEpoch  = isPause(epoFrom:epoTill);

            % Exclude short pauses sliced between conditions, and extract a
            % cell array of the pause time-stamps and the corrected pause
            % time-series:
            [pauseLatsCells, isPauseEpochCorr] = ...
                            local_extract_pause_timestamps( ...
                                  isPauseEpoch, SRATE, MIN_GAP, ...
                                        MIN_PAUSE, MAX_PAUSE);

            % Store the cell array of pause time-stamps:
            AllStoriesPauses(s).Conditions(e).pauseLatsCells  = ...
                pauseLatsCells;

            % Store the pause time-series:
            AllStoriesPauses(s).Conditions(e).isPause         = ...
                isPauseEpochCorr;
            
            % Store a time vector in ms:
            AllStoriesPauses(s).Conditions(e).times           = ...
                0 : sampInterval : COND_DUR*1000-sampInterval;

        end % End of the loop over conditions
    end % End of the conditional
end % End of the loop over stories

%% Save and notify: -------------------------------------------------------

save(fullfile(PATHOUT_PAUSE_EPOCHS,'pause_epochs.mat'), 'AllStoriesPauses');

disp('All of the pause epochs were extracted successfully');

end % End of the main function



% -------------------------------------------------------------------------
%% ********************** Local functions *********************************
% -------------------------------------------------------------------------

% Order of local functions:
% 1) local_match_pause_file_2_story_no()
% 2) local_check_num_conditions()
% 3) local_extract_pause_timestamps()


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


%% local_extract_pause_timestamps() ---------------------------------------

% This function extracts a cell array containing the time-stamps of the
% pauses in samples. The time-stamps are relative to the start of each
% epoch. Each cell corresponds to one pause. The function makes sure that
% no two consecutive pauses are too close to each other; it eliminates
% close pauses. Also, overly long pauses and short pauses sliced between
% conditions get eliminated by the function. The function also outputs a
% time-series of the corrected pauses.

function [pauseLatsCells, isPauseEpochCorr] = ...
                            local_extract_pause_timestamps( ...
                                    isPauseEpoch, SRATE, MIN_GAP, ...
                                            MIN_PAUSE, MAX_PAUSE)

pauseLats    = find(isPauseEpoch);

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

%% Get rid of long pauses and short pauses sliced between conditions:

pauseDurationsSecs   = cellfun(@length, pauseLatsCells) / SRATE;  
pauseLatsCells       = pauseLatsCells(pauseDurationsSecs <= MAX_PAUSE & ...
                                      pauseDurationsSecs >= MIN_PAUSE );

%% Create a time-series containing only the remaining pauses:
isPauseEpochCorr                      = false(size(isPauseEpoch));
isPauseEpochCorr([pauseLatsCells{:}]) = true;

end