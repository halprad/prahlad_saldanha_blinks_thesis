% ps_05_epoch_post_speaker_blinks( ...
%     FRAME_RATE, COND_DUR, MIN_GAP, POST_BLINK_SECS, ...
%     PATHIN_SPEAKER_BLINKS, PATHOUT_POST_BLINK_EPOCHS)
%
% This function uses the time-stamps of the speaker blinks to determine
% post-blink intervals. The function then epochs the post-blinks by the
% experimental conditions.
%
% Post-blinks sliced between conditions get eliminated. Additionally, for
% every pair of overlapping post-blinks, the latter one gets eliminated.
%
% The function saves a structure array with a structure for every story,
% each of which contains a structure array for every condition. In each of
% these sub-structures, a time-series of whether each sample was a
% post-blink is saved, along with a cell array of the time-stamps of
% different post-blinks segregated into different cells.
%
% Inputs:
% FRAME_RATE                         - The video frame rate in hertz.
% COND_DUR                           - The duration of each condition in 
%                                      seconds.
% MIN_GAP                            - The minimum acceptable gap between
%                                      post-blinks in seconds.
% POST_BLINK_SECS                    - The duration following a speaker
%                                      blink to consider as a post-blink in
%                                      seconds.
% STIM_KEY_TBL                       - The stimulus key table to help 
%                                      organize the data.
% PATHIN_SPEAKER_BLINKS              - The path to where the speaker blink
%                                      files are saved.
% PATHOUT_POST_BLINK_EPOCHS          - The path to where the post speaker
%                                      blink epochs should be saved.

% -------------------------------------------------------------------------
%% **************** ps_05_epoch_post_speaker_blinks ***********************
%--------------------------------------------------------------------------

function ps_05_epoch_post_speaker_blinks( ...
    FRAME_RATE, COND_DUR, MIN_GAP, POST_BLINK_SECS, STIM_KEY_TBL, ...
    PATHIN_SPEAKER_BLINKS, PATHOUT_POST_BLINK_EPOCHS)

% Determine the interval between frames:
frameInterval = 1000/FRAME_RATE; % In ms

%% Get the names of the speaker blink files and useful info about them: ---

SpeakerBlinkFiles = dir(fullfile(PATHIN_SPEAKER_BLINKS,'*labels.mat'));

storyNames        = STIM_KEY_TBL.Properties.VariableNames;

% Determine to which story number each speaker blink file corresponds:
blinkFileStoryNo        = arrayfun(@(speakerBlinkFile) ...
                          local_match_blink_file_2_story_no( ...
                                            speakerBlinkFile, storyNames), ...
                          SpeakerBlinkFiles);

% Determine how many conditions there were for each story:
numConditions           = arrayfun(@(storyNo) ...
                          local_check_num_conditions(STIM_KEY_TBL, storyNo), ...
                          1 : width(STIM_KEY_TBL));


%% Loop over all stories and epoch the post speaker blink intervals: ------

AllStoriesPostBlinks = struct();

for s = 1:length(storyNames)
    if numConditions(s) <= 0
        % There was one story which wasn't used and had 0 conditions, so:
        AllStoriesPostBlinks(s).Conditions = [];
    else
        currFileIdx = blinkFileStoryNo == s;
        load( fullfile(PATHIN_SPEAKER_BLINKS, SpeakerBlinkFiles(currFileIdx).name), ...
             'labeledFrames');

        % Preallocate: 
        isPostBlink    = false(1, FRAME_RATE*COND_DUR*numConditions(s));

        % Get the time-stamps in samples to the post-blink intervals:
        postBlinkCells = arrayfun(@(blinkLat) ...
                         blinkLat   : ...
                            blinkLat + round(FRAME_RATE*POST_BLINK_SECS) - 1, ...
                         labeledFrames, 'UniformOutput', false);        

        % Dump the time-stamps of all the intervals in one array:
        postBlinkLats  = [postBlinkCells{:}];

        % Keep only the blinks shown during the experiment:
        postBlinkLats(postBlinkLats > length(isPostBlink)) = []; 

        % Make this into a time-series:
        isPostBlink(postBlinkLats) = true;

        %% Epoch the time-series over conditions: -------------------------

        epoTill = 0;
        for e = 1:numConditions(s)
            epoFrom           = epoTill + 1;
            epoTill           = epoTill + COND_DUR*FRAME_RATE;
            isPostBlinkEpoch  = isPostBlink(epoFrom : epoTill);

            % Exclude post-blinks sliced between conditions and close
            % post-blinks. Then, extract a cell array of the post-blink
            % time-stamps and the corrected post-blink time-series:
            [postBlinkLatsCells, isPostBlinkEpochCorr] = ...
                    local_extract_post_blink_timestamps( ...
                          isPostBlinkEpoch, FRAME_RATE, MIN_GAP, ...
                                                POST_BLINK_SECS);

            % Store the cell array of post-blink time-stamps:
            AllStoriesPostBlinks(s).Conditions(e).postBlinkLatsCells  = ...
                postBlinkLatsCells;

            % Store the post-blink time-series:
            AllStoriesPostBlinks(s).Conditions(e).isPostBlink         = ...
                isPostBlinkEpochCorr;
            
            % Store a time vector in ms (the frame interval is the offset
            % of the first frame):
            AllStoriesPostBlinks(s).Conditions(e).times               = ...
                frameInterval : frameInterval : COND_DUR*1000;

        end % End of the loop over conditions
    end % End of the conditional
end % End of the loop over stories

%% Save and notify: -------------------------------------------------------

save(fullfile(PATHOUT_POST_BLINK_EPOCHS, 'post_speaker_blink_epochs.mat'), ...
    'AllStoriesPostBlinks');

disp('All of the speaker post-blink epochs were extracted successfully');

end % End of the main function



% -------------------------------------------------------------------------
%% ******************** Local functions ***********************************
% -------------------------------------------------------------------------

% Order of local functions:
% 1) local_match_blink_file_2_story_no()
% 2) local_check_num_conditions()
% 3) local_extract_post_blink_timestamps()



%% local_match_blink_file_2_story_no() ------------------------------------

% This function matches the index of a speaker blink file to the
% story-video which has the same name. As a result, the story number (in
% the experiment) the blink file corresponds to is determined.

function speakerBlinkFileStoryNo = local_match_blink_file_2_story_no( ...
                                             SpeakerBlinkFiles, storyNames)

[~,speakerBlinkFileName] = fileparts(SpeakerBlinkFiles.name);
speakerBlinkFileName     = extractBefore(speakerBlinkFileName,'_labels');
speakerBlinkFileStoryNo  = find(strcmp(speakerBlinkFileName, storyNames));
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


%% local_extract_post_blink_timestamps() ----------------------------------

% This function extracts a cell array containing the time-stamps of the
% post-blinks in samples. The time-stamps are relative to the start of each
% epoch. Each cell corresponds to one post-blink. Post-blinks too close to
% one another get eliminated. Also, post-blinks sliced between conditions
% get eliminated. The function also outputs a time-series of the corrected
% post-blinks.

function [postBlinkLatsCells, isPostBlinkEpochCorr] = ...
                    local_extract_post_blink_timestamps( ...
                               isPostBlinkEpoch, FRAME_RATE, MIN_GAP,...
                                                      POST_BLINK_SECS)

postBlinkLatsEpoch = find(isPostBlinkEpoch);

%%  Determine where the ends of the post-blinks are:

% A time-stamp more than a sample ahead of the previous one is part of a
% separate post-blink, hence:
isPostBlinkEnd   = 1 < (postBlinkLatsEpoch(2:end) - postBlinkLatsEpoch(1:end-1));

% The last stamp will always be the end of a post-blink.
isPostBlinkEnd   = [isPostBlinkEnd, true];

postBlinkEndLats = find(isPostBlinkEnd);

%% Segregate the post-blink time-stamps into different cells:
postBlinkLatsCells = {};
epoTill = 0;
for p = 1:length(postBlinkEndLats)
    epoFrom = epoTill + 1;
    epoTill = postBlinkEndLats(p);
    postBlinkLatsCells{p} = postBlinkLatsEpoch(epoFrom : epoTill);
end

%% Get rid of close post-blinks:

currIdx = 2;
numRemaining = length(postBlinkLatsCells);
while currIdx <= numRemaining - 1

    % Determine if the next post-blink is too close and eliminate accordingly:
    currentEndSecs = postBlinkLatsCells{currIdx}    (end) / FRAME_RATE;
    nextStartSecs  = postBlinkLatsCells{currIdx + 1}  (1) / FRAME_RATE;
    isNext2Close   = MIN_GAP > nextStartSecs - currentEndSecs;

    if isNext2Close
        postBlinkLatsCells(currIdx + 1) = [];
    end

    % Determine if the prior post-blink is too close and eliminate accordingly:
    priorEndSecs     = postBlinkLatsCells{currIdx - 1}(end) / FRAME_RATE;
    currentStartSecs = postBlinkLatsCells{currIdx}      (1) / FRAME_RATE;
    isPrior2Close    = MIN_GAP > currentStartSecs - priorEndSecs;

    if isPrior2Close
        postBlinkLatsCells(currIdx - 1) = [];
    end

    % Prepare to check the next post-blink:
    currIdx      =  currIdx + 1;
    numRemaining = length(postBlinkLatsCells);
end


%% Get rid of post-blinks sliced between conditions:
postBlinkDurationsSecs = cellfun(@length, postBlinkLatsCells) / FRAME_RATE; 
postBlinkLatsCells     = postBlinkLatsCells( ...
                            postBlinkDurationsSecs >= POST_BLINK_SECS);

% %% Get rid of overlapping post-blinks:
% currIdx = 1; 
% numRemaining = length(postBlinkLatsCells); 
% while currIdx <= numRemaining - 1
%     % Check whether the current post blink and the next one overlap:
%     areOverlapping = 0 < ...
%                         length(intersect(postBlinkLatsCells{currIdx}, ...
%                                          postBlinkLatsCells{currIdx+1}));
%     while areOverlapping 
%         % Delete the next post-blink 
%         postBlinkLatsCells(currIdx+1) = []; 
% 
%         % If there are no more cells left, terminate the loop:
%         numRemaining = length(postBlinkLatsCells);
%         if currIdx > numRemaining - 1
%             break;
%         end
% 
%         % Otherwise, check for further overlaps:
%         areOverlapping = 0 < ...
%                             length(intersect(postBlinkLatsCells{currIdx}, ...
%                                              postBlinkLatsCells{currIdx+1}));
%     end % End of the while loop eliminating overlaps
% 
%     currIdx = currIdx + 1;
% 
% end % End of the while loop iterating over post-blinks

%% Create a time-series containing only the remaining post-blinks:
isPostBlinkEpochCorr                          = false(size(isPostBlinkEpoch));
isPostBlinkEpochCorr([postBlinkLatsCells{:}]) = true;

end