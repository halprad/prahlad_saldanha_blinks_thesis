% ps_01_extract_pauses_change_thresholds( ...
%     FS_NEW, MIN_SPEECH, MIN_PAUSE, DEFAULT_THRESH, ...
%     PATHIN_AUDIO_FILES, PATHOUT_PAUSES, PATHOUT_THRESH_TBL)
%
% Adjust the sound thresholds and extract pauses.
%
% Inputs:
% FS_NEW             - Sample rate the sound data should be downsampled to.
%                      This should be set to be the same sample rate as
%                      that of the EEG data.
% WIN_LEN_MOV_AVG    - The length of a moving average filter window in
%                      seconds which will be used to smoothen the sound
%                      data.
% MIN_SPEECH         - The minimum acceptable duration of a speech segment
%                      in seconds.
% MIN_PAUSE          - The minimum acceptable duration of a pause segment 
%                      in seconds.
% DEFAULT_THRESH     - The default cut-off threshold in dB to extract
%                      pauses.
% PATHIN_AUDIO_FILES - Where the raw audio files are stored as .wav files.
% PATHOUT_PAUSES     - Where the extracted pause vectors, Audacity labels, 
%                      and figures should be stored.
% PATHOUT_THRESH_TBL - Where the table with the decided sound thresholds
%                      should be stored.
%
% An interactive figure will pop-up for every audio file being opened. The
% subplot on the left will display a histogram of absolute sound
% intensities. The subplot on the top right will display a snippet of the
% absolute sound intensity time-series. The subplot on the bottom right
% will display whether regions of that snippet were classified as pauses or
% not. After deciding the threshold for each audio file, the function will
% save the extracted pauses, the Audacity labels of the speech onsets, the
% decided thresholds, and the figures.
%
% Interactivity:
% <- and ->  arrowkeys  - Scrolls through the time-series.
% spacebar              - Plays the audio of the displayed snippet.
% backspace             - Creates a cross-hair on the mouse's cursor. Click
%                         within the histogram to adjust the sound
%                         threshold.
% enter/return          - Accept the current threshold and save the figure
%                         and pauses.



% -------------------------------------------------------------------------
%% ************* ps_01_extract_pauses_change_thresholds() *****************
% -------------------------------------------------------------------------

function ps_01_extract_pauses_change_thresholds( ...
    FS_NEW, WIN_LEN_MOV_AVG, MIN_SPEECH, MIN_PAUSE, DEFAULT_THRESH, ...
    PATHIN_AUDIO_FILES, PATHOUT_PAUSES, PATHOUT_THRESH_TBL)


% Get a list of all the audio files:
AudioFiles = dir(fullfile(PATHIN_AUDIO_FILES,'*.wav'));

% Preallocate space for the selected thresholds:
selectedThreshs = cell(length(AudioFiles), 2);

for snd = 1:length(AudioFiles)
    
    % Get the name of the file without the extension:
    [~,fileName] = fileparts(AudioFiles(snd).name);


    %% Process the audio data: --------------------------------------------

    [originalSound, fsOld] = audioread( fullfile(PATHIN_AUDIO_FILES,AudioFiles(snd).name) );

    % Convert the sound to mono:
    soundDat = mean(originalSound,2);

    % Transform the sound:
    soundDat = abs(soundDat);
    soundDat = movmean(soundDat, round(WIN_LEN_MOV_AVG*fsOld));
    soundDat = log10(soundDat);
    soundDat = resample(soundDat, FS_NEW, fsOld);
    soundDat = soundDat';

    % Determine the length of the data in seconds:
    datLengthSecs = length(soundDat) / FS_NEW;


    %% Visualize the audio data: ------------------------------------------

    % Create a figure:
    fig = figure;
    fig.WindowState = 'maximized';
    sgtitle(fileName, 'Interpreter', 'none');

    % Draw a histogram of the sound amplitudes:
    histPlot = subplot(2,3,[1,4]);
    histogram(soundDat);
    xlim([1.1*min(soundDat), 1.1*max(soundDat)]);
    xlabel('Absolute sound intensity [dB]');
    ylabel('Number of samples');
    title('Distribution of sound intensities');

    % Draw the amplitude threshold on the histogram:
    threshLineHist = xline(histPlot, DEFAULT_THRESH, '--', 'Label','Intensity threshold');

    % Plot the sound intensity timeseries:
    soundPlot = subplot(2,3,2:3);
    timeVec   = (0:length(soundDat)-1) / FS_NEW;
    plot(timeVec, soundDat);
    soundAxes = gca;
    xlim([0,5]);
    ylim([-5,0]);
    xlabel('Time [s]');
    ylabel('Absolute sound intensity [dB]');
    title('Sound intensity over time');

    % Draw the amplitude threshold on the sound intensity timeseries:
    threshLineSound = yline(soundPlot, DEFAULT_THRESH, '--', 'Label', 'Intensity threshold');

    % Visualize the pauses:
    subplot(2,3,5:6);
    [onsLats, ofsLats, isPause] = local_extract_pauses(soundDat, DEFAULT_THRESH, ...
                                             MIN_PAUSE, MIN_SPEECH, FS_NEW);
    pausePlot  = stem(timeVec, isPause);
    pauseAxes  = gca;
    xlim( soundAxes.XLim );
    ylim([0,1]);
    xlabel('Time [s]');
    ylabel('Is it a pause? [bool]');
    pauseAxes.YAxis.TickValues = [0,1];
    pauseTitle = title(['Pauses with thresh = ', num2str(DEFAULT_THRESH), ' DB']);

    
    %% Adjust the sound intensity threshold: -------------------------------

    % Make the figure listen for button presses:
    fig.WindowKeyPressFcn = @local_listen_4_keys;
    
    % Allow the user to adjust the sound threshold:
    sndThresh         = DEFAULT_THRESH;
    checkAdjustThresh = true;

    while checkAdjustThresh

        % Listen for button presses:
        fig.UserData.Key = '';
        waitfor(fig, 'UserData');

        % Adjust the sound threshold accordingly:
        switch fig.UserData.Key
            case 'rightarrow'
                local_scroll_right(soundAxes,pauseAxes,datLengthSecs)
            case 'leftarrow'
                local_scroll_left(soundAxes,pauseAxes)
            case 'space'
                local_play_sound(soundAxes,pauseAxes,originalSound,fsOld)
            case 'backspace'
                [onsLats, ofsLats, isPause, sndThresh] = ...
                local_adjust_thresh( soundDat, MIN_PAUSE, MIN_SPEECH, FS_NEW, ...
                                   threshLineHist, threshLineSound, pauseTitle, pausePlot);
            case 'return'
                checkAdjustThresh = false;
        end

    end % End of the while loop

    %% Save the pause data: -----------------------------------------------

    % Save the pause vector:
    save( fullfile(PATHOUT_PAUSES,[fileName,'_pauses.mat']), 'isPause');

    % Save the onset and offset timepoints (in seconds):
    onsSecs = (onsLats/FS_NEW)';
    ofsSecs = (ofsLats/FS_NEW)';

    onsLabs = repmat({'Onset'},  length(onsSecs), 1);
    ofsLabs = repmat({'Offset'}, length(ofsSecs), 1);

    onsTbl  = table(onsSecs, onsSecs, onsLabs);
    writetable(onsTbl, ...
               fullfile(PATHOUT_PAUSES, [fileName,'_onsets.txt']), ...
               'WriteVariableNames', false, ...
               'Delimiter', '\t');

    ofsTbl  = table(ofsSecs, ofsSecs, ofsLabs);
    writetable(ofsTbl, ...
               fullfile(PATHOUT_PAUSES, [fileName,'_offsets.txt']), ...
               'WriteVariableNames', false, ...
               'Delimiter', '\t');

    % Store the decided sound threshold:
    selectedThreshs(snd, :) = {fileName, sndThresh};

    % Save the figure:
    soundAxes.XLim = [0, 50];
    pauseAxes.XLim = [0, 50];
    saveas(fig, fullfile(PATHOUT_PAUSES, [fileName,'_pauses.jpg']));
    close(fig);

end % End of the for loop

%% Save the thresholds: ---------------------------------------------------

threshTbl = cell2table(selectedThreshs, 'VariableNames', {'Soundtrack','Threshold'});
writetable(threshTbl, fullfile(PATHOUT_THRESH_TBL,'speech_pause_thresholds.xlsx'));

end % End of the main function



% -------------------------------------------------------------------------
%% ************************* Local functions ******************************
% -------------------------------------------------------------------------

% Order of local functions:
% 1) local_extract_pauses()
% 2) local_correct_extremeties()
% 3) local_listen_4_keys()
% 4) local_scroll_right()
% 5) local_scroll_left()
% 6) local_play_sound()
% 7) local_adjust_thresh()


%% local_extract_pauses() -------------------------------------------------

% This function finds the onset and offset latencies of speech segments,
% and creates a logical vector which specifies whether each sample in
% soundDat is part of a pause or not.

function  [onsLats, ofsLats, isPause] = local_extract_pauses(soundDat, sndThresh, ...
                                                   MIN_PAUSE, MIN_SPEECH, FS_NEW)

% Find initial speech-edges
isLoud   = soundDat >= sndThresh;
isOnset  = ~isLoud(1:end-1) & isLoud(2:end);
isOnset  = [isOnset, false]; % So that isOnset is the same length as soundDat.
isOffset = isLoud(1:end-1) & ~isLoud(2:end);
isOffset = [isOffset, false]; 

% The preceeding code doesn't detect edges at the boundaries of soundDat
% well, but these are needed to extract pauses, hence:
[isOnset, isOffset] = local_correct_extremeties(isOnset, isOffset);

% Convert short speech segments to pauses:
onsLats = find(isOnset);
ofsLats = find(isOffset);
for ons = onsLats
    laterOfs = ofsLats(ofsLats >= ons);

    if isempty(laterOfs)
        break;
    end

    if laterOfs(1)-ons < MIN_SPEECH*FS_NEW
        isOnset(ons) = false;
        isOffset(laterOfs(1)) = false;
    end
end

% The extremeties might need corrections again, hence:
[isOnset, isOffset] = local_correct_extremeties(isOnset, isOffset);

% Convert short pauses to speech:
onsLats = find(isOnset);
ofsLats = find(isOffset);
for ofs = ofsLats
    laterOns = onsLats(onsLats >= ofs);

    if isempty(laterOns)
        break;
    end

    if laterOns(1)-ofs < MIN_PAUSE*FS_NEW
        isOffset(ofs) = false;
        isOnset(laterOns(1)) = false;
    end
end

% Again, the extremeties might need correcting:
[isOnset, isOffset] = local_correct_extremeties(isOnset, isOffset);

% Create timestamps in samples of the speech onsets and offsets, and create
% a logical vector which answers whether each sample in soundDat is part of
% a pause or not:
onsLats = find(isOnset);
ofsLats = find(isOffset);
isPause = false(size(isOnset));
for ofs = ofsLats
    laterOns = onsLats(onsLats >= ofs);

    if isempty(laterOns)
        break;
    end

    isPause(ofs:laterOns(1)) = true;
end

end


%% local_correct_extremeties() --------------------------------------------

% This function determines whether the first and last samples of soundDat
% should be classified as speech onset or offset edges to ensure that
% speech and pause segments at the tips of soundDat can be classified
% properly.

function  [isOnset, isOffset] = local_correct_extremeties(isOnset, isOffset)

% Determine whether the first sample should be a speech onset or an offset
% edge:
isFirstOnsetBeforeOffset = find(isOnset,1,"first") <= find(isOffset,1,"first");
if isFirstOnsetBeforeOffset && ~isOnset(1)
    isOffset(1) = true;
elseif ~isOffset(1)  
    isOnset(1) = true;
end

% Determine whether the last sample should be an onset or offset edge:
isLastOffsetAfterOnset = find(isOffset,1,"last") >= find(isOnset,1,"last");
if isLastOffsetAfterOnset && ~isOffset(end)
    isOnset(end) = true;
elseif ~isOnset(end)
    isOffset(end) = true;
end

end


%% local_listen_4_keys() --------------------------------------------------

% This function detects which key the user pressed while the figure was
% open and stores the key's name in the figure object:

function  local_listen_4_keys(fig,KeyDataObj)
fig.UserData.Key = KeyDataObj.Key;  
end


%% local_scroll_right() ---------------------------------------------------

% This function changes the sound and pause axes limits to show later
% timepoints:

function local_scroll_right(soundAxes, pauseAxes, datLengthSecs)
currLim = soundAxes.XLim;
newLim  = currLim + 0.5;
if any(newLim > datLengthSecs)
    return;
end
soundAxes.XLim = newLim;
pauseAxes.XLim = newLim;
end


%% local_scroll_left() ----------------------------------------------------

% This function changes the sound and pause axes limits to show earlier
% timepoints:

function local_scroll_left(soundAxes,pauseAxes)
currLim = soundAxes.XLim;
newLim  = currLim - 0.5;
if any(newLim < 0)
    return;
end
soundAxes.XLim = newLim;
pauseAxes.XLim = newLim;
end


%% local_play_sound() -----------------------------------------------------

% This function plays the sound of the segment in timeseries which is
% currently being displayed. It moves a red line across the timeseries in
% synch with the sound.

function local_play_sound(soundAxes, pauseAxes, originalSound,fsOld)
audIdx     = 1 + (soundAxes.XLim*fsOld);
audio      = originalSound(audIdx(1):audIdx(2), :);
startTime  = tic;
sound(audio, fsOld);
soundLine  = xline(soundAxes, soundAxes.XLim(1), 'r', 'LineWidth', 2);
pauseLine  = xline(pauseAxes, pauseAxes.XLim(1), 'r', 'LineWidth', 2);
dispSegDur = 5; % Duration of the sound data segment on display in seconds
while toc(startTime) < dispSegDur - 0.1
    linePosistion   = toc(startTime) + soundAxes.XLim(1);
    soundLine.Value = linePosistion;
    pauseLine.Value = linePosistion;
    drawnow limitrate;
end
delete(soundLine);
delete(pauseLine);
end


%% local_adjust_thresh() --------------------------------------------------

% This function allows the user to click on the histogram to adjust the
% pause threshold. It updates the plots in the figure to reflect the
% adjustment.

function [onsLats, ofsLats, isPause, ampThresh] = ...
         local_adjust_thresh( soundDat, MIN_PAUSE, MIN_SPEECH, FS_NEW, ...
                  threshLineHist, threshLineSound, pauseTitle, pausePlot)

[ampThresh,~] = ginput(1);
ampThresh     = round(ampThresh, 2); % Round to two decimal places

% Find the pauses with the new threshold:
[onsLats, ofsLats, isPause] = local_extract_pauses(soundDat, ampThresh, ...
                                              MIN_PAUSE, MIN_SPEECH, FS_NEW);

% Update the figure:
pausePlot.YData       = isPause;
pauseTitle.String     = ['Pauses with thresh = ', num2str(ampThresh), ' DB'];
threshLineHist.Value  = ampThresh;
threshLineSound.Value = ampThresh;
end

