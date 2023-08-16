% ps_01_extract_pauses_use_preset_thresholds( ...
%     FS_NEW, MIN_SPEECH, MIN_PAUSE, ...
%     PATHIN_AUDIO_FILES, PATHIN_THRESH_TBL, PATHOUT_PAUSES, )
%
% This function extracts pauses from the soundtracks. 
% 
% For each soundtrack the function first downsamples the sound data to the
% sample rate of the EEG data and applies transformations to it before
% extracting the pauses. The function saves a logical vector which
% specifies whether each sample was during a pause. The function also saves
% text files with the timestamps of the onsets and offsets of each pause.
% Finally, it saves a figure which display a histogram of the sound
% intensities, a plot of the time-series of the sound tracks, and a plot
% of the time-series of pauses.
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
% PATHIN_AUDIO_FILES - Where the raw audio files are stored as .wav files.
% PATHIN_THRESH_TBL -  Where the table with the decided sound thresholds
%                      is stored.
% PATHOUT_PAUSES     - Where the extracted pause vectors, Audacity labels, 
%                      and figures should be stored.



% -------------------------------------------------------------------------
%% ************ ps_01_extract_pauses_use_preset_thresholds() **************
% -------------------------------------------------------------------------

function ps_01_extract_pauses_use_preset_thresholds( ...
    FS_NEW, WIN_LEN_MOV_AVG, MIN_SPEECH, MIN_PAUSE, ...
    PATHIN_AUDIO_FILES, PATHIN_THRESH_TBL, PATHOUT_PAUSES)


% Load the table with the thresholds for each sound-track:
threshTbl  = readtable(fullfile(PATHIN_THRESH_TBL,'speech_pause_thresholds.xlsx'));

% Get a list of all the audio files:
AudioFiles = dir(fullfile(PATHIN_AUDIO_FILES,'*.wav'));


for snd = 1:length(AudioFiles)
    
    % Get the name of the file without the extension:
    [~,fileName] = fileparts(AudioFiles(snd).name);
    
    % Get sound intensity threshold for this soundtrack:
    idx          = find(strcmp(fileName, threshTbl.Soundtrack));
    sndThresh    = threshTbl.Threshold(idx);


    %% Process the audio data: ---------------------------------------------

    [originalSound, fsOld] = audioread( fullfile(PATHIN_AUDIO_FILES,AudioFiles(snd).name) );

    % Convert the sound to mono:
    soundDat = mean(originalSound,2);

    % Transform the sound:
    soundDat = abs(soundDat);
    soundDat = movmean(soundDat, round(WIN_LEN_MOV_AVG*fsOld));
    soundDat = log10(soundDat);
    soundDat = resample(soundDat, FS_NEW, fsOld);
    soundDat = soundDat';


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
    xline(histPlot, sndThresh, '--', 'Label','Intensity threshold');

    % Plot the sound intensity timeseries:
    soundPlot = subplot(2,3,2:3);
    timeVec   = (0:length(soundDat)-1) / FS_NEW;
    plot(timeVec, soundDat);
    soundAxes = gca;
    xlim([0,50]);
    ylim([-5,0]);
    xlabel('Time [s]');
    ylabel('Absolute sound intensity [dB]');
    title('Sound intensity over time');

    % Draw the amplitude threshold on the sound intensity timeseries:
    yline(soundPlot, sndThresh, '--', 'Label', 'Intensity threshold');

    % Visualize the pauses:
    subplot(2,3,5:6);
    [onsLats, ofsLats, isPause] = local_extract_pauses(soundDat, sndThresh, ...
                                             MIN_PAUSE, MIN_SPEECH, FS_NEW);
    stem(timeVec, isPause);
    pauseAxes  = gca;
    xlim( soundAxes.XLim );
    ylim([0,1]);
    xlabel('Time [s]');
    ylabel('Is it a pause? [bool]');
    pauseAxes.YAxis.TickValues = [0,1];
    title(['Pauses with thresh = ', num2str(sndThresh), ' DB']); 


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

    % Save the figure:
    saveas(fig, fullfile(PATHOUT_PAUSES, [fileName,'_pauses.jpg']));
    close(fig);

end % End of the for loop

end % End of the main function



% -------------------------------------------------------------------------
%% ************************* Local functions ******************************
% -------------------------------------------------------------------------

% Order of local functions:
% 1) local_extract_pauses()
% 2) local_correct_extremeties()



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


