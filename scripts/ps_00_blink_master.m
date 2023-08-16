% ps_00_blink_master.m
%
% Author  - Prahlad Saldanha 
%
% Credits - Much of these scripts were based on code from Holtze et al.
%           2023, which is available at
%           https://doi.org/10.5281/ZENODO.7612544. The video labeler was
%           based on an app developed by Ibrahim Celik. The ICA
%           decomposition was based on code from Stropahl et al. 2018,
%           which is available at
%           https://figshare.com/s/48f8d9de715bafa5811b
%
% This script will run all of the code required to reproduce the results
% for my master's thesis. It uses separate function scripts for each step
% of the process. Many of these functions use local functions.
%
% This script is structured as follows:
%
% Part 0: Preliminaries
%
% Part 1: Extract variables
%
% Part 2: Extract statistics
%
% Part 3: Plot results
%
% The naming conventions are as follows:
%
% Main functions  -> They are all prefixed with 'ps' and are named with
%                    lower case letters in snake_case.
% Local functions -> They are all prefixed with 'local' and are named with
%                    lower case letters in snake_case.
% Paths           -> They are in ALLCAPS SNAKE_CASE with PATH, PATHIN
%                    or PATHOUT as prefixes.
% Function inputs -> The inputs to all main functions are in ALLCAPS
%                    SNAKE_CASE. 
% Variables       -> Other variables are all in lowerCamelCase except
%                    for structures - they are in UpperCamelCase.



% -------------------------------------------------------------------------
%% ****************** Part 0: Preliminaries *******************************
% -------------------------------------------------------------------------

close all; clear all; clc;


% Get the paths to all the folders needed for the analysis: ---------------

PATH_RAW        = uigetdir('','Select the raw data folder');

PATH_PROCESSED  = uigetdir('','Select the processed data folder');

PATH_STIMULI    = uigetdir('','Select the folder with the stimuli');

PATH_TOOLBOX    = uigetdir('','Select the folder with the toolboxes');


% Add the toolboxes to MATLAB's search path: ------------------------------

% EEGLAB version 2023.0 (modify this line if using a different version):
addpath(fullfile(PATH_TOOLBOX,'eeglab2023.0/'));

% Blinker version 1.1.3 (modify this line if using a different version):
addpath(genpath(fullfile(PATH_TOOLBOX,'eeglab2023.0/plugins/Blinker1.1.3/')));

% EyeCatch:
addpath(fullfile(PATH_TOOLBOX,'eye_catch/'));


% Load the stimulus-key table: --------------------------------------------

STIM_KEY_TBL = readtable( fullfile(PATH_STIMULI,'stimulus_key.xlsx'), ...
                          'NumHeaderLines', 1, 'ReadRowNames', true);


% Get the names of the variables in the workspace: ------------------------

% (This is for tidyness and modularity later)
keepVars = who; 
keepVars = [keepVars; {'keepVars'}];


% -------------------------------------------------------------------------
%% ****************** Part 1: Extract variables ***************************
% -------------------------------------------------------------------------

%% 1-a Extract pauses: ----------------------------------------------------

% The pauses can either be extracted by using a table with the pre-existing
% sound intensity thresholds for each of the soundtracks, or by manually
% changing the thresholds with the help of interactive figures.

% Function inputs (read the function documentations for details):
FS_NEW          = 500;   % In Hz
WIN_LEN_MOV_AVG = 0.01;  % In seconds
MIN_SPEECH      = 0.015; % In seconds
MIN_PAUSE       = 0.2;   % In seconds
DEFAULT_THRESH  = -3;    % In dB (only useful when manually adjusting 
                         % thresholds)

% Relevant paths:
PATHIN_AUDIO_FILES = fullfile(PATH_STIMULI, 'sound_extracts');
PATHOUT_PAUSES     = fullfile(PATH_PROCESSED, '01_speech_pauses');
PATHIN_THRESH_TBL  = fullfile(PATH_PROCESSED, '00_speech_pause_thresholds');
PATHOUT_THRESH_TBL = fullfile(PATH_PROCESSED, '00_speech_pause_thresholds_adjusted');

if ~exist(PATHOUT_PAUSES, "dir")
    mkdir(PATHOUT_PAUSES);
end

if ~exist(PATHOUT_THRESH_TBL, "dir")
    mkdir(PATHIN_THRESH_TBL);
end

% Use the preset sound intensity thresholds to extract pauses? 
usePresetThresholds = true; % "true" will reproduce the reported results

% Run the main function:

if usePresetThresholds

    ps_01_extract_pauses_use_preset_thresholds( ...
        FS_NEW, WIN_LEN_MOV_AVG, MIN_SPEECH, MIN_PAUSE, ...
        PATHIN_AUDIO_FILES, PATHIN_THRESH_TBL, PATHOUT_PAUSES);

else    

    ps_01_extract_pauses_change_thresholds( ...
        FS_NEW, WIN_LEN_MOV_AVG, MIN_SPEECH, MIN_PAUSE, DEFAULT_THRESH, ...
        PATHIN_AUDIO_FILES, PATHOUT_PAUSES, PATHOUT_THRESH_TBL);

end

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});

%% 1-b Extract speaker blinks: --------------------------------------------

% This was done manually by using 'ps_video_labeler.mlapp'. 
% 
% The labels can be inspected by opening one of the raw videos from the
% stimulus folder and loading the labels from its corresponding file in
% processed data -> 00_speaker_blink_frames.


%% 1-c Extract listener blinks: -------------------------------------------

%% 1-c-i  Compute the ICA weights:

% Settings for the EEG data preparation before the ICA decomposition:
LOW_CUT_OFF  = 1;  % High-pass filter cut-off in Hz.
HIGH_CUT_OFF = 20; % Low-pass filter cut-off in Hz.
PRUNE        = 3;  % Threshold in standard deviations for rejecting regular
                   % epochs based on joint-probability.

% Relevant paths (one renamed to clarify its purpose):
PATH_TOOLBOX;
PATHIN_RAW   = PATH_RAW;
PATHOUT_ICA  = fullfile(PATH_PROCESSED, '02_eeg_sessions_ica');

if ~exist(PATHOUT_ICA, "dir")
    mkdir(PATHOUT_ICA);
end

% Run the main function:
ps_02_do_ica( ...
    LOW_CUT_OFF, HIGH_CUT_OFF, PRUNE, ...
    PATH_TOOLBOX, PATHIN_RAW, PATHOUT_ICA);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});


%% 1-c-ii Extract the listener blinks with the Blinker toolbox:

% The Blinker toolbox is made to extract eye-blinks by using the IC
% activations of the merged-by-session EEG datasets.

% Relevant paths:
PATHIN_ICA              = fullfile(PATH_PROCESSED, '02_eeg_sessions_ica/');
PATHOUT_LISTENER_BLINKS = fullfile(PATH_PROCESSED, '03_listener_blinks');

if ~exist(PATHOUT_LISTENER_BLINKS, "dir")
    mkdir(PATHOUT_LISTENER_BLINKS);
end

% Run the main function:
ps_03_extract_listener_blinks( ...
    PATHIN_ICA, PATHOUT_LISTENER_BLINKS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});

%% 1-d Epoch pauses: ------------------------------------------------------

% The pause get epoched into conditions for each story. The epochs get
% stored in a strucure array.

% Information to epoch pauses:
SRATE     = 500;  % In Hz - the EEG sample rate.
COND_DUR  = 30;   % In secs - determines the length of epochs.
MIN_GAP   = 1;    % In secs - pauses closer together than this will be 
                  % eliminated.
MIN_PAUSE = 0.2;  % In secs - pauses shorter than this will be eliminated.
MAX_PAUSE = 1;    % In secs - pauses longer than this will be eliminated.
STIM_KEY_TBL;

% Relevant paths:
PATHIN_PAUSES         = fullfile(PATH_PROCESSED, '01_speech_pauses/');
PATHOUT_PAUSE_EPOCHS  = fullfile(PATH_PROCESSED, '04_pause_epochs');

if ~exist(PATHOUT_PAUSE_EPOCHS, 'dir')
    mkdir(PATHOUT_PAUSE_EPOCHS);
end

% Run the main function:
ps_04_epoch_pauses(...
    SRATE, COND_DUR, MIN_GAP, MIN_PAUSE, MAX_PAUSE, STIM_KEY_TBL, ...
    PATHIN_PAUSES, PATHOUT_PAUSE_EPOCHS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});


%% 1-e Epoch post-speaker-blink intervals: --------------------------------

% The intervals following speaker-blinks get epoched into conditions for
% each story which was presented during the experiment. The epochs get
% stored within a structure array.

% Information to epoch post-speaker-blinks:
FRAME_RATE      = 25; % Number of video per second
COND_DUR        = 30; % Number of seconds in each of the conditions
MIN_GAP         = 1;  % Number of seconds that should be between 
                      % post-blinks.
POST_BLINK_SECS = 1;  % Number of seconds after a blink to consider the 
                      % post-blink
STIM_KEY_TBL;

% Relevant paths:
PATHIN_SPEAKER_BLINKS              = ...
            fullfile(PATH_PROCESSED, '00_speaker_blink_frames/');
PATHOUT_POST_BLINK_EPOCHS  = ...
            fullfile(PATH_PROCESSED, '05_post_speaker_blink_epochs');

if ~exist(PATHOUT_POST_BLINK_EPOCHS, "dir")
    mkdir(PATHOUT_POST_BLINK_EPOCHS);
end

% Run the main function:
ps_05_epoch_post_speaker_blinks( ...
    FRAME_RATE, COND_DUR, MIN_GAP, POST_BLINK_SECS, STIM_KEY_TBL, ...
    PATHIN_SPEAKER_BLINKS, PATHOUT_POST_BLINK_EPOCHS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});

%% 1-f Epoch listener blinks: ---------------------------------------------

% The listener-blink time-stamps get epoched into conditions for each story
% which was presented during the experiment. The epochs get stored within a
% structure array - one for every participant.

% Information to epoch listener-blinks:
COND_DUR = 30; % In secs
SRATE    = 500;% In hertz - EEG sample rate

% Relevant paths:
PATHIN_LISTENER_BLINKS          = fullfile(PATH_PROCESSED, '03_listener_blinks');
PATHOUT_LISTENER_BLINK_EPOCHS   = fullfile(PATH_PROCESSED, '06_listener_blink_epochs');

if ~exist(PATHOUT_LISTENER_BLINK_EPOCHS, "dir")
    mkdir(PATHOUT_LISTENER_BLINK_EPOCHS);
end

% Run the main function:
ps_06_epoch_listener_blinks(...
    COND_DUR, SRATE, ...
    PATHIN_LISTENER_BLINKS, PATHOUT_LISTENER_BLINK_EPOCHS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});


% -------------------------------------------------------------------------
%% ****************** Part 2: Extract statistics **************************
% -------------------------------------------------------------------------

%% 2-a Calculate proportions of pauses containing listener blinks: --------

% The proportions of pauses that contain listener blinks get calculated for
% each type of condition. This produces two tables, one of which can be
% imported into R for the regressions.

% Information to calculate separate proportions for each type of condition:
CONDITION_TYPES = {'AV', 'V', 'A', 'AV-nolips'};
STIM_KEY_TBL;

% Relevant paths:
PATHIN_LISTENER_BLINK_EPOCHS = fullfile(PATH_PROCESSED, '06_listener_blink_epochs');
PATHIN_PAUSE_EPOCHS          = fullfile(PATH_PROCESSED, '04_pause_epochs/');
PATHOUT_PROPORTIONS          = fullfile(PATH_PROCESSED, '07_listener_blink_pause_proportions');

if ~exist(PATHOUT_PROPORTIONS, "dir")
    mkdir(PATHOUT_PROPORTIONS);
end

% Run the main function:
ps_07_calculate_listener_blink_pause_proportions(...
    CONDITION_TYPES, STIM_KEY_TBL, ...
    PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_PAUSE_EPOCHS, PATHOUT_PROPORTIONS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});


%% 2-b Shuffle pauses and calculate listener blink-pause proportions: -----

% The pauses are shuffled several times and each time the proportion of one
% type of condition containing at least one speaker blink is calculated.
% Vectors of these proportions are calculated for each participant.

% Information for the shuffling and calculating of proportions:
CONDITION_TYPE = 'AV';  % Type of condition for which to calculate proportions.
N_ITERATIONS   = 10000; % Number of times to shuffle.
STIM_KEY_TBL;

% Relevant paths:
PATHIN_LISTENER_BLINK_EPOCHS = fullfile(PATH_PROCESSED, '06_listener_blink_epochs/');
PATHIN_PAUSE_EPOCHS          = fullfile(PATH_PROCESSED, '04_pause_epochs/');
PATHOUT_SHUFFLED_PROPORTIONS = fullfile(PATH_PROCESSED, '08_shuffled_listener_blink_pause_proportions');

if ~exist(PATHOUT_SHUFFLED_PROPORTIONS, "dir")
    mkdir(PATHOUT_SHUFFLED_PROPORTIONS);
end

% Run the main function:
ps_08_shuffle_pauses_calculate_listener_blink_pause_proportions( ...
    CONDITION_TYPE, STIM_KEY_TBL, N_ITERATIONS, ...
    PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_PAUSE_EPOCHS, PATHOUT_SHUFFLED_PROPORTIONS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});


%% 2-c Calculate proportions of post-blinks containing listener blinks: ---

% Information to calculate separate proportions for each type of condition:
FRAME_RATE      = 25;  % Video frame rate in Hz
SRATE           = 500; % EEG sample rate in Hz
CONDITION_TYPES = {'AV', 'V', 'A', 'AV-nolips'};
STIM_KEY_TBL; 

% Relevant paths:
PATHIN_LISTENER_BLINK_EPOCHS = fullfile(PATH_PROCESSED, '06_listener_blink_epochs');
PATHIN_POST_BLINK_EPOCHS     = fullfile(PATH_PROCESSED, '05_post_speaker_blink_epochs');
PATHOUT_PROPORTIONS          = fullfile(PATH_PROCESSED, '09_synched_blink_proportions');

if ~exist(PATHOUT_PROPORTIONS, "dir")
    mkdir(PATHOUT_PROPORTIONS);
end

% Run the main function:
ps_09_calculate_synched_blink_proportions( ...
    FRAME_RATE, SRATE, CONDITION_TYPES, STIM_KEY_TBL, ...
    PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_POST_BLINK_EPOCHS, ...
    PATHOUT_PROPORTIONS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});


%% 2-d Shuffle post-blinks and calculate synched blink proportions: -------

FRAME_RATE = 25;
SRATE = 500; 
CONDITION_TYPE = 'AV';
STIM_KEY_TBL;
N_ITERATIONS = 10000;


PATHIN_LISTENER_BLINK_EPOCHS = fullfile(PATH_PROCESSED, '06_listener_blink_epochs/');
PATHIN_POST_BLINK_EPOCHS     = fullfile(PATH_PROCESSED, '05_post_speaker_blink_epochs/');
PATHOUT_SHUFFLED_PROPORTIONS = fullfile(PATH_PROCESSED, '10_shuffled_synched_blink_proportions');

if ~exist(PATHOUT_SHUFFLED_PROPORTIONS, "dir")
    mkdir(PATHOUT_SHUFFLED_PROPORTIONS);
end

ps_10_shuffle_post_blinks_calculate_synched_blink_proportions(...
    FRAME_RATE, SRATE, CONDITION_TYPE, STIM_KEY_TBL, N_ITERATIONS, ...
    PATHIN_LISTENER_BLINK_EPOCHS, PATHIN_POST_BLINK_EPOCHS, ...
    PATHOUT_SHUFFLED_PROPORTIONS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});


%% 2-e Calculate proportion of pauses containing speaker blinks:-----------

SPEAKERS = {'Aaron', 'Anna', 'Janto', 'Jupiter', 'Laura', 'Mara'};
FRAME_RATE = 25;
SRATE = 500;
COND_DUR = 30;
MIN_GAP = 1;
MIN_PAUSE = 0.2;
MAX_PAUSE = 1;

PATHIN_SPEAKER_BLINKS = fullfile(PATH_PROCESSED, '00_speaker_blink_frames/');
PATHIN_PAUSES = fullfile(PATH_PROCESSED,'01_speech_pauses');
PATHOUT_PROPORTIONS = fullfile(PATH_PROCESSED, '11_speaker_blink_pause_proportions');

if ~exist(PATHOUT_PROPORTIONS, "dir")
    mkdir(PATHOUT_PROPORTIONS);
end

ps_11_calculate_speaker_blink_pause_proportions( ...
    SPEAKERS, STIM_KEY_TBL, FRAME_RATE, SRATE, COND_DUR, ...
    MIN_GAP, MIN_PAUSE, MAX_PAUSE, ...
    PATHIN_SPEAKER_BLINKS, PATHIN_PAUSES, ...
    PATHOUT_PROPORTIONS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});


%% 2-f Shuffle pauses and calculate speaker blink-pause proportions: ------

SPEAKERS = {'Aaron', 'Anna', 'Janto', 'Jupiter', 'Laura', 'Mara'};
N_ITERATIONS = 10000;
FRAME_RATE = 25;
SRATE = 500;
COND_DUR = 30;
MIN_GAP = 1;
MIN_PAUSE = 0.2;
MAX_PAUSE = 1;

PATHIN_SPEAKER_BLINKS = fullfile(PATH_PROCESSED, '00_speaker_blink_frames/');
PATHIN_PAUSES = fullfile(PATH_PROCESSED,'01_speech_pauses');
PATHOUT_SHUFFLED_PROPORTIONS = fullfile(PATH_PROCESSED, '12_shuffled_speaker_blink_pause_proportions');

if ~exist(PATHOUT_SHUFFLED_PROPORTIONS, "dir")
    mkdir(PATHOUT_SHUFFLED_PROPORTIONS);
end

ps_12_shuffle_pauses_calculate_speaker_blink_pause_proportions( ...
    SPEAKERS, STIM_KEY_TBL, N_ITERATIONS, FRAME_RATE, SRATE, COND_DUR, ...
    MIN_GAP, MIN_PAUSE, MAX_PAUSE, ...
    PATHIN_SPEAKER_BLINKS, PATHIN_PAUSES, ...
    PATHOUT_SHUFFLED_PROPORTIONS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});

%% 2-g Do binomial tests: -------------------------------------------------

% Check whether the observed proportions of the listener blinks in pauses
% are significant each of the listeners. Determine whether the number of
% significant cases is significant using a binomial test. Then repeat the
% process to determine whether the number listener whose blinks following
% speaker blinks was significant, and whether the number of speakers who
% blinked in speech-pauses was significant.

CONDITION               = 'AV';
THRESH_INDIVIDUAL_SUBS  = [0.05, 0.05, 0.05];
THRESH_BINOMIAL_TESTS   = [0.05, 0.05, 0.05];

PATHIN_LISTENER_PROPORTIONS                 = fullfile(PATH_PROCESSED, '07_listener_blink_pause_proportions');
PATHIN_SHUFFLED_LISTENER_PROPORTIONS        = fullfile(PATH_PROCESSED, '08_shuffled_listener_blink_pause_proportions');
PATHIN_SYNCHED_BLINKS_PROPROTIONS           = fullfile(PATH_PROCESSED, '09_synched_blink_proportions');
PATHIN_SHUFFLED_SYNCHED_BLINKS_PROPORTIONS  = fullfile(PATH_PROCESSED, '10_shuffled_synched_blink_proportions');
PATHIN_SPEAKER_PROPORTIONS                  = fullfile(PATH_PROCESSED, '11_speaker_blink_pause_proportions');
PATHIN_SHUFFLED_SPEAKER_PROPORTIONS         = fullfile(PATH_PROCESSED, '12_shuffled_speaker_blink_pause_proportions');
PATHOUT_BINOMIAL_TEST_RESULTS               = fullfile(PATH_PROCESSED, '13_binomial_test_results');

if ~exist(PATHOUT_BINOMIAL_TEST_RESULTS, 'dir')
    mkdir(PATHOUT_BINOMIAL_TEST_RESULTS);
end

ps_13_do_binomial_tests(...
    CONDITION, THRESH_INDIVIDUAL_SUBS, THRESH_BINOMIAL_TESTS, ...
    PATHIN_LISTENER_PROPORTIONS, PATHIN_SHUFFLED_LISTENER_PROPORTIONS, ...
    PATHIN_SYNCHED_BLINKS_PROPROTIONS, PATHIN_SHUFFLED_SYNCHED_BLINKS_PROPORTIONS, ...
    PATHIN_SPEAKER_PROPORTIONS, PATHIN_SHUFFLED_SPEAKER_PROPORTIONS, ...
    PATHOUT_BINOMIAL_TEST_RESULTS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});

% -------------------------------------------------------------------------
%% ****************** Part 3: Plot results ********************************
% -------------------------------------------------------------------------

%% 3-a Plot comparisons between conditions of the observed proportions: ---

% Info for the function:
COND_CONTRASTS_BLINK_PAUSES   = {'V', 'AV'; ...
                                 'A', 'AV'; ...
                                 'AV-nolips', 'AV'};
COND_CONTRASTS_SYNCHED_BLINKS = {'V', 'AV'; ...
                                 'A', 'AV'};
HALF_A4_SIZE                  = [21-3.18, (29.7-2.54)/2]; % In cm

% Relevant paths:
LISTENER_BLINK_PAUSE_PROPORTIONS_PATH = fullfile(PATH_PROCESSED,'07_listener_blink_pause_proportions/');
SYNCHED_BLINK_PROPORTIONS_PATH        = fullfile(PATH_PROCESSED,'09_synched_blink_proportions/');
PATHOUT_PLOTS                         = fullfile(PATH_PROCESSED,'14_observed_proportions_plots');

if ~exist(PATHOUT_PLOTS, 'dir')
    mkdir(PATHOUT_PLOTS);
end

% Run the main function:
ps_14_plot_observed_proportions( ...
    COND_CONTRASTS_BLINK_PAUSES, COND_CONTRASTS_SYNCHED_BLINKS, ...
    HALF_A4_SIZE,...
    LISTENER_BLINK_PAUSE_PROPORTIONS_PATH, SYNCHED_BLINK_PROPORTIONS_PATH, ...
    PATHOUT_PLOTS);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});

%% 3-b Plot individual histograms of the proportions after shuffling: -----

% Info for the function:
ALPHA               = 0.95;
FOUR_FIFTHS_A4_SIZE = [21-3.18, (29.7-2.54)*(4/5)]; % In cm

% Relevant paths:
PATHIN_LISTENER                 = fullfile(PATH_PROCESSED,'07_listener_blink_pause_proportions/');
PATHIN_SYNCHED_BLINKS           = fullfile(PATH_PROCESSED,'09_synched_blink_proportions/');
PATHIN_SPEAKER                  = fullfile(PATH_PROCESSED,'11_speaker_blink_pause_proportions/');
PATHIN_LISTENER_SHUFFLED        = fullfile(PATH_PROCESSED,'08_shuffled_listener_blink_pause_proportions/');
PATHIN_SYNCHED_BLINKS_SHUFFLED  = fullfile(PATH_PROCESSED,'10_shuffled_synched_blink_proportions/');
PATHIN_SPEAKER_SHUFFLED         = fullfile(PATH_PROCESSED,'12_shuffled_speaker_blink_pause_proportions/');
PATHOUT_PLOTS                   = fullfile(PATH_PROCESSED,'15_individual_shuffled_proportions_plots');

if ~exist(PATHOUT_PLOTS, 'dir')
    mkdir(PATHOUT_PLOTS);
end

% Run the main function:
ps_15_plot_shuffled_proportions_4_individuals(...
    ALPHA, FOUR_FIFTHS_A4_SIZE, ...
    PATHIN_LISTENER, PATHIN_SYNCHED_BLINKS, PATHIN_SPEAKER, ...
    PATHIN_LISTENER_SHUFFLED, PATHIN_SYNCHED_BLINKS_SHUFFLED, PATHIN_SPEAKER_SHUFFLED,...
    PATHOUT_PLOTS)

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});

%% 3-c Plot pie-charts to visualize the results of the binomial tests: ----

% Info for function:
HALF_A4_SIZE                  = [21-3.18, (29.7-2.54)/2]; % In cm

% Relevant paths:
PATHIN_BINOMIAL_RESULTS  = fullfile(PATH_PROCESSED,'13_binomial_test_results/');
PATHOUT_PLOT             = fullfile(PATH_PROCESSED,'16_binomial_pie_charts');

if ~exist(PATHOUT_PLOT, 'dir')
    mkdir(PATHOUT_PLOT);
end

% Run the main function:
ps_16_plot_binomial_pie_charts( ...
    HALF_A4_SIZE, PATHIN_BINOMIAL_RESULTS, PATHOUT_PLOT);

% Clear all variables defined in this section:
clearvars('-except', keepVars{:});

