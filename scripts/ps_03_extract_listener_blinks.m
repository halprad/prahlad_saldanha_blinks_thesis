% ps_03_extract_listener_blinks( ...
%     PATHIN_ICA, PATHOUT_LISTENER_BLINKS)
%
% This function extracts listener eye-blinks from the IC activations of the
% merged-by-session EEG datasets. The function saves the blink timestamps
% along with other helful functions including ones to help with epoching
% later in the pipeline.
%
% Inputs:
% PATHIN_ICA               - The folder in which the merged EEG datasets
%                            with the ICA weights attached are stored.
% PATHOUT_LISTENER_BLINKS  - The folder in which the blinks should be
%                            stored.


function ps_03_extract_listener_blinks( ...
    PATHIN_ICA, PATHOUT_LISTENER_BLINKS)


%% Get the names of the EEG datasets with the ICA weights attached: -------

IcaFiles = dir(fullfile(PATHIN_ICA,'*.set'));


%% For every dataset, extract listener eye-blinks by using the ICs: -------

ErrorLog = [];
badSets  = {};
wb = waitbar(0, 'Finding blinks...');

eeglab nogui;
for f = 1:length(IcaFiles)
    try

        subSessId = extractBefore(IcaFiles(f).name, '_ICA');


        %% Load a dataset and extract eye-blinks: -------------------------

        Eeg       = pop_loadset('filename', IcaFiles(f).name, 'filepath', PATHIN_ICA);
        Params    = checkBlinkerDefaults(struct(), getBlinkerDefaults(Eeg));
        Params.signalTypeIndicator = 'UseICs';
        Params.showMaxDistribution = false;
        [~, ~, Blinks, BlinkFits, BlinkProperties, BlinkStatistics, Params] = ...
                                           pop_blinker(Eeg, Params);


        %% Extract useful variables and save them with the blinks: --------

        eegTimes  = Eeg.times;
        EegEvents = Eeg.event;
        
        save(fullfile(PATHOUT_LISTENER_BLINKS, [subSessId,'_blinks.mat']), ...
            'Blinks', 'BlinkFits', 'BlinkProperties', ...
            'BlinkStatistics', 'Params', ...
            'eegTimes', 'EegEvents');

        waitbar(f/length(IcaFiles), wb);

    catch Me

        ErrorLog = [ErrorLog, Me];
        badSets  = [badSets; {subSessId}];

    end
end


%% Save the error messages and celebrate:

delete(wb);
save(fullfile(PATHOUT_LISTENER_BLINKS,'error_log'), "ErrorLog", "badSets");
disp('The listener eye-blinks were successfuly extracted.')
disp(['There were ', num2str(length(badSets)), ' bad datasets']);
load('splat.mat');
sound(y,Fs);

clear global; % To get rid of EEGLAB variables

end
