%function ps_02_do_ica( ...
%    LOW_CUT_OFF, HIGH_CUT_OFF, PRUNE, ...
%    PATH_TOOLBOX, PATHIN_RAW, PATHOUT_PROCESSED)
%
% This function merges EEG datasets for every session for every participant
% and calculates the ICA weights. It saves one dataset per session per
% participant with the ICA weights attached.
%
% Inputs:
% LOW_CUT_OFF       - High-pass filter order cut-off in Hz (for filtering
%                     prior to the ICA decomposition).
% HIGH_CUT_OFF      - Low-pass filter order cut-off in Hz (for filtering
%                     prior to the ICA decomposition).
% PRUNE             - Threshold in standard deviations beyong which
%                     artifactual epochs should be rejected based on their
%                     joint-probabilities.
% PATH_TOOLBOX      - The path to where the toolboxes are kept.
% PATHIN_RAW        - The path to where the raw EEG data is stored.
% PATHOUT_ICA       - The path to where the merged datasets with their ICA
%                     weights attached should be stored.


function ps_02_do_ica( ...
    LOW_CUT_OFF, HIGH_CUT_OFF, PRUNE, ...
    PATH_TOOLBOX, PATHIN_RAW, PATHOUT_PROCESSED)

%% Get information about the subjects: ------------------------------------

SubjFolders  = dir(fullfile(PATHIN_RAW,'Sub*'));
subjects     = {SubjFolders.name};
sessions     = {'A', 'B'};


%% Compute the ICA weight per session per subject: ------------------------

% In each session there were multiple EEG data files for each subject; one
% per 30 second experimental condition. For every subject and for every
% session, one long dataset is produced by merging the the datasets of all
% the conditions. The ICA weights of this dataset are then computed, and
% it is stored with the weights attached.

ErrorLog = [];
badSets  = {};
wb       = waitbar(0, 'Computing ICs...');

eeglab nogui;
for s = 1:length(subjects)
    for ss = 1:length(sessions)
        try
            AllEeg = [];

            %% Get the paths to each of the datasets to be merged: --------

            subSessPath  = fullfile(PATHIN_RAW, subjects{s}, ...
                                    ['Session_', sessions{ss}]);
            SubSessFiles = dir(fullfile(subSessPath,'*.set'));


            %% Load and merge the datasets: -------------------------------

            Eeg           = pop_loadset('filename', {SubSessFiles.name}, ...
                                        'filepath', subSessPath);
            [AllEeg, Eeg] = eeg_store(AllEeg, Eeg);
            Eeg           = pop_mergeset(AllEeg, 1:length(AllEeg));


            %% Prepare the merged dataset: --------------------------------

            Eeg    = pop_chanedit   (Eeg,    'lookup', fullfile(PATH_TOOLBOX,'eeglab2023.0\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc'));
            EegIca = pop_eegfiltnew (Eeg,    'locutoff', LOW_CUT_OFF);
            EegIca = pop_eegfiltnew (EegIca, 'hicutoff', HIGH_CUT_OFF);
            EegIca = eeg_regepochs  (EegIca);
            EegIca = pop_jointprob  (EegIca, 1, 1:24, PRUNE, PRUNE, 1, 1);


            %% Compute the ICA weights: -----------------------------------

            EegIca         = pop_runica(EegIca, 'icatype', 'runica', 'extended', 1);
            Eeg.icawinv    = EegIca.icawinv;
            Eeg.icasphere  = EegIca.icasphere;
            Eeg.icaweights = EegIca.icaweights;
            Eeg            = eeg_checkset(Eeg);


            %% Save the merged dataset with its ICA weights: --------------

            Eeg.setname    = [subjects{s}, '_', sessions{ss},'_ICA'];
            Eeg            = pop_saveset(Eeg, [Eeg.setname, '.set'], PATHOUT_PROCESSED);

            progress = ((s-1)*length(sessions) + ss) / length(subjects)*length(sessions);
            waitbar(progress, wb);


        catch Me

            %% Store any error messages: ----------------------------------

            ErrorLog = [ErrorLog, Me];
            badSets  = [badSets, {subjects{s},sessions{ss}}];

        end
    end
end

%% Save the error messages and celebrate: ---------------------------------

save(fullfile(PATHOUT_PROCESSED, 'error_log.mat'), "ErrorLog", "badSets");
delete(wb); 
disp('The ICA decompositions completed successfully.');
disp(['There were ', num2str(length(badSets)), ' bad sets']);
load('handel.mat'); sound(y,Fs);

clear global; % To get rid of global EEGLAB variables.

end
