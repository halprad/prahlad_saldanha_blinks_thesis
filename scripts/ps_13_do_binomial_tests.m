function ps_13_do_binomial_tests(...
    CONDITION, THRESH_INDIVIDUAL_SUBS, THRESH_BINOMIAL_TESTS, ...
    PATHIN_LISTENER_PROPORTIONS, PATHIN_SHUFFLED_LISTENER_PROPORTIONS, ...
    PATHIN_SYNCHED_BLINKS_PROPROTIONS, PATHIN_SHUFFLED_SYNCHED_BLINKS_PROPORTIONS, ...
    PATHIN_SPEAKER_PROPORTIONS, PATHIN_SHUFFLED_SPEAKER_PROPORTIONS, ...
    PATHOUT_BINOMIAL_TEST_RESULTS)

%% Load the tables of the actual proportions, get the filenames of the proportions after shuffling:


% Proportions of listener blinks in pauses:

load(fullfile(PATHIN_LISTENER_PROPORTIONS,'listener_blink_pause_proportions.mat'), ...
    'listenerBlinkPauseProportionsTbl');

rows = strcmp(listenerBlinkPauseProportionsTbl.Condition, CONDITION);
listenerProportions = listenerBlinkPauseProportionsTbl.Proportion(rows, :);
listenerIds = listenerBlinkPauseProportionsTbl.Subject(rows, :);

ShuffledListenerPropFiles = dir(...
    fullfile(PATHIN_SHUFFLED_LISTENER_PROPORTIONS, '*shuffled*.mat') );


% Proportions of listener blinks in post-speaker-blink intervals:

load(fullfile(PATHIN_SYNCHED_BLINKS_PROPROTIONS,'synched_blink_proportions.mat'), ...
    'synchedBlinkProportionsTbl');

rows = strcmp(synchedBlinkProportionsTbl.Condition, CONDITION);
synchedBlinkProportions = synchedBlinkProportionsTbl.Proportion(rows, :);
synchedBlinkIds = synchedBlinkProportionsTbl.Subject(rows, :);

ShuffledSynchedBlinkPropFiles = dir(...
    fullfile(PATHIN_SHUFFLED_SYNCHED_BLINKS_PROPORTIONS, '*shuffled*.mat') );


% Proportions of speaker blinks in pauses:

load(fullfile(PATHIN_SPEAKER_PROPORTIONS,'speaker_blink_pause_proportions.mat'), ...
    'speakerBlinkPauseProportionsTbl');

speakerProportions = speakerBlinkPauseProportionsTbl.Proportion;
speakerIds = speakerBlinkPauseProportionsTbl.Speaker;

ShuffledSpeakerPropFiles = dir(...
    fullfile(PATHIN_SHUFFLED_SPEAKER_PROPORTIONS, '*shuffled*.mat') );


%% Count the number of listeners whose blinks in pause proportions were significant:

listenerStatsCells = {};
numListeners = length(listenerProportions);
numSigListeners = 0;

for l = 1:length(listenerProportions)

    % Load the proportions after shuffling:
    load(fullfile(PATHIN_SHUFFLED_LISTENER_PROPORTIONS,...
                  ShuffledListenerPropFiles(l).name), ...
         'proportionsAfterShuffling');

    % Count the number below the actual proportion and the total:
    numBelow = sum(listenerProportions(l) > proportionsAfterShuffling);
    numTotal = length(proportionsAfterShuffling);

    % Calculate the p-value
    pVal = 1 - numBelow/numTotal;

    % Increment the count of significant cases accordingly:
    if pVal <= THRESH_INDIVIDUAL_SUBS(1)
        numSigListeners = numSigListeners + 1;
    end

    % Store the results:
    listenerStatsCells = [listenerStatsCells; ...
                          {listenerIds{l}, pVal}];

end

% Calculate the p-value that the number of listeners with significant
% results was significant and store the information:
binomPVal = binocdf(numSigListeners-1, numListeners, ...
                    THRESH_BINOMIAL_TESTS(1), 'upper');

binomResultsCells(1,:) = {'ListenerBlinkPauses', numSigListeners, ...
                           numListeners, binomPVal};


%% Count the number of listeners whose blinks in post-blinks were significant:

synchedBlinksStatsCells = {};
numSynchedBlinks = length(synchedBlinkProportions);
numSigSynchedBlinks = 0;

for b = 1:length(synchedBlinkProportions)

    % Load the proportions after shuffling:
    load(fullfile(PATHIN_SHUFFLED_SYNCHED_BLINKS_PROPORTIONS,...
                  ShuffledSynchedBlinkPropFiles(b).name), ...
         'proportionsAfterShuffling');

    % Count the number below the actual proportion and the total:
    numBelow = sum(synchedBlinkProportions(b) > proportionsAfterShuffling);
    numTotal = length(proportionsAfterShuffling);

    % Calculate the p-value
    pVal = 1 - numBelow/numTotal;

    % Increment the count of significant cases accordingly:
    if pVal <= THRESH_INDIVIDUAL_SUBS(2)
        numSigSynchedBlinks = numSigSynchedBlinks + 1;
    end

    % Store the results:
    synchedBlinksStatsCells = [synchedBlinksStatsCells; ...
                               {synchedBlinkIds{b}, pVal}];

end

% Calculate the p-value that the number of listeners with significant
% results was significant and store the information:
binomPVal = binocdf(numSigSynchedBlinks-1, numSynchedBlinks, ...
                    THRESH_BINOMIAL_TESTS(2), 'upper');

binomResultsCells(2,:) = {'SynchedBlinks', numSigSynchedBlinks, ...
                           numSynchedBlinks, binomPVal};

%% Count the number of listeners whose blinks in post-blinks were significant:

speakerBlinksStatsCells = {};
numSpeakers = length(speakerProportions);
numSigSpeakers = 0;

for s = 1:length(speakerProportions)

    % Load the proportions after shuffling:
    load(fullfile(PATHIN_SHUFFLED_SPEAKER_PROPORTIONS,...
                  ShuffledSpeakerPropFiles(s).name), ...
         'proportionsAfterShuffling');

    % Count the number below the actual proportion and the total:
    numBelow = sum(speakerProportions(s) > proportionsAfterShuffling);
    numTotal = length(proportionsAfterShuffling);

    % Calculate the p-value
    pVal = 1 - numBelow/numTotal;

    % Increment the count of significant cases accordingly:
    if pVal <= THRESH_INDIVIDUAL_SUBS(3)
        numSigSpeakers = numSigSpeakers + 1;
    end

    % Store the results:
    speakerBlinksStatsCells = [speakerBlinksStatsCells; ...
                               {speakerIds{s}, pVal}];

end

% Calculate the p-value that the number of speakers with significant
% results was significant and store the information:
binomPVal = binocdf(numSigSpeakers-1, numSpeakers, ...
                    THRESH_BINOMIAL_TESTS(3), 'upper');

binomResultsCells(3,:) = {'SpeakerBlinkPauses', numSigSpeakers, ...
                           numSpeakers, binomPVal};

%% Save the results: ------------------------------------------------------

% Save the statistics of the listener blinks in pauses:

listenerStatsTbl = cell2table(listenerStatsCells, ...
                             'VariableNames', {'Listener','PValue'});

save(fullfile(PATHOUT_BINOMIAL_TEST_RESULTS, ...
              'listener_blink_pause_proportions_pvalues.mat'), ...
     'listenerStatsTbl');

writetable(listenerStatsTbl, ...
           fullfile(PATHOUT_BINOMIAL_TEST_RESULTS, ...
                    'listener_blink_pause_proportions_pvalues.xlsx'));


% Save the statistics of the listener blinks in post-speaker-blink
% intervals:

synchedBlinkStatsTbl = cell2table(synchedBlinksStatsCells, ...
                             'VariableNames', {'Listener','PValue'});

save(fullfile(PATHOUT_BINOMIAL_TEST_RESULTS, ...
              'synched_blink_proportions_pvalues.mat'), ...
     'synchedBlinkStatsTbl');

writetable(synchedBlinkStatsTbl, ...
           fullfile(PATHOUT_BINOMIAL_TEST_RESULTS, ...
                    'synched_blink_proportions_pvalues.xlsx'));


% Save the statistics of the speaker blinks in pauses:

speakerStatsTbl = cell2table(speakerBlinksStatsCells, ...
                             'VariableNames', {'Speaker','PValue'});

save(fullfile(PATHOUT_BINOMIAL_TEST_RESULTS, ...
              'speaker_blink_pause_proportions_pvalues.mat'), ...
     'speakerStatsTbl');

writetable(speakerStatsTbl,...
            fullfile(PATHOUT_BINOMIAL_TEST_RESULTS, ...
                    'speaker_blink_pause_proportions_pvalues.xlsx'));

% Save the results of the binomial tests:

binomResultsTbl = cell2table(binomResultsCells, ...
                             'VariableNames', {'Test','NSignificant', ...
                                               'NTotal','PValue'});

save(fullfile(PATHOUT_BINOMIAL_TEST_RESULTS, ...
              'binomial_test_results.mat'),...
     'binomResultsTbl');

writetable(binomResultsTbl,...
           fullfile(PATHOUT_BINOMIAL_TEST_RESULTS, ...
           'binomal_test_results.xlsx'));

disp('All binomial tests were completed and the results were saved successfully.');
