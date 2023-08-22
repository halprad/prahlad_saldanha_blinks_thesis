function ps_16_plot_shuffled_proportions_4_individuals(...
    ALPHA, FOUR_FIFTHS_A4_SIZE, ...
    PATHIN_LISTENER, PATHIN_SYNCHED_BLINKS, PATHIN_SPEAKER, ...
    PATHIN_LISTENER_SHUFFLED, PATHIN_SYNCHED_BLINKS_SHUFFLED, PATHIN_SPEAKER_SHUFFLED,...
    PATHOUT_PLOTS)

for t = 1:3 % For 3 types of tests

    % Set the variables depending on the test:
    switch t
        case 1
            load(fullfile(PATHIN_LISTENER, ...
                'listener_blink_pause_proportions.mat'), ...
                'listenerBlinkPauseProportionsTbl');
            ogTbl           = listenerBlinkPauseProportionsTbl;
            pathin          = PATHIN_LISTENER_SHUFFLED;
            grandTitleText  = 'Speech pauses with listener blinks after shuffling';
            plotFileName    = 'individual_shuffled_listener_blink_pauses';
            hasParticipants = true;
        case 2
            load(fullfile(PATHIN_SYNCHED_BLINKS, ...
                'synched_blink_proportions.mat'), ...
                'synchedBlinkProportionsTbl');
            ogTbl           = synchedBlinkProportionsTbl;
            pathin          = PATHIN_SYNCHED_BLINKS_SHUFFLED;
            grandTitleText  = 'Intervals after speaker blinks with listener blinks after shuffling';
            plotFileName    = 'individual_shuffled_synched_blinks';
            hasParticipants = true;
        case 3
            load(fullfile(PATHIN_SPEAKER, ...
                'speaker_blink_pause_proportions.mat'), ...
                'speakerBlinkPauseProportionsTbl');
            ogTbl           = speakerBlinkPauseProportionsTbl;
            pathin          = PATHIN_SPEAKER_SHUFFLED;
            grandTitleText  = 'Speech pauses with speaker blinks after shuffling';
            plotFileName    = 'individual_shuffled_speaker_blink_pauses';
            hasParticipants = false;
    end

    % Get the list of files:
    files           = dir(fullfile(pathin, '*shuffled*.mat'));

    % Subset the data and set the plotting parameters:
    if hasParticipants
        avRows          = strcmp('AV', ogTbl.Condition);
        ogProportions   = ogTbl.Proportion(avRows);
        ids             = cellfun(@(subs) extract(subs,digitsPattern), ...
                          ogTbl.Subject(avRows));
        titleText       = "Participant: ";
        nPlotRows       = ceil(length(ids)/3);
        nPlotColums     = floor(length(ids)/5);
    else
        ogProportions   = ogTbl.Proportion;
        ids             = ogTbl.Speaker;
        ids             = cellfun(@(id) id(1:2), ...
                          ids, 'UniformOutput', false);
        titleText       = "Speaker: ";
        nPlotRows       = ceil(length(ids)/2);
        nPlotColums     = floor(length(ids)/3);
    end

    % Create a full page figure and a variable for storing the axis
    % handles:
    fig     = figure('Units','centimeters','InnerPosition',[1,1,FOUR_FIFTHS_A4_SIZE]);
    allAxes = [];

    % For all participants/speakers, make histograms:
    for f = 1:length(files)

        % Load the data:
        load(fullfile(pathin,files(f).name), 'proportionsAfterShuffling');
        
        % Prepare to plot
        subplot(nPlotRows, nPlotColums, f);
        sortedShuffledProps = sort(proportionsAfterShuffling);
        confidenceInterval  = sortedShuffledProps( ...
                              ceil(ALPHA * length(sortedShuffledProps)));

        % Plot a histogram with xlines for the observed proportion and the
        % one-sided confidence interval:
        histogram(proportionsAfterShuffling);
        xline(ogProportions(f), ...
              'Color', 'red', ...
              'LineWidth', 1); 
        xline(confidenceInterval, ...
              'LineStyle', '--');

        % Add a simple title:
        title(titleText + ids{f}, 'FontWeight', 'normal');

        % Set the font, store the axis:
        ax                  =  gca();
        ax.FontName         = 'Times New Roman';
        ax.FontSize         = 9;
        allAxes             = [allAxes, ax];

    end % End of the loop over subplots

    % Make all of the subplots have the same axis limits:
    allAxesXlims    = cat(1, allAxes.XLim);
    allAxesYlims    = cat(1, allAxes.YLim);
    bestXlim        = [min(allAxesXlims(:,1)), max(allAxesXlims(:,2))];
    bestYlim        = [min(allAxesYlims(:,1)), max(allAxesYlims(:,2))];
    [allAxes.XLim]  = deal(bestXlim);
    [allAxes.YLim]  = deal(bestYlim);

    % Add the title:
    sgtitle(grandTitleText, 'FontName', 'Times New Roman', 'FontSize', 13);

    % Set the axis labels:
    height          = allAxes(1).Position(2) + allAxes(1).Position(4) - ...
                      allAxes(end).Position(2);
    width           = allAxes(nPlotColums).Position(1) + allAxes(nPlotColums).Position(3) - ...
                      allAxes(1).Position(1);
    lowerLeftIdx    = nPlotColums * (nPlotRows-1) + 1;
    commonAllAxes   = axes('position', ...
                           [allAxes(lowerLeftIdx).Position(1), ...
                            allAxes(lowerLeftIdx).Position(2), ...
                            width, height], ...
                           'Visible', 'off');
    commonAllAxes.XLabel.Visible = 'on';
    commonAllAxes.YLabel.Visible = 'on';
    commonAllAxes.FontName       = 'Times New Roman';
    commonAllAxes.FontSize       = 12;
    xlabel(commonAllAxes, 'Proportion');
    ylabel(commonAllAxes,'Number of permutations');

    % Save the plot as a full page tiff and a full page pdf:
    print(fig, fullfile(PATHOUT_PLOTS,plotFileName), '-dpdf', '-fillpage');
    print(fig, fullfile(PATHOUT_PLOTS,plotFileName), '-dtiffn');
    close;

end % End of the loop over types of tests