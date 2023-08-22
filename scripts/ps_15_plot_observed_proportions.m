function ps_15_plot_observed_proportions( ...
    COND_CONTRASTS_BLINK_PAUSES, COND_CONTRASTS_SYNCHED_BLINKS, ...
    FOUR_FIFTHS_A4_SIZE,...
    LISTENER_BLINK_PAUSE_PROPORTIONS_PATH, SYNCHED_BLINK_PROPORTIONS_PATH, ...
    PATHOUT_PLOTS)

for f = 1:2 % For two figures:

    % Make a figure that's half a page big:
    fig = figure('Units','centimeters','InnerPosition',[1,1,FOUR_FIFTHS_A4_SIZE]);

    % Depending on the figure, use the appropriate variables:
    switch f
        case 1
            load(fullfile(LISTENER_BLINK_PAUSE_PROPORTIONS_PATH, ...
                'listener_blink_pause_proportions.mat'), ...
                'listenerBlinkPauseProportionsTbl');
            proportionsTbl     = listenerBlinkPauseProportionsTbl;
            currCondContrasts  = COND_CONTRASTS_BLINK_PAUSES;
            titleText          = 'Speech pauses with listener blinks';
            plotFileName       = 'blink_pause_comparisons';
        case 2
            load(fullfile(SYNCHED_BLINK_PROPORTIONS_PATH, ...
                'synched_blink_proportions.mat'), ...
                'synchedBlinkProportionsTbl');
            proportionsTbl     = synchedBlinkProportionsTbl;
            currCondContrasts  = COND_CONTRASTS_SYNCHED_BLINKS;
            titleText          = 'Intervals after speaker blinks with listener blinks';
            plotFileName       = 'synched_blink_comparisions';
    end

    PropAxes  = [];
    LogitAxes = [];

    nComparisons = size(currCondContrasts, 1);

    for c = 1:nComparisons

        % Segregate conditions:
        cond1      = currCondContrasts{c, 1};
        cond2      = currCondContrasts{c, 2};

        cond1Rows  = strcmp(cond1, proportionsTbl.Condition);
        cond2Rows  = strcmp(cond2, proportionsTbl.Condition);

        cond1Props = proportionsTbl.Proportion(cond1Rows);
        cond2Props = proportionsTbl.Proportion(cond2Rows);

        % Plot the proportions:
        subplot(nComparisons, 2, c*2-1);
        scatter(zeros(15,1), cond1Props, ...
                50, 1:length(cond1Props), ...
                'filled', 'o', 'MarkerFaceAlpha', 0.5);
        hold on;
        scatter(ones(15,1),  cond2Props, ...
                50, 1:length(cond1Props), ...
               'filled', 'o', 'MarkerFaceAlpha', 0.5);
        plot([0,1], [cond1Props,cond2Props],  'Color', [0.8,0.8,0.8], 'LineStyle',':');
        ax            = gca();
        ax.XTick      = [0,1];
        ax.XLim       = [-0.5,1.5];
        ax.XTickLabel = {cond1, cond2};
        ax.TickDir    = 'none';
        ax.FontName   = 'Times New Roman';
        ax.FontSize   = 11;
        PropAxes      = [PropAxes, ax];

        % Plot the logit of the proportions:
        subplot(nComparisons, 2, c*2);
        scatter(zeros(15,1), local_convert_2_logits(cond1Props), ...
                50, 1:length(cond1Props), ...
                'filled', 'o', 'MarkerFaceAlpha', 0.5);
        hold on;
        scatter(ones(15,1),  local_convert_2_logits(cond2Props), ...
                50, 1:length(cond1Props), ...
                'filled', 'o', 'MarkerFaceAlpha', 0.5);
        plot([0,1], local_convert_2_logits([cond1Props,cond2Props]),  'Color', [0.8,0.8,0.8], 'LineStyle',':');
        ax            = gca();
        ax.XTick      = [0,1];
        ax.XLim       = [-0.5,1.5];
        ax.XTickLabel = {cond1, cond2};
        ax.TickDir    = 'none';
        ax.FontName   = 'Times New Roman';
        ax.FontSize   = 11;
        LogitAxes     = [LogitAxes, ax];

    end % End of loop over condition comparisions

    % Add a title:
    sgtitle(titleText, 'FontName', 'Times New Roman');

    % Set a color map to make individual participants stand out
    colormap('lines');

    % Make all of the proportion plots have the same y limits:
    allPropYlims    = cat(1, PropAxes.YLim);
    bestYlim        = [min(allPropYlims(:,1)), max(allPropYlims(:,2))];
    [PropAxes.YLim] = deal(bestYlim);

    % Make all of the logit plots have the same y-limits:
    allLogitYlims    = cat(1, LogitAxes.YLim);
    bestYlim         = [min(allLogitYlims(:,1)), max(allLogitYlims(:,2))];
    [LogitAxes.YLim] = deal(bestYlim);

    % Add a y-axis label to the left column:
    height          = PropAxes(1).Position(2) + PropAxes(1).Position(4) - ...
                      PropAxes(end).Position(2);
    width           = PropAxes(1).Position(3);
    commonPropAxes  = axes('position', ...
                           [PropAxes(end).Position(1), PropAxes(end).Position(2), ...
                            width, height], ...
                           'Visible', 'off');
    commonPropAxes.YLabel.Visible = 'on';
    commonPropAxes.FontName       = 'Times New Roman';
    commonPropAxes.FontSize       = 11;
    ylabel(commonPropAxes,'Proportion');

    % Add a y-axis label to the right column:
    height           = LogitAxes(1).Position(2) + LogitAxes(1).Position(4) - ...
                       LogitAxes(end).Position(2);
    width            = LogitAxes(1).Position(3);
    commonLogitAxes  = axes('position', ...
                           [LogitAxes(end).Position(1), LogitAxes(end).Position(2), ...
                            width, height], ...
                           'Visible', 'off');
    commonLogitAxes.YLabel.Visible = 'on';
    commonLogitAxes.FontName       = 'Times New Roman';
    commonLogitAxes.FontSize       = 11;
    ylabel(commonLogitAxes,'Logit of proportion');

    % Save a full page pdf and a 4/5 page tiff:
    print(fig, fullfile(PATHOUT_PLOTS,plotFileName), '-dpdf', '-fillpage');
    print(fig, fullfile(PATHOUT_PLOTS,plotFileName), '-dtiffn');
    close;

end % End of loop over figures

end % End of the main function



% -------------------------------------------------------------------------
%% ****************** Local function **************************************
% -------------------------------------------------------------------------

%% local_convert_2_logits()
% This functions converts proportions to logits of the proportions.

function propsInLogits = local_convert_2_logits(proportions)
propsInOdds   = proportions ./ (1-proportions);
propsInLogits = log(propsInOdds);
end

