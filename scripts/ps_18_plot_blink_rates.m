function ps_18_plot_blink_rates( ...
    COND_CONTRASTS, FOUR_FIFTHS_A4_SIZE,...
    PATHIN_RATES, PATHOUT_PLOTS)

load(fullfile(PATHIN_RATES,'listener_blink_rates.mat'), 'listenerBlinkRatesTbl');

% Get the condition names:
conditions = unique(listenerBlinkRatesTbl.Condition);
    
% Prepare a figure for the plots:
fig    = figure('Units','centimeters','InnerPosition',[1,1,FOUR_FIFTHS_A4_SIZE]);
hold on;
PropAxes = [];
nComparisons = size(COND_CONTRASTS, 1);

for c = 1:nComparisons

    % Segregate conditions:
    cond1      = COND_CONTRASTS{c, 1};
    cond2      = COND_CONTRASTS{c, 2};

    cond1Rows  = strcmp(cond1, listenerBlinkRatesTbl.Condition);
    cond2Rows  = strcmp(cond2, listenerBlinkRatesTbl.Condition);

    cond1Props = listenerBlinkRatesTbl.Proportion(cond1Rows);
    cond2Props = listenerBlinkRatesTbl.Proportion(cond2Rows);

    % Plot the proportions:
    subplot(1, nComparisons, c);
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

end


% Add a title:
sgtitle('Listener blink rates', 'FontName', 'Times New Roman');

% Set a color map to make individual participants stand out
colormap('lines');

% Make all of the proportion plots have the same y limits:
allPropYlims    = cat(1, PropAxes.YLim);
bestYlim        = [min(allPropYlims(:,1)), max(allPropYlims(:,2))];
[PropAxes.YLim] = deal(bestYlim);

% Add a y-axis label:
height          = PropAxes(1).Position(2) + PropAxes(1).Position(4) - ...
                  PropAxes(end).Position(2);
width           = PropAxes(end).Position(1) + PropAxes(end).Position(3) - ...
                  PropAxes(1).Position(1);
commonPropAxes  = axes('position', ...
                  [PropAxes(1).Position(1), PropAxes(1).Position(2), ...
                  width, height], ...
                  'Visible', 'off');
commonPropAxes.YLabel.Visible = 'on';
commonPropAxes.FontName       = 'Times New Roman';
commonPropAxes.FontSize       = 11;
ylabel(commonPropAxes,'Blinks per minute');


print(fig, fullfile(PATHOUT_PLOTS,'listener_blink_rates'), '-dpdf', '-fillpage');
print(fig, fullfile(PATHOUT_PLOTS,'listener_blink_rates'), '-dtiffn');
close;


end