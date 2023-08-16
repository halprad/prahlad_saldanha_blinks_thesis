function ps_16_plot_binomial_pie_charts( ...
    HALF_A4_SIZE, PATHIN_BINOMIAL_RESULTS, PATHOUT_PLOT)


load(fullfile(PATHIN_BINOMIAL_RESULTS,'binomial_test_results.mat'), 'binomResultsTbl');

subplotTitles  = ["Listeners who blinked during pauses", ...
                  "Listeners who blinked after speakers did"...
                  "Speakers who blinked during pauses"];
subplotPos     = {[1,2], [7,8,9], [4,5]};

% Make a figure that's half a page big:
fig = figure('Units','centimeters','InnerPosition',[4,4,HALF_A4_SIZE]);

for t = 1:length(binomResultsTbl.Test)
    counts              = [binomResultsTbl.NSignificant(t), ...
                           binomResultsTbl.NTotal(t) - binomResultsTbl.NSignificant(t)];
    standOut            = [true, false];
    numSignificant      = counts(1);
    percentNot          = counts(2) - numSignificant;
    labels              = {num2str(counts(1)), num2str(counts(2))};

    subplot(2, 5, subplotPos{t});
    pObj = pie(counts, standOut, labels);

    % Adjust the text around the pie chart:
    [pObj([2,4]).FontName] = deal('Times New Roman');
    [pObj([2,4]).FontSize] = deal(9);

    title(gca, subplotTitles(t), 'FontName', 'Times New Roman','FontSize', 12);

end

colormap('summer');
legend(["Significant", "Not significant"], 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'Position', [0.68,0.24,0.19,0.07]);


% Save a full page pdf and a half page tiff:
print(fig, fullfile(PATHOUT_PLOT,'binomial_pie_charts'), '-dpdf', '-fillpage');
print(fig, fullfile(PATHOUT_PLOT,'binomial_pie_charts'), '-dtiffn');
close;

end % End of main function