function formatLatexYY(ax)
% Applies LaTeX formatting to both axes and labels (including yyaxis)
% Usage: formatLatex(gca)

    set(ax, 'TickLabelInterpreter', 'latex', ...
            'FontSize', 16, ...
            'FontWeight', 'bold');

    grid(ax, 'on');
    grid(ax, 'minor');
    ax.GridAlpha = 0.5;
    ax.GridLineWidth = 0.5;

    % Format X label
    if ~isempty(get(get(ax, 'XLabel'), 'String'))
        set(get(ax, 'XLabel'), 'Interpreter', 'latex', ...
                               'FontSize', 16, ...
                               'FontWeight', 'bold');
    end

    % Format both YLabels if using yyaxis
    try
        yyaxis(ax, 'left');
        ylabL = get(ax, 'YLabel');
        if ~isempty(get(ylabL, 'String'))
            set(ylabL, 'Interpreter', 'latex', ...
                       'FontSize', 16, ...
                       'FontWeight', 'bold');
        end

        yyaxis(ax, 'right');
        ylabR = get(ax, 'YLabel');
        if ~isempty(get(ylabR, 'String'))
            set(ylabR, 'Interpreter', 'latex', ...
                       'FontSize', 16, ...
                       'FontWeight', 'bold');
        end
    catch
        % Fallback: standard single y-axis
        ylab = get(ax, 'YLabel');
        if ~isempty(get(ylab, 'String'))
            set(ylab, 'Interpreter', 'latex', ...
                      'FontSize', 16, ...
                      'FontWeight', 'bold');
        end
    end

    % Title
    if ~isempty(get(get(ax, 'Title'), 'String'))
        set(get(ax, 'Title'), 'Interpreter', 'latex', ...
                              'FontSize', 18, ...
                              'FontWeight', 'bold');
    end

    % Subtitle (optional, R2020b+)
    if isprop(ax, 'Subtitle')
        if ~isempty(get(get(ax, 'Subtitle'), 'String'))
            set(get(ax, 'Subtitle'), 'Interpreter', 'latex', ...
                                     'FontSize', 14, ...
                                     'FontWeight', 'bold');
        end
    end
end