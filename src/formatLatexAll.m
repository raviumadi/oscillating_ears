function formatLatexAll(fig)
% Applies LaTeX formatting to all axes in the given figure (or gcf if not provided)

    if nargin < 1
        fig = gcf;
    end

    axAll = findall(fig, 'type', 'axes');
    for ax = axAll(:)'  % row vector for consistency
        formatLatexYY(ax);
    end
end