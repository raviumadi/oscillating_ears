T = readtable('binaural_summary.csv');

% Round columns
T.TipVel = round(T.TipVel, 2);
T.DeltaF_BW = round(T.DeltaF_BW, 1);
T.DeltaF_Std = round(T.DeltaF_Std, 1);
T.BLD_Std = round(T.BLD_Std, 1);

% Rename for clarity
T.Properties.VariableNames = {'Run','EarFreq_Hz','ThetaDeg','TipVel', 'TipDisp_mm', 'DeltaF_BW','DeltaF_Std','BLD_Std'};

fid = fopen('binaural_summary_table.tex', 'w');

fprintf(fid, '%% LaTeX Table: Binaural Doppler Summary\n');
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\scriptsize\n');
fprintf(fid, '\\renewcommand{\\arraystretch}{1.15}\n');
fprintf(fid, '\\begin{tabularx}{\\linewidth}{c c c c c c c c}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '\\textbf{Trial} & \\textbf{Ear Freq(Hz)} & \\textbf{Angle($\\degree$)} & \\textbf{Tip Vel.(m/s)} & \\textbf{Excur.(mm)} & \\textbf{$\\Delta f$ BW(Hz)} & \\textbf{$\\Delta f$ $\\sigma$(Hz)} & \\textbf{BLD $\\sigma$(dB)} \\\\\n');
fprintf(fid, '\\midrule\n');

for i = 1:height(T)
    fprintf(fid, '%d & %d & %d & %.2f & %.1f & %.1f & %.1f & %.1f \\\\\n', ...
        T.Run(i), T.EarFreq_Hz(i), T.ThetaDeg(i), ...
        T.TipVel(i), T.TipDisp_mm(i), T.DeltaF_BW(i), T.DeltaF_Std(i), T.BLD_Std(i));
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabularx}\n');
fprintf(fid, '\\caption{Summary of simulated binaural Doppler runs showing ear oscillation frequency, angular amplitude, computed tip velocity, frequency bandwidth ($\\Delta f$), its standard deviation, and binaural level difference (BLD) standard deviation.}\n');
fprintf(fid, '\\label{tab:binaural_summary}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);
disp('LaTeX table written to binaural_summary_table.tex');