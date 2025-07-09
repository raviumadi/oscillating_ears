files = {
    'recorded_response_0.wav'
    'recorded_response.wav'
    'recorded_response_5.wav'
    'recorded_response_10.wav'
    'recorded_response_15.wav'
    'recorded_response_20.wav'
    'recorded_response_25.wav'
    'recorded_response_30.wav'
    };

trials = ["Control", "Manual Fanning", "5Hz", "10Hz", "15Hz", "20Hz", "25Hz", "30Hz"];

% Parameters
winlen = 2048;
overlap = 1800;
nfft = 4096;
fs = 192e3;
win_samples = round(0.1 * fs);
freq_range = [75e3 85e3];

figure('Position', [100 100 1200 1000]);

for k = 1:length(files)
    [y, fs] = audioread(files{k});
    audio = y(:,1);
    trigger = y(:,2);

    trigger_thresh = 0.5 * max(trigger);
    trigger_edges = find(diff(trigger > trigger_thresh) == 1);
    N = length(trigger_edges);

    deltaFs = [];

    % === Compute Δf for each segment
    for i = 1:N-3
        idx = trigger_edges(i);
        if idx + win_samples - 1 > length(audio), continue; end

        segment = audio(idx : idx + win_samples - 1);
        [S,F,T] = spectrogram(segment, winlen, overlap, nfft, fs);
        S_dB = 20*log10(abs(S) + eps);

        f_idx = F >= freq_range(1) & F <= freq_range(2);
        f_sub = F(f_idx);
        spec_slice = mean(S_dB(f_idx,:), 2);
        peak_dB = max(spec_slice);
        mask = spec_slice >= (peak_dB - 20);
        f_band = f_sub(mask);

        if ~isempty(f_band)
            deltaFs(end+1) = max(f_band) - min(f_band);
        end
    end

    % === Use first segment for plotting
    idx = trigger_edges(2);
    segment = audio(idx : idx + win_samples - 1);
    [S1,F1,T1] = spectrogram(segment, winlen, overlap, nfft, fs);
    S1_dB = 20*log10(abs(S1) + eps);
    S1_dB = max(min(S1_dB, -20), -80);  % Clip for white background

    % Compute Δf for visual overlay
    f_idx = F1 >= freq_range(1) & F1 <= freq_range(2);
    f_sub = F1(f_idx);
    spec_slice = mean(S1_dB(f_idx,:), 2);
    peak_dB = max(spec_slice);
    mask = spec_slice >= (peak_dB - 20);
    f_band = f_sub(mask);
    if ~isempty(f_band)
        f_min = min(f_band);
        f_max = max(f_band);
    else
        f_min = NaN;
        f_max = NaN;
    end

    % === Plot
    subplot(4, 2, k);
    imagesc(T1, F1/1000, S1_dB); axis xy;
    ylim(freq_range / 1000);
    caxis([-80 -20]);
    cmap = turbo(256); cmap(1,:) = [1 1 1]; colormap(cmap);

    % Use average Δf over all segments for title
    mean_df = mean(deltaFs, 'omitnan');
    title(sprintf('%s --- mean $\\Delta f$ = %.0f\\ Hz', trials(k), mean_df), ...
        'Interpreter', 'latex');
    xlabel('Time (s)');
    ylabel('Freq (kHz)');
    set(gca, 'FontSize', 10);
    formatLatex(gca)
    % Overlay white dashed lines
    hold on;
    if ~isnan(f_min) && ~isnan(f_max)
        plot([T1(1), T1(end)], [f_min, f_min]/1000, 'w--', 'LineWidth', 1.2);
        plot([T1(1), T1(end)], [f_max, f_max]/1000, 'w--', 'LineWidth', 1.2);
    end
    hold off;
end

sgtitle('Frequency Spread from 80~kHz Tone Reflected off Oscillating Reflector', 'Interpreter', 'latex', 'FontSize', 18);
% exportgraphics(gcf, '/Users/ravi/Documents/projects/doppler/binaural_doppler/manuscript/fig/rec_spec_mean.pdf', 'Resolution', 300, 'Append', false)
%% === Figure: Plot Waveforms of Selected Segments ===
figure('Name', 'Waveforms of Selected Segments', 'Position', [100 100 1200 1000]);

for k = 1:length(files)
    [y, fs] = audioread(files{k});
    audio = y(:,1);
    trigger = y(:,2);

    trigger_thresh = 0.5 * max(trigger);
    trigger_edges = find(diff(trigger > trigger_thresh) == 1);
    if length(trigger_edges) < 2, continue; end

    idx = trigger_edges(2);  % segment used in previous spectrograms
    segment = audio(idx : idx + win_samples - 1);
    t = (0:win_samples-1) / fs * 1000;  % time in ms

    subplot(4, 2, k);
    plot(t, segment, 'k');
    ylim([-0.1 0.1])
    title(sprintf('%s', trials(k)), 'Interpreter', 'none');
    xlabel('Time (ms)');
    ylabel('Amplitude');
    set(gca, 'FontSize', 10);
    formatLatex(gca)
end
sgtitle('Waveforms of 80 kHz Tone after Bouncing off an Oscillating Reflector', 'FontSize', 18, 'Interpreter', 'latex');
% exportgraphics(gcf, '/Users/ravi/Documents/projects/doppler/binaural_doppler/manuscript/fig/rec_waveform.pdf', 'Resolution', 300, 'Append', false)
