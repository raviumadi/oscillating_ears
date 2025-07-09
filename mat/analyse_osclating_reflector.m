clear; close all;

%% === Load Data ===
[signal, fs] = audioread('recorded_response_15.wav');
audio = signal(:,1);          % Main signal
trigger = signal(:,2);        % Trigger signal

%% === Parameters ===
win_dur = 0.1;                         % 100 ms
win_samples = round(win_dur * fs);
trigger_thresh = 0.5 * max(trigger);
freq_range = [75e3 85e3];
nfft = 4096;
overlap = 1800;

%% === Detect Trigger Onsets ===
trigger_bin = trigger > trigger_thresh;
trigger_edges = find(diff(trigger_bin) == 1);

fprintf('Found %d segments\n', length(trigger_edges));

% Prepare Δf storage
deltaFs = zeros(length(trigger_edges),1);

% Layout grid for spectrograms
numPlots = length(trigger_edges);
cols = ceil(sqrt(numPlots));
rows = ceil(numPlots / cols);

figure('Name','All Spectrograms','NumberTitle','off','Position',[100 100 1400 900]);

for i = 1:length(trigger_edges)
    idx = trigger_edges(i);
    if idx + win_samples - 1 > length(audio)
        continue;
    end

    % === Extract Segment ===
    segment = audio(idx : idx + win_samples - 1);

        % === Compute Spectrogram ===
    [S,F,T] = spectrogram(segment, 2048, overlap, nfft, fs);
    S_mag = abs(S);
    S_dB = 20*log10(S_mag + eps);

    % Instantaneous freq = peak frequency in each column
    [~, peak_idx] = max(S_mag, [], 1);
    inst_freq = F(peak_idx);  % in Hz

    % === -20 dB Band Δf Estimation ===
    f_idx = F >= freq_range(1) & F <= freq_range(2);
    f_sub = F(f_idx);
    spec_slice = mean(S_dB(f_idx,:), 2);
    peak_dB = max(spec_slice);
    mask = spec_slice >= (peak_dB - 20);
    f_passband = f_sub(mask);

    if isempty(f_passband)
        deltaFs(i) = 0;
        f_min = NaN;
        f_max = NaN;
    else
        f_min = min(f_passband);
        f_max = max(f_passband);
        deltaFs(i) = f_max - f_min;
    end

    % === Plot ===
    subplot(rows, cols, i);
    imagesc(T, F/1000, S_dB); axis xy;
    ylim(freq_range / 1000);
    colormap(turbo);
    title(sprintf('Seg %d | Δf = %.0f Hz', i, deltaFs(i)));
    xlabel('Time (s)');
    ylabel('Freq (kHz)');
    set(gca,'FontSize',8);

    % === Overlay lines marking f_min and f_max ===
    hold on;
    if ~isnan(f_min)
        plot([T(1), T(end)], [f_min, f_min]/1000, 'w--', 'LineWidth', 1.2);
        plot([T(1), T(end)], [f_max, f_max]/1000, 'w--', 'LineWidth', 1.2);
    end
    hold off;
end

%% === Δf Summary Plot ===
figure('Name','Δf Summary','NumberTitle','off');
plot(deltaFs, 'o-', 'LineWidth', 1.5);
xlabel('Segment #');
ylabel('Estimated Δf (Hz)');
title('Doppler Spread (Δf) Across Segments (−3 dB Band)');
grid on;