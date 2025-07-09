clear; clc;

%% Parameters
duration = 100;                % Duration in milliseconds
fs = 192e3;                    % Sampling rate (Hz)
cf_freq = 80e3;                % CF bat call frequency (Hz)
bat_velocity = 0;              % Forward flight velocity (m/s)
sound_speed = 343;             % Speed of sound (m/s)

% Ear motion parameters
ear_height = 0.02;             % Ear height in metres (e.g. 3 cm)
theta_max_deg = 30;            % Max angular displacement (degrees)
theta_max = deg2rad(theta_max_deg/2);  % Convert to radians
ear_freq = 20;                 % Ear oscillation frequency (Hz)
omega_max = 2 * pi * ear_freq;         % Max angular velocity (rad/s)

% Compute max linear velocity at tip
v_tip = theta_max * omega_max * ear_height;
fprintf('Computed max tip velocity: %.3f m/s\n', v_tip);

%% 1. Generate CF Call
[cf_call] = generateCFBatCall(cf_freq, duration, fs, 0, bat_velocity);
t = (0:length(cf_call)-1) / fs;
cf_shifted = cf_call;


%% 2. Simulate Echoes with Time-Varying Angular Ear Motion (with phase warping)

n_segments = 100;
segment_positions = linspace(0, ear_height, n_segments);  % base to tip
segment_velocity = linspace(0, v_tip, n_segments);         % scaled velocity
delays = segment_positions / sound_speed;  % one-way delay only
dt = 1/fs;

echo_left = zeros(size(cf_call));
echo_right = zeros(size(cf_call));

% Use this switch to toggle averaging (true) or summation (false)
useAveraging = false;

% Optional: store how many contributions per sample if averaging
if useAveraging
    contrib_count_L = zeros(size(cf_call));
    contrib_count_R = zeros(size(cf_call));
end

for i = 1:n_segments
    delay_samples = round(delays(i) * fs);
    v_seg = segment_velocity(i);

    phiL = 0;  % cumulative phase for left
    phiR = 0;  % cumulative phase for right

    for j = 1:length(t)
        % Oscillating linear velocity at time j
        vL = v_seg * sin(2*pi*ear_freq*t(j));
        vR = -v_seg * sin(2*pi*ear_freq*t(j));

        % Effective receiver velocity including forward flight
        v_eff_L = bat_velocity + vL;
        v_eff_R = bat_velocity + vR;

        % Instantaneous Doppler-shifted frequency
        fL = cf_freq * (1 + v_eff_L / sound_speed);
        fR = cf_freq * (1 + v_eff_R / sound_speed);

        % Integrate phase
        phiL = phiL + 2*pi*fL*dt;
        phiR = phiR + 2*pi*fR*dt;

        idxL = j + delay_samples;
        idxR = j + delay_samples;

        if idxL <= length(cf_call)
            echo_left(idxL) = echo_left(idxL) + sin(phiL);
            if useAveraging
                contrib_count_L(idxL) = contrib_count_L(idxL) + 1;
            end
        end
        if idxR <= length(cf_call)
            echo_right(idxR) = echo_right(idxR) + sin(phiR);
            if useAveraging
                contrib_count_R(idxR) = contrib_count_R(idxR) + 1;
            end
        end
    end
end

% Apply averaging if enabled
if useAveraging
    echo_left(contrib_count_L > 0) = echo_left(contrib_count_L > 0) ./ contrib_count_L(contrib_count_L > 0);
    echo_right(contrib_count_R > 0) = echo_right(contrib_count_R > 0) ./ contrib_count_R(contrib_count_R > 0);
end

% Normalize echoes
echo_left = echo_left / max(abs(echo_left));
echo_right = echo_right / max(abs(echo_right));

%% 3. Time-Varying Cross-Correlation

window_size = round(0.005 * fs);
step_size = round(0.001 * fs);
n_steps = floor((length(cf_call) - window_size) / step_size);

corr_L = zeros(n_steps, 1);
corr_R = zeros(n_steps, 1);

for i = 1:n_steps
    idx = (1:window_size) + (i-1)*step_size;
    x = cf_shifted(idx);
    yL = echo_left(idx);
    yR = echo_right(idx);

    [cL, ~] = xcorr(x, yL, 'coeff');
    [cR, ~] = xcorr(x, yR, 'coeff');
    corr_L(i) = max(cL);
    corr_R(i) = max(cR);
end

time_axis = ((1:n_steps) * step_size) / fs * 1e3;

%% 4. Plot Time-Domain, Cross-Correlation, and Instantaneous Frequency

v_tip_L = v_tip * sin(2*pi*ear_freq*t);
v_tip_R = -v_tip * sin(2*pi*ear_freq*t);
v_norm_L = v_tip_L / max(abs(v_tip_L));
v_norm_R = v_tip_R / max(abs(v_tip_R));

% Time vector
time_ms = t * 1e3;

% Ensure row vectors
time_ms = time_ms(:)';
echo_left = echo_left(:)';
echo_right = echo_right(:)';

figure;

subplot(4,1,1);
plot(t*1e3, cf_call, 'color', [0.5 0.5 0.5]);
title('CF Call Waveform -- 80~kHz');
xlabel('Time (ms)'); ylabel('Amplitude');
formatLatex(gca)

subplot(4,1,2);
% yyaxis right;
plot(time_ms, v_norm_L, 'Color', 'b', 'LineWidth', 1.2); hold on
plot(time_ms, v_norm_R, 'Color', [1, 0.5, 0], 'LineStyle', '-', 'LineWidth', 1.2);
ylabel('Norm. Velocity');

% yyaxis left;
% Plot left ear with 50% opacity (black)
% Patch Left (Blue)
        patch([time_ms, fliplr(time_ms)], ...
            [echo_left, zeros(size(echo_left))], ...
            'blue', 'FaceAlpha', 0.9, 'EdgeColor', 'none'); hold on;

        % Patch Right (Orange)
        patch([time_ms, fliplr(time_ms)], ...
            [echo_right, zeros(size(echo_right))], ...
            [1 0.5 0], 'FaceAlpha', 0.9, 'EdgeColor', 'none');
xlabel('Time (ms)'); ylabel('Amplitude');
title(['Synth. Echo and Velocity -- ' num2str(ear_freq) '~Hz, '...
    num2str(theta_max_deg) '$^o$ \& ' '$v_{tip}$ = ' num2str(round(v_tip,2)) '~m/s']);

formatLatex(gca)

subplot(4,1,3);
plot(time_axis, corr_L, 'b', 'LineWidth', 1.5); hold on;
plot(time_axis, corr_R, 'color', [1, 0.5, 0,], 'LineWidth', 1.5);
ylim([0 1])
title('Time-Varying Cross-Correlation');
xlabel('Time (ms)'); ylabel('Max Correlation');
legend('Left', 'Right', 'Location', 'northwest',  'Orientation', 'Horizontal', 'interpreter', 'latex');
formatLatex(gca)

% 5. Clean Instantaneous Frequency Estimate (Hilbert + Envelope Masking)

% Hilbert analytic signals
analytic_left = hilbert(echo_left');
analytic_right = hilbert(echo_right');

% Instantaneous phase → frequency (Hz)
phase_left = unwrap(angle(analytic_left));
phase_right = unwrap(angle(analytic_right));
instfreq_left_raw = [0; diff(phase_left)] * fs / (2 * pi);
instfreq_right_raw = [0; diff(phase_right)] * fs / (2 * pi);

% Envelopes
env_left = abs(analytic_left);
env_right = abs(analytic_right);

% Envelope threshold
threshold = 0.01 * max([env_left; env_right]);

% Apply masking based on envelope
instfreq_left_masked = instfreq_left_raw;
instfreq_right_masked = instfreq_right_raw;
instfreq_left_masked(env_left < threshold) = NaN;
instfreq_right_masked(env_right < threshold) = NaN;

% Smooth masked signals
instfreq_left_clean = movmean(instfreq_left_masked, 50, 'omitnan');
instfreq_right_clean = movmean(instfreq_right_masked, 50, 'omitnan');


% Plot
subplot(4,1,4);
plot(time_ms, instfreq_left_clean, 'b', 'LineWidth', 1.5); hold on;
plot(time_ms, instfreq_right_clean, 'color', [1, 0.5, 0,], 'LineWidth', 1.5);
yline(cf_freq, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
title('Instantaneous Frequency of Echo vs. CF Frequency');
xlabel('Time (ms)');
ylabel('Frequency (kHz)');
legend('Left Ear', 'Right Ear', 'Location', 'northwest', 'Orientation', 'Horizontal', 'interpreter', 'latex');
ylim([78000 82000])
yticklabels(yticks./1000)
formatLatex(gca)
sgtitle('Binaural Simulation: Oscillating Ear Motion and Doppler Effects', 'FontSize', 18, 'interpreter', 'latex');

%% export
% exportgraphics(gcf, '/Users/ravi/Documents/projects/doppler/binaural_doppler/manuscript/fig/sim_method.pdf', 'Resolution', 300, 'Append', false)
%% 6. Spectrograms of Both Ears

figure;
subplot(2,1,1);
spectrogram(echo_left, 1024, 800, 1024, fs, 'yaxis');
title('Spectrogram of Left Ear Echo');

subplot(2,1,2);
spectrogram(echo_right, 1024, 800, 1024, fs, 'yaxis');
title('Spectrogram of Right Ear Echo');

%% 8. Ear Tip Velocity, Δf and Instantaneous Binaural Level Difference (BLD)

% Compute tip velocity profiles (sinusoidal over time)
v_tip_L =  v_tip * sin(2*pi*ear_freq*t);   % left tip velocity
v_tip_R = -v_tip * sin(2*pi*ear_freq*t);   % right tip velocity

% Instantaneous frequency difference between ears
delta_f = instfreq_right_clean - instfreq_left_clean;

% Binaural level difference (BLD in dB)
env_left = abs(hilbert(echo_left));
env_right = abs(hilbert(echo_right));
BLD_dB = 20 * log10(env_right ./ env_left);
BLD_dB(~isfinite(BLD_dB)) = NaN; % suppress log(0) warnings


        
% Plot
figure('Color','w','Position',[100 100 1000 600]);

% Subplot 1: Tip velocity
subplot(3,1,1);
plot(t*1e3, v_tip_L*100, 'k', 'LineWidth', 1.2); hold on;
plot(t*1e3, v_tip_R*100, 'r', 'LineWidth', 1.2);
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('Time (ms)');
ylabel('Tip Velocity (cm/s)');
legend('Left Ear', 'Right Ear', 'Location', 'northeast');
title('Ear Tip Velocity Profiles (Angular Motion Derived)');
grid on;

% Subplot 2: Δf
subplot(3,1,2);
plot(t*1e3, delta_f, 'm', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('\Deltaf (Hz)');
title('Binaural Instantaneous Frequency Difference (\Deltaf = f_R - f_L)');
grid on;

% Subplot 3: BLD
subplot(3,1,3);
plot(t*1e3, BLD_dB, 'b', 'LineWidth', 1.3);
xlabel('Time (ms)');
ylabel('BLD (dB)');
title('Instantaneous Binaural Level Difference (Right - Left)');
grid on;
%% Stat
% Ensure column vectors
delta_f_col = delta_f(:);
BLD_dB_col = BLD_dB(:);

% Mask out NaNs (where either is missing or invalid)
valid_idx = isfinite(delta_f_col) & isfinite(BLD_dB_col);

% Compute Pearson correlation
r = corr(delta_f_col(valid_idx), BLD_dB_col(valid_idx));

fprintf("Correlation between Δf and BLD: %.3f\n", r);

%% Optional plot
% figure;
% scatter(delta_f_col(valid_idx), BLD_dB_col(valid_idx), 10, 'filled');
% xlabel('\Deltaf (Hz)');
% ylabel('BLD (dB)');
% title('Δf vs. Binaural Level Difference');
% grid on;

%%
% === Additional Preprocessing ===

% Spectrogram parameters
win = 1024; overlap = 800; nfft = 1024;

% Envelope and BLD
env_left = abs(hilbert(echo_left));
env_right = abs(hilbert(echo_right));
BLD_dB = 20 * log10(env_right ./ env_left);
BLD_dB(~isfinite(BLD_dB)) = NaN;

% Δf
delta_f = instfreq_right_clean - instfreq_left_clean;

% Normalised velocity profiles
v_norm_L = v_tip_L / max(abs(v_tip_L));
v_norm_R = v_tip_R / max(abs(v_tip_R));

% === PLOTTING ===

figure('Color','w','Position',[100 100 1000 800]);

%% Top Subplot: RGB Overlay of Left (green) and Right (red) Spectrograms
% Spectrogram Overlay with White Background (Left = Black, Right = Red Only)
subplot(3,1,1);

% Spectrogram settings
win = 2048; overlap = 1900; nfft = 4096;

[S_L,F,T,P_L] = spectrogram(echo_left, win, overlap, nfft, fs);
[S_R,~,~,P_R] = spectrogram(echo_right, win, overlap, nfft, fs);

% Convert to dB
P_L_dB = 10*log10(abs(P_L));
P_R_dB = 10*log10(abs(P_R));
floor_dB = -90; ceil_dB = -40;

% Clip and normalize
P_L_dB = max(P_L_dB, floor_dB);
P_R_dB = max(P_R_dB, floor_dB);
P_L_norm = (P_L_dB - floor_dB) / (ceil_dB - floor_dB);
P_R_norm = (P_R_dB - floor_dB) / (ceil_dB - floor_dB);

% Init white background
img = ones([size(P_L_norm,1), size(P_L_norm,2), 3]);

% Subtract black (L) from all channels
img(:,:,1) = img(:,:,1) - 0.7 * P_L_norm;  % R
img(:,:,2) = img(:,:,2) - 0.7 * P_L_norm;  % G
img(:,:,3) = img(:,:,3) - 0.7 * P_L_norm;  % B

% Subtract green & blue from R ear to get red
img(:,:,2) = img(:,:,2) - 0.7 * P_R_norm;  % remove green
img(:,:,3) = img(:,:,3) - 0.7 * P_R_norm;  % remove blue

% Clamp values
img = min(max(img, 0), 1);

% Plot
imagesc(T*1e3, F/1e3, img);
axis xy;
ylim([50 95]);
xlabel('Time (ms)');
ylabel('Frequency (kHz)');
title('Spectrograms: Left = Black, Right = Red');
grid on;
ax = gca;
ax.Layer = 'top';

%% Middle Subplot: Echoes + Velocity
subplot(3,1,2);
yyaxis left;
plot(t*1e3, echo_left, 'k'); hold on;
plot(t*1e3, echo_right, 'r');
ylabel('Echo Amplitude');

yyaxis right;
plot(t*1e3, v_norm_L, '--k');
plot(t*1e3, v_norm_R, '--r');
ylabel('Normalised Tip Velocity');

xlabel('Time (ms)');
title('Echo Waveforms \& Ear Tip Velocities');
grid on;

%% Bottom Subplot: Δf and BLD with Coloured Axes
subplot(3,1,3);
yyaxis left
plot(t*1e3, delta_f, 'm', 'LineWidth', 1.2);
ylabel('$\Delta$f (kHz)', 'Color', 'm');
yticklabels(yticks./1000)
set(gca, 'ycolor', 'm');

yyaxis right
plot(t*1e3, BLD_dB, 'b', 'LineWidth', 1.2);
ylabel('BLD (dB)', 'Color', 'b');
set(gca, 'ycolor', 'b');

xlabel('Time (ms)');
title('Binaural Frequency and Level Difference');
grid on;
formatLatexAll(gcf)