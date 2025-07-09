% ear movement freq is included in this version.
clear; clc;

%% Parameters
duration = 100;              % Duration in milliseconds
fs = 192e3;                     % Sampling rate (Hz)
cf_freq = 80e3;                 % CF bat call frequency (Hz)
bat_velocity = 2;               % Forward flight velocity (m/s)
sound_speed = 343;              % Speed of sound (m/s)
v_base = 0; v_tip = 3;         % Ear velocity profile (max tip speed)
ear_freq = 10;                   % Ear oscillation frequency (Hz)

doppler_factor = (1 + bat_velocity / sound_speed); % Not used directly here

%% 1. Generate CF Call
[cf_call] = generateCFBatCall(cf_freq, duration, fs, 0, bat_velocity);
t = (0:length(cf_call)-1) / fs;
cf_shifted = cf_call;

%% 2. Simulate Time-Varying Echoes for Both Ears with Ear Motion Oscillation

n_segments = 100;
ear_length = 0.03;
segment_positions = linspace(0, ear_length, n_segments);
segment_velocity = linspace(v_base, v_tip, n_segments); 
delays = 2 * segment_positions / sound_speed;

echo_left = zeros(size(cf_call));
echo_right = zeros(size(cf_call));

for i = 1:n_segments
    delay_samples = round(delays(i) * fs);
    
    % Segment velocity weight (scaled per pinna segment)
    v_seg = segment_velocity(i);
    
    for j = 1:length(t)
        % Sinusoidal velocity at time j
        vL = v_seg * sin(2*pi*ear_freq*t(j));   % Left ear
        vR = -v_seg * sin(2*pi*ear_freq*t(j));  % Right ear (antiphase)

        % Doppler-shifted instantaneous frequencies
        fL = cf_freq * (1 + vL / sound_speed);
        fR = cf_freq * (1 + vR / sound_speed);

        % Echo samples with delay
        if j + delay_samples <= length(cf_call)
            echo_left(j + delay_samples) = echo_left(j + delay_samples) + sin(2 * pi * fL * t(j));
            echo_right(j + delay_samples) = echo_right(j + delay_samples) + sin(2 * pi * fR * t(j));
        end
    end
end

% Normalize
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

%% 4. Plot Time-Domain, Cross-Correlation, and Power Spectrum

figure;

subplot(4,1,1);
plot(t*1e3, cf_call, 'k');
title('Original CF Call'); 
xlabel('Time (ms)'); ylabel('Amplitude');

subplot(4,1,2);
plot(t*1e3, echo_left, 'k'); hold on;
plot(t*1e3, echo_right, 'r'); 
title('Echo from Left (black) and Right (red) Ears');
xlabel('Time (ms)'); ylabel('Amplitude');

subplot(4,1,3);
plot(time_axis, corr_L, 'k'); hold on;
plot(time_axis, corr_R, 'r');
title('Time-Varying Cross-Correlation');
xlabel('Time (ms)'); ylabel('Max Correlation');
legend('Left', 'Right');

% subplot(4,1,4);
% pspec(echo_left, fs, 1); hold on;
% pspec(echo_right, fs, 1);
% title('Power Spectrum: Left (black) and Right (red) Echo');
% xlabel('Frequency (Hz)'); ylabel('Power');

sgtitle('Binaural Simulation: Oscillating Ear Motion and Doppler Effects');
%% 7. Clean Full-Scale Plot of Instantaneous Frequency vs CF

% Hilbert analytic signals
analytic_left = hilbert(echo_left);
analytic_right = hilbert(echo_right);

% Instantaneous phase and frequency (Hz)
phase_left = unwrap(angle(analytic_left));
phase_right = unwrap(angle(analytic_right));
instfreq_left = [0; diff(phase_left)] * fs / (2 * pi);
instfreq_right = [0; diff(phase_right)] * fs / (2 * pi);

% Envelope threshold masking
env_left = abs(analytic_left);
env_right = abs(analytic_right);
threshold = 0.1 * max([env_left; env_right]);

% Skip masking — show all data
instfreq_left_clean = movmean(instfreq_left, 50);
instfreq_right_clean = movmean(instfreq_right, 50);

% Time in ms
time_ms = t * 1e3;

% Plot full scale
% figure('Color','w','Position',[100 100 1000 400]);
subplot(4,1,4);
plot(time_ms, instfreq_left_clean, 'k', 'LineWidth', 1.3); hold on;
plot(time_ms, instfreq_right_clean, 'r', 'LineWidth', 1.3);
yline(cf_freq, '--', 'CF Frequency', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
title('Instantaneous Frequency of Echo vs. CF Frequency');
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
legend('Left Ear', 'Right Ear', 'Location', 'northwest');
xlim([0 100]);
ylim([7.0e4 9.0e4]);  % Adjust this based on spread
grid on;

%% 6. Spectrograms of Both Ears

figure;
subplot(2,1,1);
spectrogram(echo_left, 1024, 800, 1024, fs, 'yaxis');
title('Spectrogram of Left Ear Echo');

subplot(2,1,2);
spectrogram(echo_right, 1024, 800, 1024, fs, 'yaxis');
title('Spectrogram of Right Ear Echo');