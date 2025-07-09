clear; clc;

%% Fixed Parameters
duration = 100; fs = 192e3; cf_freq = 80e3;
bat_velocity = 2; sound_speed = 343; ear_height = 0.02;

ear_freqs = [10, 20, 30, 40, 50];
theta_degs = [15, 30];
save_path = ''; % define your path

% Subplot grid
nrows = numel(ear_freqs);
ncols = numel(theta_degs);

% Create figures
fig_spec = figure('Name','Spectrogram Overlay','Color','w','Position',[100 100 1500 2000]);
fig_echo = figure('Name','Echo and Velocity','Color','w','Position',[100 100 1500 2000]);
fig_delta = figure('Name','Δf and BLD','Color','w','Position',[100 100 1500 2000]);

run_id = 1;

for i = 1:length(ear_freqs)
    for j = 1:length(theta_degs)
        ear_freq = ear_freqs(i);
        theta_max_deg = theta_degs(j);

        %% Derived parameters
        theta_max = deg2rad(theta_max_deg/2);
        omega_max = 2 * pi * ear_freq;
        v_tip = theta_max * omega_max * ear_height;
        s_max = ear_height * deg2rad(theta_max_deg);
        disp_mm = s_max * 1000;

        %% Generate call
        [cf_call] = generateCFBatCall(cf_freq, duration, fs, 0, bat_velocity);
        t = (0:length(cf_call)-1) / fs;

        n_segments = 100;
        segment_positions = linspace(0, ear_height, n_segments);
        segment_velocity = linspace(0, v_tip, n_segments);
        delays = segment_positions / sound_speed;
        dt = 1/fs;

        echo_left = zeros(size(cf_call));
        echo_right = zeros(size(cf_call));

        for seg = 1:n_segments
            delay_samples = round(delays(seg) * fs);
            v_seg = segment_velocity(seg);
            phiL = 0; phiR = 0;

            for k = 1:length(t)
                vL = v_seg * sin(2*pi*ear_freq*t(k));
                vR = -v_seg * sin(2*pi*ear_freq*t(k));
                v_eff_L = bat_velocity + vL;
                v_eff_R = bat_velocity + vR;
                fL = cf_freq * (1 + v_eff_L / sound_speed);
                fR = cf_freq * (1 + v_eff_R / sound_speed);
                phiL = phiL + 2*pi*fL*dt;
                phiR = phiR + 2*pi*fR*dt;
                idxL = k + delay_samples;
                idxR = k + delay_samples;

                if idxL <= length(cf_call)
                    echo_left(idxL) = echo_left(idxL) + sin(phiL);
                end
                if idxR <= length(cf_call)
                    echo_right(idxR) = echo_right(idxR) + sin(phiR);
                end
            end
        end

        echo_left = echo_left / max(abs(echo_left));
        echo_right = echo_right / max(abs(echo_right));

        analytic_left = hilbert(echo_left);
        analytic_right = hilbert(echo_right);

        phase_left = unwrap(angle(analytic_left));
        phase_right = unwrap(angle(analytic_right));
        instfreq_left = [0; diff(phase_left)] * fs / (2 * pi);
        instfreq_right = [0; diff(phase_right)] * fs / (2 * pi);
        env_left = abs(analytic_left);
        env_right = abs(analytic_right);
        threshold = 0.01 * max([env_left; env_right]);
        instfreq_left(env_left < threshold) = NaN;
        instfreq_right(env_right < threshold) = NaN;

        instfreq_left = movmean(instfreq_left, 50, 'omitnan');
        instfreq_right = movmean(instfreq_right, 50, 'omitnan');

        delta_f = instfreq_right - instfreq_left;
        BLD_dB = 20 * log10(env_right ./ env_left);
        BLD_dB(~isfinite(BLD_dB)) = NaN;

        v_tip_L = v_tip * sin(2*pi*ear_freq*t);
        v_tip_R = -v_tip * sin(2*pi*ear_freq*t);
        v_norm_L = v_tip_L / max(abs(v_tip_L));
        v_norm_R = v_tip_R / max(abs(v_tip_R));
        time_ms = t * 1e3;

        % -- Clean Δf: discard noisy tail end (last 20 ms)
        discard_ms = 20;
        delta_f_clean = delta_f(isfinite(delta_f));
        discard_idx = round(length(delta_f_clean) * (discard_ms / duration));
        df_valid = delta_f_clean(discard_idx:end - discard_idx);

        % Compute bandwidth and standard deviation
        bandwidth_df = max(abs(df_valid), [], 'omitnan');

        % ===== Plot 1: Spectrogram Overlay =====
        figure(fig_spec);
        subplot(nrows, ncols, run_id);
        [S_L, F, Tspec, P_L] = spectrogram(echo_left, 2048, 1900, 4096, fs);
        [S_R, ~, ~, P_R] = spectrogram(echo_right, 2048, 1900, 4096, fs);
        P_L_dB = 10 * log10(abs(P_L) + eps);
        P_R_dB = 10 * log10(abs(P_R) + eps);
        floor_dB = -80; ceil_dB = -40;
        P_L_norm = (max(min(P_L_dB, ceil_dB), floor_dB) - floor_dB) / (ceil_dB - floor_dB);
        P_R_norm = (max(min(P_R_dB, ceil_dB), floor_dB) - floor_dB) / (ceil_dB - floor_dB);
        img = ones([size(P_L_norm), 3]);
        img(:,:,1) = img(:,:,1) - P_L_norm;
        img(:,:,2) = img(:,:,2) - P_L_norm - 0.5 * P_R_norm;
        img(:,:,3) = img(:,:,3) - P_R_norm;
        imagesc(Tspec * 1e3, F / 1e3, min(max(img,0),1));
        axis xy; ylim([77 83]);
        title([ num2str(theta_max_deg) '$^o$ ' num2str(ear_freq) ' Hz'], 'Interpreter', 'latex');
        xlabel('ms'); ylabel('kHz');
        formatLatex(gca)

        % ===== Plot 2: Echo + Velocity =====
        figure(fig_echo);

        % Ensure row vectors
        time_ms = time_ms(:)';
        echo_left = echo_left(:)';
        echo_right = echo_right(:)';

        subplot(nrows, ncols, run_id);
        yyaxis left;

        % Patch Left (Blue)
        patch([time_ms, fliplr(time_ms)], ...
            [echo_left, zeros(size(echo_left))], ...
            'blue', 'FaceAlpha', 0.95, 'EdgeColor', 'none'); hold on;

        % Patch Right (Orange)
        patch([time_ms, fliplr(time_ms)], ...
            [echo_right, zeros(size(echo_right))], ...
            [1 0.5 0], 'FaceAlpha', 0.95, 'EdgeColor', 'none');
        ylabel('Amplitude')
        formatLatex(gca)

        yyaxis right;
        plot(time_ms, v_norm_L, 'b--');
        plot(time_ms, v_norm_R, 'Color', [1 0.5 0], 'LineStyle', '--');
        title([num2str(theta_max_deg) '$^o$ ' num2str(ear_freq) '~Hz'...
         ' \& ' '$v_{tip}$ = ' num2str(round(v_tip,2)) '~m/s'])

        formatLatex(gca)

        % ===== Plot 3: Δf and BLD =====
        figure(fig_delta);
        subplot(nrows, ncols, run_id);
        yyaxis left;
        plot(time_ms, delta_f, 'm');
        ylim([-2000 2000])
        ylabel('Hz')
        yyaxis right;
        plot(time_ms, BLD_dB, 'b');
        ylim([-25 25])
        title([num2str(theta_max_deg) '$^o$ ' num2str(ear_freq) '~Hz ' '$\Delta f$=' num2str(round(max(abs(bandwidth_df)))) '~Hz' ' BLD=' num2str(round(max(abs(BLD_dB))),2) '~dB'], 'Interpreter', 'latex');
        xlabel('Time (ms)');
        ylabel('dB re. 20uPa')
        formatLatexYY(gca)

        run_id = run_id + 1;
    end
end

% Add shared titles
figure(fig_spec);
sgtitle('Spectrogram Overlay: Left (Blue), Right (Orange)', 'FontSize', 18, 'Interpreter', 'latex');

figure(fig_echo);
sgtitle('Echo Amplitude and Tip Velocity: Left (Blue), Right (Orange)', 'FontSize', 18, 'Interpreter', 'latex');

figure(fig_delta);
sgtitle('Binaural $\Delta f$ and BLD (Level Difference)', 'FontSize', 18, 'Interpreter', 'latex');

%% Save all figures
savefig(fig_spec, fullfile(save_path, 'all_spectrograms.fig'));
savefig(fig_echo, fullfile(save_path, 'all_echo_velocity.fig'));
savefig(fig_delta, fullfile(save_path, 'all_deltf_bld.fig'));

exportgraphics(fig_spec, fullfile(save_path, 'all_spectrograms.pdf'), 'Resolution', 300, 'Append', false)
exportgraphics(fig_echo, fullfile(save_path, 'all_echo_velocity.pdf'), 'Resolution', 300, 'Append', false)
exportgraphics(fig_delta, fullfile(save_path, 'all_deltf_bld.pdf'), 'Resolution', 300, 'Append', false)





