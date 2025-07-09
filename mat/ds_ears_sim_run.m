clear; clc;

%% Fixed Parameters
duration = 100;        % ms
fs = 192e3;
cf_freq = 80e3;
bat_velocity = 2;
sound_speed = 343;
ear_height = 0.02;

% Parameter ranges
ear_freqs = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50];
theta_degs = [15, 30, 45];

save_path = ''; % define your path

run_id = 1;

for ear_freq = ear_freqs
    for theta_max_deg = theta_degs

        %% Derived parameters
        theta_max = deg2rad(theta_max_deg/2);
        omega_max = 2 * pi * ear_freq;
        v_tip = theta_max * omega_max * ear_height;
        fprintf('Run %02d: freq = %d Hz, angle = %d°, tip velocity = %.2f m/s\n', ...
            run_id, ear_freq, theta_max_deg, v_tip);
        % Compute one-way and peak-to-peak tip displacement in mm
        s_max = ear_height * deg2rad(theta_max_deg);        % max arc displacement (one-way, in meters)
        disp_mm = s_max * 1000;                             % convert to mm

        %% Generate CF Call
        [cf_call] = generateCFBatCall(cf_freq, duration, fs, 0, bat_velocity);
        t = (0:length(cf_call)-1) / fs;

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

        %% Hilbert transform: Instantaneous frequency
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

        %% === PLOT Final Summary ===
        fig = figure('Color','w','Position',[100 100 1000 800]);

        % Spectrogram Overlay
        subplot(3,1,1);
        % Compute both spectrograms
        [S_L, F, Tspec, P_L] = spectrogram(echo_left, 2048, 1900, 4096, fs);
        [S_R, ~, ~, P_R] = spectrogram(echo_right, 2048, 1900, 4096, fs);

        % Convert to dB
        P_L_dB = 10 * log10(abs(P_L) + eps);
        P_R_dB = 10 * log10(abs(P_R) + eps);

        % Clip dynamic range
        floor_dB = -80; ceil_dB = -40;
        P_L_dB = max(min(P_L_dB, ceil_dB), floor_dB);
        P_R_dB = max(min(P_R_dB, ceil_dB), floor_dB);

        % Normalize to [0,1]
        P_L_norm = (P_L_dB - floor_dB) / (ceil_dB - floor_dB);
        P_R_norm = (P_R_dB - floor_dB) / (ceil_dB - floor_dB);

        % Define white base background
        img = ones([size(P_L_norm), 3]);  % white background = [1 1 1]

        % Overlays
        % Subtract blue for left ear
        img(:,:,1) = img(:,:,1) - P_L_norm;           % Reduce red
        img(:,:,2) = img(:,:,2) - P_L_norm;           % Reduce green

        % Subtract blue & green for right ear (to make orange = [1 0.5 0])
        img(:,:,2) = img(:,:,2) - 0.5 * P_R_norm;     % Reduce green
        img(:,:,3) = img(:,:,3) - P_R_norm;           % Reduce blue

        % Plot RGB spectrogram image
        subplot(3,1,1);  % ensure we’re in the right subplot
        imagesc(Tspec * 1e3, F / 1e3, img);
        axis xy;
        ylim([77 83]);  % adjust to your frequency range
        xlabel('Time (ms)');
        ylabel('Frequency (kHz)');
        title('Spectrogram Overlay: Left (Blue), Right (Orange)');
        formatLatex(gca);

        % Clamp
        img = min(max(img, 0), 1);

        % Plot RGB image
        % Ensure row vectors
        time_ms = time_ms(:)';
        echo_left = echo_left(:)';
        echo_right = echo_right(:)';

        subplot(3,1,2);
        yyaxis left;

        % Patch Left (Blue)
        patch([time_ms, fliplr(time_ms)], ...
            [echo_left, zeros(size(echo_left))], ...
            'blue', 'FaceAlpha', 0.9, 'EdgeColor', 'none'); hold on;

        % Patch Right (Orange)
        patch([time_ms, fliplr(time_ms)], ...
            [echo_right, zeros(size(echo_right))], ...
            [1 0.5 0], 'FaceAlpha', 0.9, 'EdgeColor', 'none');

        ylabel('Echo Amp.');

        yyaxis right;
        plot(time_ms, v_norm_L, 'Color', [0, 0, 1], 'LineWidth', 1.2); hold on
        plot(time_ms, v_norm_R, 'Color', [1, 0.5, 0], 'LineStyle', '-', 'LineWidth', 1.2);
        ylabel('Norm. Velocity');
        xlabel('Time (ms)');
        title(['Synth. Echo and Velocity -- ' num2str(ear_freq) '~Hz, '...
            num2str(theta_max_deg) '$^o$ \& ' '$v_{tip}$ = ' num2str(round(v_tip,2)) '~m/s']);
        formatLatexYY(gca)

        % Δf and BLD
        subplot(3,1,3);
        yyaxis left;
        plot(time_ms, delta_f, 'm', 'LineWidth', 1.2);

        % Clean data
        delta_f_clean = delta_f(isfinite(delta_f));
        ylim([-2000 2000]);
        ylabel('$\Delta$f (Hz)', 'Color', 'm');
        set(gca, 'ycolor', 'm');

        yyaxis right;
        plot(time_ms, BLD_dB, 'b', 'LineWidth', 1.2);
        ylabel('BLD (dB)', 'Color', 'b');
        set(gca, 'ycolor', 'b');
        xlabel('Time (ms)');
        title('Binaural $\Delta$f and BLD');
        formatLatexYY(gca)

        % Save Figure
        filename = sprintf('%s/binaural_run%02d_%dhz_%ddeg.fig', ...
            save_path, run_id, ear_freq, theta_max_deg);
        savefig(fig, filename);

        close(fig);
        run_id = run_id + 1;

        % -- Clean Δf: discard noisy tail end (last 20 ms)
        discard_ms = 20;
        discard_idx = round(length(delta_f_clean) * (discard_ms / duration));
        df_valid = delta_f_clean(discard_idx:end - discard_idx);

        % Compute bandwidth and standard deviation
        bandwidth_df = max(abs(df_valid), [], 'omitnan');
        std_df = std(df_valid, 'omitnan');
        std_bld = std(BLD_dB, 'omitnan');

        % Save row into struct
        summary(run_id).Run = run_id;
        summary(run_id).EarFreqHz = ear_freq;
        summary(run_id).ThetaDeg = theta_max_deg;
        summary(run_id).TipVel = v_tip;
        summary(run_id).TipDisp_mm = disp_mm;
        summary(run_id).DeltaF_BW = bandwidth_df;
        summary(run_id).DeltaF_Std = std_df;
        summary(run_id).BLD_Std = std_bld;
    end
end

% === Save Summary Table ===
summaryTable = struct2table(summary);
disp(summaryTable)

% Save to file
writetable(summaryTable, 'binaural_summary.csv');
save(fullfile(save_path, 'binaural_summary.mat'), 'summaryTable');