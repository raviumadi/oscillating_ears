function call = generateCFBatCall(f0, dur_ms, fs, tail, velocity)
% Generate a constant-frequency bat call with optional Doppler shift via resampling
% Inputs:
%   f0        - frequency (Hz)
%   dur_ms    - duration in ms (emitted duration)
%   fs        - sampling rate
%   tail      - % of duration to pad with zeros
%   velocity  - relative velocity (m/s). Negative if approaching reflector (echo resample)

    c = 343;                        % Speed of sound
    dur_s = dur_ms / 1000;         % Convert ms to s
    t = 0:1/fs:dur_s - 1/fs;       % Time vector for original call
    tone = sin(2*pi*f0*t);         % Pure tone

    % Apply Hanning window
    tone = tone .* hann(length(tone))';

    % Doppler resample if velocity is not zero
    if velocity ~= 0
        % Resampling factor from Doppler equation (exact for time & frequency)
        doppler_ratio = sqrt((c - velocity) / (c + velocity));  % < 1 if approaching
        [p, q] = rat(doppler_ratio, 1e-6);  % Rational approx for resample
        tone = resample(tone, p, q);
    end

    % Normalize
    tone = tone / max(abs(tone));

    % Zero-padding
    tail_samples = round(tail / 100 * fs * dur_s);
    call = [zeros(1, tail_samples), tone, zeros(1, tail_samples)]';
end