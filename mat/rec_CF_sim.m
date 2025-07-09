%  play stimulus and record full session
clear; close all;

% Parameters
fs = 192e3;
dur = 0.1; % sec
c = 343;
asyncBuff = dsp.AsyncBuffer(20*fs);  % generous buffer for 10s

% === Generate stimulus with trigger and silence between ===
tone = generateCFBatCall(80e3, dur*1000, fs, 1, 0);
padding = zeros(fs-dur*fs, 2);  % ~1 second silence
trig = zeros(length(tone), 1); trig(1) = 0.6;
block = [tone trig];     % tone + trigger
block = [block; padding]; % add silence
out = repmat(block, 10, 1); % repeat 10 times
audiowrite("DaubCallSim.wav", out, fs); % The name is indicative of the source script developed simulating M. daubentonii calls.

% === Set up audio I/O ===
deviceReader = audioDeviceReader(...
    'SampleRate', fs, ...
    'Device', "Babyface (23461881)", ...
    'SamplesPerFrame', 1024, ...
    'BitDepth', '16-bit integer', ...
    'NumChannels', 2);
setup(deviceReader);

fileReader = dsp.AudioFileReader('DaubCallSim.wav');
fileInfo = audioinfo('DaubCallSim.wav');

deviceWriter = audioDeviceWriter(...
    'Device', "Babyface (23461881)", ...
    'SampleRate', fileReader.SampleRate);
setup(deviceWriter, zeros(fileReader.SamplesPerFrame, fileInfo.NumChannels));
disp("Press any key to start...")
pause
% === Playback and recording loop ===
disp("Playing and recording...");

frq = [30, 25, 20, 15, 10, 5, 0];

for i = 1: length(frq)
    speaker_freq = frq(i);
    disp(['Speaker Frq_' num2str(speaker_freq)])
    pause(8);
    disp('Recording now...')
    while ~isDone(fileReader)
        audioToPlay = fileReader();
        deviceWriter(audioToPlay);
        [audioRecorded, ~] = deviceReader();
        write(asyncBuff, audioRecorded);
    end

    % === Save full recorded signal ===
    audioData = read(asyncBuff);
    audiowrite(['recorded_response_' num2str(speaker_freq) '.wav'], audioData, fs);
    reset(fileReader);
    reset(asyncBuff);
end
disp("Recording complete");

% === Clean up ===
release(fileReader);
release(deviceReader);
release(deviceWriter);