clear;
clc;
close all;

% ----------------------------
% Step 1: Select and read file
% ----------------------------
[file, path] = uigetfile('*.txt', 'Select your PPG text file');

if isequal(file, 0)
    error('No file was selected.');
end

fullFileName = fullfile(path, file);
data = readmatrix(fullFileName);

if isempty(data)
    error('The file is empty or could not be read.');
end

if size(data, 2) < 2
    error('The file must have 2 columns: time_ms and ir_value');
end

time_ms = data(:, 1);
ir_raw  = data(:, 2);

% ----------------------------------------
% Step 2: Remove zero values (no finger)
% ----------------------------------------
valid_idx = ir_raw > 0;
time_ms = time_ms(valid_idx);
ir_raw  = ir_raw(valid_idx);

if isempty(ir_raw)
    error('All IR values are zero. Recording is not usable.');
end

% ----------------------------------------
% Step 3: Convert time to seconds
% ----------------------------------------
time_s = (time_ms - time_ms(1)) / 1000;

% ----------------------------------------
% Step 4: Estimate sampling frequency
% ----------------------------------------
dt = diff(time_s);
Fs = 1 / median(dt);

disp(['Estimated Sampling Frequency = ', num2str(Fs), ' Hz']);

% ----------------------------------------
% Step 5: Smooth raw PPG slightly
% ----------------------------------------
smooth_window_sec = 0.15;
smooth_window_samples = round(smooth_window_sec * Fs);
if smooth_window_samples < 3
    smooth_window_samples = 3;
end
if mod(smooth_window_samples, 2) == 0
    smooth_window_samples = smooth_window_samples + 1;
end

ir_smooth = movmean(ir_raw, smooth_window_samples);

% ----------------------------------------
% Step 6: Extract slow respiratory component
% Respiration usually ~0.1 to 0.5 Hz
% (6 to 30 breaths/min)
% ----------------------------------------
low_cutoff_resp  = 0.10;   % Hz
high_cutoff_resp = 0.50;   % Hz

[b_resp, a_resp] = butter(3, [low_cutoff_resp high_cutoff_resp] / (Fs/2), 'bandpass');
resp_signal = filtfilt(b_resp, a_resp, ir_smooth);

% ----------------------------------------
% Step 7: Normalize respiratory signal
% ----------------------------------------
resp_norm = (resp_signal - mean(resp_signal)) / std(resp_signal);

% ----------------------------------------
% Step 8: Estimate respiration rate using FFT
% ----------------------------------------
N = length(resp_norm);
Y = fft(resp_norm);
P2 = abs(Y / N);
P1 = P2(1:floor(N/2)+1);
f = Fs * (0:floor(N/2)) / N;

% Restrict to respiration band only
resp_band = (f >= 0.10 & f <= 0.50);

if ~any(resp_band)
    error('No valid respiration frequency range found.');
end

[~, max_idx] = max(P1(resp_band));
f_resp_candidates = f(resp_band);
resp_freq_hz = f_resp_candidates(max_idx);

resp_rate_bpm = resp_freq_hz * 60;

disp(['Estimated Respiration Rate = ', num2str(resp_rate_bpm), ' breaths/min']);

% ----------------------------------------
% Step 9: Optional peak-based respiration estimate
% ----------------------------------------
minBreathDistance_sec = 1.5;   % avoids unrealistically close breaths
minBreathDistance_samples = round(minBreathDistance_sec * Fs);

[resp_peaks, resp_locs] = findpeaks(resp_norm, ...
    'MinPeakDistance', minBreathDistance_samples, ...
    'MinPeakProminence', 0.2);

resp_peak_times = time_s(resp_locs);

if numel(resp_peak_times) >= 2
    breath_intervals = diff(resp_peak_times);
    breath_rate_bpm_peaks = 60 / mean(breath_intervals);
    disp(['Peak-based Respiration Rate = ', num2str(breath_rate_bpm_peaks), ' breaths/min']);
else
    breath_rate_bpm_peaks = NaN;
    disp('Peak-based Respiration Rate could not be estimated reliably.');
end

% ----------------------------------------
% Step 10: Plot raw PPG
% ----------------------------------------
figure;
plot(time_s, ir_raw, 'b');
grid on;
xlabel('Time (seconds)');
ylabel('IR Value');
title('Raw PPG Signal');

% ----------------------------------------
% Step 11: Plot smoothed PPG
% ----------------------------------------
figure;
plot(time_s, ir_raw, 'Color', [0.7 0.7 0.7]); hold on;
plot(time_s, ir_smooth, 'b', 'LineWidth', 1.2);
grid on;
xlabel('Time (seconds)');
ylabel('IR Value');
title('Raw and Smoothed PPG');
legend('Raw PPG', 'Smoothed PPG');

% ----------------------------------------
% Step 12: Plot extracted respiration signal
% ----------------------------------------
figure;
plot(time_s, resp_norm, 'm', 'LineWidth', 1.2); hold on;
if ~isempty(resp_locs)
    plot(resp_peak_times, resp_peaks, 'ko', 'MarkerFaceColor', 'g');
end
grid on;
xlabel('Time (seconds)');
ylabel('Normalized Amplitude');
title('Extracted Respiratory Modulation from PPG');
legend('Respiratory Signal', 'Detected Breaths');

% ----------------------------------------
% Step 13: Plot respiration spectrum
% ----------------------------------------
figure;
plot(f, P1, 'b', 'LineWidth', 1.2); hold on;
xline(resp_freq_hz, '--r', 'LineWidth', 1.2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Respiration Frequency Spectrum');
xlim([0 1.0]);
legend('Spectrum', 'Estimated Respiration Frequency');