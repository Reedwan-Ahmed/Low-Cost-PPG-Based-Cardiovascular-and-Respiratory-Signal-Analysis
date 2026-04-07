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
% Step 5: Remove slow baseline drift
% ----------------------------------------
baseline_window_sec = 1.0;
baseline_window_samples = round(baseline_window_sec * Fs);

if mod(baseline_window_samples, 2) == 0
    baseline_window_samples = baseline_window_samples + 1;
end

baseline = movmean(ir_raw, baseline_window_samples);
ir_detrended = ir_raw - baseline;

% ----------------------------------------
% Step 6: Bandpass filter for pulse signal
% ----------------------------------------
low_cutoff  = 0.5;   % Hz
high_cutoff = 5.0;   % Hz

[b, a] = butter(3, [low_cutoff high_cutoff] / (Fs/2), 'bandpass');
ppg_filtered = filtfilt(b, a, ir_detrended);

% ----------------------------------------
% Step 7: Normalize filtered signal
% ----------------------------------------
ppg_norm = (ppg_filtered - mean(ppg_filtered)) / std(ppg_filtered);

% ----------------------------------------
% Step 8: Detect pulse peaks
% ----------------------------------------
minPeakDistance_sec = 0.4;
minPeakDistance_samples = round(minPeakDistance_sec * Fs);

[peak_values, peak_locs] = findpeaks(ppg_norm, ...
    'MinPeakDistance', minPeakDistance_samples, ...
    'MinPeakHeight', 0.2, ...
    'MinPeakProminence', 0.4);

peak_times = time_s(peak_locs);

if numel(peak_times) < 3
    error('Too few peaks detected. Recording is not enough for PRV.');
end

% ----------------------------------------
% Step 9: Pulse intervals and pulse rate
% ----------------------------------------
pulse_intervals = diff(peak_times);      % seconds
pulse_rate_bpm = 60 ./ pulse_intervals;  % BPM

mean_bpm = mean(pulse_rate_bpm);
disp(['Average Pulse Rate = ', num2str(mean_bpm), ' BPM']);

% ----------------------------------------
% Step 10: Remove abnormal intervals (simple cleanup)
% ----------------------------------------
valid_pi = pulse_intervals > 0.4 & pulse_intervals < 1.5;
pulse_intervals_clean = pulse_intervals(valid_pi);
pi_times_clean = peak_times(2:end);
pi_times_clean = pi_times_clean(valid_pi);

if numel(pulse_intervals_clean) < 10
    error('Too few clean pulse intervals for PRV analysis.');
end

% ----------------------------------------
% Step 11: Time-domain PRV features
% ----------------------------------------
SDNN  = std(pulse_intervals_clean);                          % seconds
RMSSD = sqrt(mean(diff(pulse_intervals_clean).^2));         % seconds
MeanPI = mean(pulse_intervals_clean);                       % seconds
MeanHR = 60 / MeanPI;                                       % BPM

disp(' ');
disp('Time-Domain PRV Features');
disp(['Mean Pulse Interval = ', num2str(MeanPI), ' s']);
disp(['Mean Pulse Rate     = ', num2str(MeanHR), ' BPM']);
disp(['SDNN                = ', num2str(SDNN), ' s']);
disp(['RMSSD               = ', num2str(RMSSD), ' s']);

% ----------------------------------------
% Step 12: Interpolate PRV for frequency analysis
% ----------------------------------------
prv_fs = 4;   % Hz, standard resampling rate for interval series
t_prv = pi_times_clean(1):1/prv_fs:pi_times_clean(end);

prv_interp = interp1(pi_times_clean, pulse_intervals_clean, t_prv, 'pchip');

% Remove mean and slow trend
prv_interp_detrended = detrend(prv_interp);

% ----------------------------------------
% Step 13: PSD using Welch method
% ----------------------------------------
window_length = min(length(prv_interp_detrended), 256);

[pxx, f] = pwelch(prv_interp_detrended, window_length, [], [], prv_fs);

% Frequency bands
lf_band = (f >= 0.04 & f < 0.15);
hf_band = (f >= 0.15 & f <= 0.40);

LF_power = trapz(f(lf_band), pxx(lf_band));
HF_power = trapz(f(hf_band), pxx(hf_band));

if HF_power == 0
    LF_HF_ratio = NaN;
else
    LF_HF_ratio = LF_power / HF_power;
end

disp(' ');
disp('Frequency-Domain PRV Features');
disp(['LF Power   = ', num2str(LF_power)]);
disp(['HF Power   = ', num2str(HF_power)]);
disp(['LF/HF      = ', num2str(LF_HF_ratio)]);

% ----------------------------------------
% Step 14: Simple stress interpretation
% ----------------------------------------
disp(' ');
disp('Simple PRV-based Interpretation');

if isnan(LF_HF_ratio)
    disp('LF/HF ratio could not be computed reliably.');
else
    if LF_HF_ratio < 1
        disp('Lower LF/HF: relatively relaxed autonomic trend.');
    elseif LF_HF_ratio >= 1 && LF_HF_ratio < 2
        disp('Moderate LF/HF: balanced to mild stress trend.');
    else
        disp('Higher LF/HF: possible sympathetic dominance / higher stress trend.');
    end
end

% ----------------------------------------
% Step 15: Plots
% ----------------------------------------
figure;
plot(time_s, ppg_norm, 'm'); hold on;
plot(peak_times, peak_values, 'ko', 'MarkerFaceColor', 'g');
grid on;
xlabel('Time (seconds)');
ylabel('Normalized Amplitude');
title('Normalized Filtered PPG with Detected Peaks');
legend('Filtered PPG', 'Detected Peaks');

figure;
subplot(2,1,1);
plot(pi_times_clean, pulse_intervals_clean, 'b-o', 'LineWidth', 1);
grid on;
xlabel('Time (seconds)');
ylabel('Pulse Interval (s)');
title('Clean Pulse Intervals');

subplot(2,1,2);
plot(pi_times_clean, 60 ./ pulse_intervals_clean, 'r-o', 'LineWidth', 1);
grid on;
xlabel('Time (seconds)');
ylabel('Pulse Rate (BPM)');
title('Pulse Rate from Clean Pulse Intervals');

figure;
plot(t_prv, prv_interp_detrended, 'k');
grid on;
xlabel('Time (seconds)');
ylabel('Detrended PRV');
title('Interpolated PRV Signal');

figure;
plot(f, pxx, 'b', 'LineWidth', 1.2); hold on;
xline(0.04, '--k');
xline(0.15, '--r');
xline(0.40, '--g');
grid on;
xlabel('Frequency (Hz)');
ylabel('PSD');
title('PRV Power Spectral Density');
legend('PSD', '0.04 Hz', '0.15 Hz', '0.40 Hz');
xlim([0 0.5]);