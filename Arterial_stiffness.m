clear;
clc;
close all;

% ============================
% Step 1: Select and read file
% ============================
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

% ==================================
% Step 2: Remove zero values/no finger
% ==================================
valid_idx = ir_raw > 0;
time_ms = time_ms(valid_idx);
ir_raw  = ir_raw(valid_idx);

if isempty(ir_raw)
    error('All IR values are zero. Recording is not usable.');
end

% ============================
% Step 3: Convert time to sec
% ============================
time_s = (time_ms - time_ms(1)) / 1000;

% ============================
% Step 4: Estimate Fs
% ============================
dt = diff(time_s);
Fs = 1 / median(dt);

disp(['Estimated Sampling Frequency = ', num2str(Fs), ' Hz']);

% ==========================================
% Step 5: Slight smoothing + pulse band filter
% ==========================================
smooth_window_sec = 0.10;
smooth_window_samples = round(smooth_window_sec * Fs);

if smooth_window_samples < 3
    smooth_window_samples = 3;
end
if mod(smooth_window_samples, 2) == 0
    smooth_window_samples = smooth_window_samples + 1;
end

ir_smooth = movmean(ir_raw, smooth_window_samples);

low_cutoff  = 0.5;   % Hz
high_cutoff = 8.0;   % Hz

[b, a] = butter(3, [low_cutoff high_cutoff] / (Fs/2), 'bandpass');
ppg_filtered = filtfilt(b, a, ir_smooth);

ppg_norm = (ppg_filtered - mean(ppg_filtered)) / std(ppg_filtered);

% ============================
% Step 6: Detect pulse peaks
% ============================
minPeakDistance_sec = 0.4;
minPeakDistance_samples = round(minPeakDistance_sec * Fs);

[peak_values, peak_locs] = findpeaks(ppg_norm, ...
    'MinPeakDistance', minPeakDistance_samples, ...
    'MinPeakHeight', 0.2, ...
    'MinPeakProminence', 0.4);

if numel(peak_locs) < 5
    error('Too few peaks detected for averaged pulse analysis.');
end

% ==========================================
% Step 7: Extract beat-to-beat pulses
% from one local minimum to next local minimum
% ==========================================
M = 200;   % resampled points per pulse

pulse_matrix = [];
pulse_durations = [];

for k = 2:(numel(peak_locs)-1)
    prev_peak = peak_locs(k-1);
    curr_peak = peak_locs(k);
    next_peak = peak_locs(k+1);

    % local minimum before current peak
    [~, prev_min_rel] = min(ir_smooth(prev_peak:curr_peak));
    start_idx = prev_peak + prev_min_rel - 1;

    % local minimum after current peak
    [~, next_min_rel] = min(ir_smooth(curr_peak:next_peak));
    end_idx = curr_peak + next_min_rel - 1;

    if end_idx <= start_idx
        continue;
    end

    pulse_duration = time_s(end_idx) - time_s(start_idx);

    % keep only physiologically reasonable pulse durations
    if pulse_duration < 0.4 || pulse_duration > 1.5
        continue;
    end

    pulse = ir_smooth(start_idx:end_idx);

    % normalize each pulse from 0 to 1
    if max(pulse) == min(pulse)
        continue;
    end
    pulse_norm_single = (pulse - min(pulse)) / (max(pulse) - min(pulse));

    % original local pulse time
    t_local = linspace(0, pulse_duration, numel(pulse_norm_single));
    t_uniform = linspace(0, pulse_duration, M);

    % resample to fixed length
    pulse_resampled = interp1(t_local, pulse_norm_single, t_uniform, 'pchip');

    pulse_matrix = [pulse_matrix; pulse_resampled];
    pulse_durations = [pulse_durations; pulse_duration];
end

if isempty(pulse_matrix)
    error('No valid pulses extracted for averaging.');
end

% ============================
% Step 8: Average pulse
% ============================
avg_pulse = mean(pulse_matrix, 1);
avg_pulse = movmean(avg_pulse, 7);

mean_pulse_duration = mean(pulse_durations);
t_avg = linspace(0, mean_pulse_duration, M);
dt_avg = median(diff(t_avg));

% ============================
% Step 9: Derivatives
% ============================
dppg  = gradient(avg_pulse, dt_avg);
sdppg = gradient(dppg, dt_avg);

% optional small smoothing for sdPPG
sdppg = movmean(sdppg, 5);

% ============================
% Step 10: Detect candidate a,b,c,d,e
% ============================
N = numel(sdppg);

% search windows by fraction of pulse length
idx_a_range = 1 : round(0.25*N);
[a_amp, a_rel] = max(sdppg(idx_a_range));
a_idx = idx_a_range(a_rel);

idx_b_range = a_idx : round(0.45*N);
[b_amp_neg, b_rel] = min(sdppg(idx_b_range));
b_idx = idx_b_range(b_rel);
b_amp = b_amp_neg;

idx_c_range = b_idx : round(0.65*N);
[c_amp, c_rel] = max(sdppg(idx_c_range));
c_idx = idx_c_range(c_rel);

idx_d_range = c_idx : round(0.85*N);
[d_amp_neg, d_rel] = min(sdppg(idx_d_range));
d_idx = idx_d_range(d_rel);
d_amp = d_amp_neg;

idx_e_range = d_idx : N;
[e_amp, e_rel] = max(sdppg(idx_e_range));
e_idx = idx_e_range(e_rel);

% ============================
% Step 11: Ratios and aging index
% ============================
if abs(a_amp) < eps
    error('a-wave amplitude is too small for ratio calculation.');
end

b_a = b_amp / a_amp;
c_a = c_amp / a_amp;
d_a = d_amp / a_amp;
e_a = e_amp / a_amp;

AGI = (b_amp - c_amp - d_amp - e_amp) / a_amp;

disp(' ');
disp('Second-Derivative PPG Candidate Points');
disp(['a = ', num2str(a_amp)]);
disp(['b = ', num2str(b_amp)]);
disp(['c = ', num2str(c_amp)]);
disp(['d = ', num2str(d_amp)]);
disp(['e = ', num2str(e_amp)]);

disp(' ');
disp('Normalized Ratios');
disp(['b/a = ', num2str(b_a)]);
disp(['c/a = ', num2str(c_a)]);
disp(['d/a = ', num2str(d_a)]);
disp(['e/a = ', num2str(e_a)]);

disp(' ');
disp(['Aging Index (AGI) = ', num2str(AGI)]);

% ============================
% Step 12: Plots
% ============================

% all extracted pulses + average pulse
figure;
plot(t_avg, pulse_matrix', 'Color', [0.8 0.8 0.8]); hold on;
plot(t_avg, avg_pulse, 'b', 'LineWidth', 2);
grid on;
xlabel('Time (seconds)');
ylabel('Normalized Amplitude');
title('Extracted Pulses and Averaged Pulse');
legend('Individual Pulses', 'Averaged Pulse');

% averaged pulse + derivatives
figure;
subplot(3,1,1);
plot(t_avg, avg_pulse, 'b', 'LineWidth', 1.8);
grid on;
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Averaged PPG Pulse');

subplot(3,1,2);
plot(t_avg, dppg, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time (seconds)');
ylabel('dPPG');
title('First Derivative of Averaged PPG');

subplot(3,1,3);
plot(t_avg, sdppg, 'k', 'LineWidth', 1.5); hold on;
plot(t_avg(a_idx), a_amp, 'ro', 'MarkerFaceColor', 'r');
plot(t_avg(b_idx), b_amp, 'bo', 'MarkerFaceColor', 'b');
plot(t_avg(c_idx), c_amp, 'go', 'MarkerFaceColor', 'g');
plot(t_avg(d_idx), d_amp, 'mo', 'MarkerFaceColor', 'm');
plot(t_avg(e_idx), e_amp, 'co', 'MarkerFaceColor', 'c');
grid on;
xlabel('Time (seconds)');
ylabel('sdPPG');
title('Second Derivative of Averaged PPG with a,b,c,d,e');
legend('sdPPG','a','b','c','d','e');