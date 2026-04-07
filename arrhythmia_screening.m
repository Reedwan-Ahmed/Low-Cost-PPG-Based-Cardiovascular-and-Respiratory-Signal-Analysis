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
% Step 5: Smooth + pulse band filtering
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
high_cutoff = 5.0;   % Hz

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

peak_times = time_s(peak_locs);

if numel(peak_times) < 4
    error('Too few peaks detected for irregular pulse analysis.');
end

% ============================
% Step 7: Pulse intervals
% ============================
ppi = diff(peak_times);          % seconds
ppi_ms = ppi * 1000;             % milliseconds
ppi_times = peak_times(2:end);

% Gross artifact cleanup only
valid_ppi = ppi > 0.4 & ppi < 1.5;
ppi = ppi(valid_ppi);
ppi_ms = ppi_ms(valid_ppi);
ppi_times = ppi_times(valid_ppi);

if numel(ppi_ms) < 10
    error('Too few valid pulse intervals after cleanup.');
end

% ============================
% Step 8: Poincare plot values
% ============================
ppi_n  = ppi_ms(1:end-1);
ppi_n1 = ppi_ms(2:end);

if numel(ppi_n) < 5
    error('Too few interval pairs for Poincare plot.');
end

% ============================
% Step 9: Irregularity metrics
% ============================
meanPPI   = mean(ppi_ms);
medianPPI = median(ppi_ms);
stdPPI    = std(ppi_ms);
cvPPI     = 100 * stdPPI / meanPPI;   % coefficient of variation %

successive_diff_ms = diff(ppi_ms);
abs_successive_diff_ms = abs(successive_diff_ms);

RMSSD_ms = sqrt(mean(successive_diff_ms.^2));

% Poincare metrics
SD1 = std(ppi_n1 - ppi_n) / sqrt(2);
SD2 = sqrt(2 * std(ppi_ms)^2 - SD1^2);

SD1_SD2_ratio = SD1 / SD2;

% Simple irregular interval percentage
deviation_percent = abs(ppi_ms - medianPPI) / medianPPI * 100;
irregular_count = sum(deviation_percent > 20);
irregular_percent = 100 * irregular_count / numel(ppi_ms);

% Successive large changes
large_jump_count = sum(abs_successive_diff_ms > 80);
large_jump_percent = 100 * large_jump_count / numel(abs_successive_diff_ms);

% ============================
% Step 10: Simple project screening flag
% ============================
screen_flag = false;

if irregular_percent > 10
    screen_flag = true;
end

if large_jump_percent > 10
    screen_flag = true;
end

if cvPPI > 8
    screen_flag = true;
end

% ============================
% Step 11: Display results
% ============================
disp(' ');
disp('Irregular Pulse Screening Metrics');
disp(['Mean PPI              = ', num2str(meanPPI), ' ms']);
disp(['Median PPI            = ', num2str(medianPPI), ' ms']);
disp(['STD of PPI            = ', num2str(stdPPI), ' ms']);
disp(['CV of PPI             = ', num2str(cvPPI), ' %']);
disp(['RMSSD                 = ', num2str(RMSSD_ms), ' ms']);
disp(['SD1                   = ', num2str(SD1), ' ms']);
disp(['SD2                   = ', num2str(SD2), ' ms']);
disp(['SD1/SD2               = ', num2str(SD1_SD2_ratio)]);
disp(['Irregular Intervals   = ', num2str(irregular_percent), ' %']);
disp(['Large Successive Jumps= ', num2str(large_jump_percent), ' %']);

disp(' ');
disp('Project Screening Interpretation');
if screen_flag
    disp('Irregular pulse trend detected in this recording (screening result).');
else
    disp('No strong irregular pulse trend detected in this recording.');
end

% ============================
% Step 12: Plots
% ============================

% Detected peaks
figure;
plot(time_s, ppg_norm, 'm'); hold on;
plot(peak_times, peak_values, 'ko', 'MarkerFaceColor', 'g');
grid on;
xlabel('Time (seconds)');
ylabel('Normalized Amplitude');
title('Normalized Filtered PPG with Detected Peaks');
legend('Filtered PPG', 'Detected Peaks');

% PPI over time
figure;
subplot(2,1,1);
plot(ppi_times, ppi_ms, 'b-o', 'LineWidth', 1);
grid on;
xlabel('Time (seconds)');
ylabel('PPI (ms)');
title('Pulse-to-Pulse Intervals');

subplot(2,1,2);
plot(ppi_times(2:end), abs_successive_diff_ms, 'r-o', 'LineWidth', 1);
grid on;
xlabel('Time (seconds)');
ylabel('|ΔPPI| (ms)');
title('Absolute Successive Pulse Interval Changes');

% Poincare plot
figure;
scatter(ppi_n, ppi_n1, 35, 'filled'); hold on;
min_val = min([ppi_n ppi_n1]);
max_val = max([ppi_n ppi_n1]);
plot([min_val max_val], [min_val max_val], 'k--', 'LineWidth', 1.2);
grid on;
xlabel('PPI_n (ms)');
ylabel('PPI_{n+1} (ms)');
title('Poincare Plot for Irregular Pulse Screening');
legend('Interval Pairs', 'Line of Identity', 'Location', 'best');
axis equal;

% Histogram of PPI
figure;
histogram(ppi_ms, 12);
grid on;
xlabel('PPI (ms)');
ylabel('Count');
title('Distribution of Pulse Intervals');