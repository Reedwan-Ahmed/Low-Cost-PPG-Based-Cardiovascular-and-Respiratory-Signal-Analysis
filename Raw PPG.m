clear;
clc;
close all;

% Select the saved txt file
[file, path] = uigetfile('*.txt', 'Select your PPG text file');

if isequal(file, 0)
    error('No file was selected.');
end

fullFileName = fullfile(path, file);

% Read the text file: column 1 = time in ms, column 2 = IR value
data = readmatrix(fullFileName);

if isempty(data)
    error('The file is empty or could not be read.');
end

if size(data, 2) < 2
    error('The file must have at least 2 columns: time_ms, ir_value');
end

time_ms = data(:, 1);
ir_raw  = data(:, 2);

% Remove rows where finger was not detected (IR = 0)
valid_idx = ir_raw > 0;
time_ms = time_ms(valid_idx);
ir_raw  = ir_raw(valid_idx);

if isempty(ir_raw)
    error('All values are zero. Finger contact was not recorded properly.');
end

% Convert time to seconds and start from zero
time_s = (time_ms - time_ms(1)) / 1000;

% Estimate sampling frequency from timestamp difference
dt = diff(time_s);
Fs = 1 / median(dt);

disp(['Estimated Sampling Frequency = ', num2str(Fs), ' Hz']);

% Plot raw IR signal
figure;
plot(time_s, ir_raw, 'b');
grid on;
xlabel('Time (seconds)');
ylabel('IR Value');
title('Raw PPG Signal from MAX30102');

% Zoomed view for easier inspection
figure;
plot(time_s, ir_raw, 'b');
grid on;
xlabel('Time (seconds)');
ylabel('IR Value');
title('Raw PPG Signal (Zoom View)');
xlim([0 min(10, time_s(end))]);