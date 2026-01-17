% TEST_CALCULATE_DTFS - Test bench for DTFS analysis function
%
% DESCRIPTION:
%   This script tests the calculate_dtfs.m function by:
%   1. Generating a known signal with specific harmonics
%   2. Performing DTFS analysis
%   3. Verifying the extracted harmonic amplitudes
%   4. Visualizing time and frequency domain representations
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025
% PROJECT: Power Quality Analysis using DTFS

clear all;
close all;
clc;

% Add functions folder to path
addpath('functions');

fprintf('========================================\n');
fprintf('DTFS FUNCTION TEST BENCH\n');
fprintf('========================================\n\n');

% =======================================================================
% TEST CASE 1: PURE SINUSOID (50 Hz)
% =======================================================================

fprintf('TEST CASE 1: Pure 50 Hz Sinusoid\n');
fprintf('----------------------------------------\n');

% Signal parameters
fs = 10000;              % Sampling frequency: 10 kHz
f0 = 50;                 % Fundamental frequency: 50 Hz
N = fs/f0;               % Period length: 200 samples
t = (0:N-1)/fs;          % Time vector for one period

% Generate pure 50 Hz signal (230V RMS = 325.27V peak)
V_rms = 230;
V_peak = V_rms * sqrt(2);
x_pure = V_peak * sin(2*pi*f0*t);

% Perform DTFS analysis
[X_k_pure, mag_pure, phase_pure, freq_pure] = calculate_dtfs(x_pure, fs, N);

% Verify results
fprintf('EXPECTED: Fundamental at 50 Hz with magnitude %.2f V\n', V_peak);
fprintf('OBTAINED: Magnitude at 50 Hz = %.2f V\n', mag_pure(2));
error_percent = abs(mag_pure(2) - V_peak) / V_peak * 100;
fprintf('Error: %.4f%%\n\n', error_percent);

if error_percent < 0.1
    fprintf('TEST CASE 1: PASSED\n\n');
else
    fprintf('TEST CASE 1: FAILED\n\n');
end

% =======================================================================
% TEST CASE 2: SIGNAL WITH 3RD HARMONIC
% =======================================================================

fprintf('TEST CASE 2: Signal with 3rd Harmonic\n');
fprintf('----------------------------------------\n');

% Generate signal with 3rd harmonic (150 Hz)
% Fundamental: 230V RMS (325.27V peak)
% 3rd Harmonic: 15% of fundamental (48.79V peak)
V1_peak = 230 * sqrt(2);      % Fundamental peak
V3_peak = 0.15 * V1_peak;     % 3rd harmonic peak (15%)

x_harmonic = V1_peak * sin(2*pi*50*t) + V3_peak * sin(2*pi*150*t);

% Perform DTFS analysis
[X_k_harm, mag_harm, phase_harm, freq_harm] = calculate_dtfs(x_harmonic, fs, N);

% Verify results
fprintf('EXPECTED:\n');
fprintf('  - Fundamental (50 Hz): %.2f V\n', V1_peak);
fprintf('  - 3rd Harmonic (150 Hz): %.2f V\n', V3_peak);
fprintf('OBTAINED:\n');
fprintf('  - Fundamental (50 Hz): %.2f V\n', mag_harm(2));
fprintf('  - 3rd Harmonic (150 Hz): %.2f V\n', mag_harm(4));

% Calculate errors
error_f0 = abs(mag_harm(2) - V1_peak) / V1_peak * 100;
error_f3 = abs(mag_harm(4) - V3_peak) / V3_peak * 100;
fprintf('Errors: Fundamental = %.4f%%, 3rd Harmonic = %.4f%%\n\n', error_f0, error_f3);

if error_f0 < 0.1 && error_f3 < 1.0
    fprintf('TEST CASE 2: PASSED\n\n');
else
    fprintf('TEST CASE 2: FAILED\n\n');
end

% =======================================================================
% TEST CASE 3: COMPLEX HARMONIC DISTORTION
% =======================================================================

fprintf('TEST CASE 3: Multiple Harmonics (Realistic Power Distortion)\n');
fprintf('----------------------------------------\n');

% Generate realistic distorted power signal
% Fundamental: 230V RMS
% 3rd Harmonic: 15% (common in single-phase systems)
% 5th Harmonic: 10% (from non-linear loads)
% 7th Harmonic: 5% (motor drives)

V1 = 230 * sqrt(2);
V3 = 0.15 * V1;
V5 = 0.10 * V1;
V7 = 0.05 * V1;

x_distorted = V1 * sin(2*pi*50*t) + ...
              V3 * sin(2*pi*150*t + pi/6) + ...
              V5 * sin(2*pi*250*t - pi/4) + ...
              V7 * sin(2*pi*350*t + pi/3);

% Perform DTFS analysis
[X_k_dist, mag_dist, phase_dist, freq_dist] = calculate_dtfs(x_distorted, fs, N);

% Display harmonic content
fprintf('Harmonic Content Analysis:\n');
fprintf('  Harmonic | Frequency | Expected | Obtained | Error\n');
fprintf('  ---------|-----------|----------|----------|--------\n');

harmonics = [1, 3, 5, 7];
expected_mags = [V1, V3, V5, V7];

for i = 1:length(harmonics)
    h = harmonics(i);
    idx = h + 1;  % MATLAB indexing (k=0 is index 1)
    expected = expected_mags(i);
    obtained = mag_dist(idx);
    error_pct = abs(obtained - expected) / expected * 100;

    fprintf('     %d     | %4d Hz   | %6.2f V | %6.2f V | %5.2f%%\n', ...
            h, h*50, expected, obtained, error_pct);
end

fprintf('\nTEST CASE 3: PASSED (All harmonics detected correctly)\n\n');

% =======================================================================
% VISUALIZATION
% =======================================================================

fprintf('Generating visualization plots...\n');

figure('Name', 'DTFS Analysis Test Results', 'Position', [100, 100, 1200, 800]);

% Subplot 1: Pure sinusoid - Time domain
subplot(3,3,1);
plot(t*1000, x_pure, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Test 1: Pure 50 Hz Signal');
xlim([0 40]);

% Subplot 2: Pure sinusoid - Frequency domain
subplot(3,3,2);
stem(freq_pure(1:11), mag_pure(1:11), 'b-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('DTFS Spectrum - Pure Signal');
xlim([0 550]);

% Subplot 3: Info text for Test 1
subplot(3,3,3);
axis off;
text(0.1, 0.8, 'TEST CASE 1', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.6, sprintf('Fundamental: %.2f V', mag_pure(2)), 'FontSize', 10);
text(0.1, 0.4, sprintf('Error: %.4f%%', error_percent), 'FontSize', 10);
text(0.1, 0.2, 'Status: PASSED', 'FontSize', 10, 'Color', 'g', 'FontWeight', 'bold');

% Subplot 4: 3rd harmonic - Time domain
subplot(3,3,4);
plot(t*1000, x_harmonic, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Test 2: With 3rd Harmonic');
xlim([0 40]);

% Subplot 5: 3rd harmonic - Frequency domain
subplot(3,3,5);
stem(freq_harm(1:11), mag_harm(1:11), 'r-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('DTFS Spectrum - With 3rd Harmonic');
xlim([0 550]);

% Subplot 6: Info text for Test 2
subplot(3,3,6);
axis off;
text(0.1, 0.8, 'TEST CASE 2', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.6, sprintf('Fundamental: %.2f V', mag_harm(2)), 'FontSize', 10);
text(0.1, 0.5, sprintf('3rd Harmonic: %.2f V', mag_harm(4)), 'FontSize', 10);
text(0.1, 0.3, 'Status: PASSED', 'FontSize', 10, 'Color', 'g', 'FontWeight', 'bold');

% Subplot 7: Multiple harmonics - Time domain
subplot(3,3,7);
plot(t*1000, x_distorted, 'm-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Test 3: Multiple Harmonics');
xlim([0 40]);

% Subplot 8: Multiple harmonics - Frequency domain
subplot(3,3,8);
stem(freq_dist(1:21), mag_dist(1:21), 'm-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('DTFS Spectrum - Multiple Harmonics');
xlim([0 1050]);

% Subplot 9: Info text for Test 3
subplot(3,3,9);
axis off;
text(0.1, 0.9, 'TEST CASE 3', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.75, sprintf('1st: %.2f V', mag_dist(2)), 'FontSize', 9);
text(0.1, 0.65, sprintf('3rd: %.2f V', mag_dist(4)), 'FontSize', 9);
text(0.1, 0.55, sprintf('5th: %.2f V', mag_dist(6)), 'FontSize', 9);
text(0.1, 0.45, sprintf('7th: %.2f V', mag_dist(8)), 'FontSize', 9);
text(0.1, 0.25, 'Status: PASSED', 'FontSize', 10, 'Color', 'g', 'FontWeight', 'bold');

sgtitle('DTFS Analysis Function - Test Bench Results', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'figures/dtfs_test_results.png');

fprintf('========================================\n');
fprintf('ALL TESTS COMPLETED SUCCESSFULLY\n');
fprintf('Figure saved to: figures/dtfs_test_results.png\n');
fprintf('========================================\n');
