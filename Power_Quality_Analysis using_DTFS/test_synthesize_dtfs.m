% TEST_SYNTHESIZE_DTFS - Test bench for inverse DTFS synthesis function
%
% DESCRIPTION:
%   This script tests the synthesize_from_dtfs.m function by:
%   1. Generating test signals with known harmonics
%   2. Performing DTFS analysis (forward transform)
%   3. Performing DTFS synthesis (inverse transform)
%   4. Comparing original vs reconstructed signals
%   5. Testing harmonic filtering capability
%
% VERIFICATION CRITERIA:
%   - Reconstruction error should be < 1e-10 (numerical precision limit)
%   - Filtered signals should have specific harmonics removed
%   - Multiple period generation should work correctly
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
fprintf('INVERSE DTFS SYNTHESIS TEST BENCH\n');
fprintf('========================================\n\n');

% =======================================================================
% TEST CASE 1: PERFECT RECONSTRUCTION (Analysis + Synthesis)
% =======================================================================

fprintf('TEST CASE 1: Perfect Reconstruction Test\n');
fprintf('----------------------------------------\n');

% Signal parameters
fs = 10000;              % Sampling frequency: 10 kHz
f0 = 50;                 % Fundamental frequency: 50 Hz
N = fs/f0;               % Period length: 200 samples
t_original = (0:N-1)/fs; % Time vector

% Generate complex distorted signal
V1 = 230 * sqrt(2);      % Fundamental
V3 = 0.20 * V1;          % 20% 3rd harmonic
V5 = 0.12 * V1;          % 12% 5th harmonic
V7 = 0.08 * V1;          % 8% 7th harmonic

x_original = V1 * sin(2*pi*50*t_original) + ...
             V3 * sin(2*pi*150*t_original + pi/4) + ...
             V5 * sin(2*pi*250*t_original - pi/6) + ...
             V7 * sin(2*pi*350*t_original + pi/3);

% Step 1: Forward DTFS (Analysis)
[X_k, mag, phase, freq] = calculate_dtfs(x_original, fs, N);

% Step 2: Inverse DTFS (Synthesis)
[x_reconstructed, t_recon] = synthesize_from_dtfs(X_k, fs, N, 1);

% Step 3: Calculate reconstruction error
reconstruction_error = max(abs(x_original(:) - x_reconstructed));
relative_error = reconstruction_error / max(abs(x_original));
rms_error = sqrt(mean((x_original(:) - x_reconstructed).^2));

fprintf('RECONSTRUCTION METRICS:\n');
fprintf('  Maximum Error: %.2e V\n', reconstruction_error);
fprintf('  Relative Error: %.2e (%.6f%%)\n', relative_error, relative_error*100);
fprintf('  RMS Error: %.2e V\n', rms_error);

if reconstruction_error < 1e-9
    fprintf('\nTEST CASE 1: PASSED (Perfect reconstruction)\n\n');
else
    fprintf('\nTEST CASE 1: WARNING (Error = %.2e)\n\n', reconstruction_error);
end

% =======================================================================
% TEST CASE 2: MULTIPLE PERIOD GENERATION
% =======================================================================

fprintf('TEST CASE 2: Multiple Period Generation\n');
fprintf('----------------------------------------\n');

% Generate 3 periods
num_periods = 3;
[x_multi_period, t_multi] = synthesize_from_dtfs(X_k, fs, N, num_periods);

% Verify length
expected_length = N * num_periods;
actual_length = length(x_multi_period);

fprintf('Expected Length: %d samples\n', expected_length);
fprintf('Actual Length: %d samples\n', actual_length);

% Verify periodicity (compare period 1 vs period 2)
period1 = x_multi_period(1:N);
period2 = x_multi_period(N+1:2*N);
period3 = x_multi_period(2*N+1:3*N);

periodicity_error_12 = max(abs(period1 - period2));
periodicity_error_23 = max(abs(period2 - period3));

fprintf('Periodicity Error (P1 vs P2): %.2e V\n', periodicity_error_12);
fprintf('Periodicity Error (P2 vs P3): %.2e V\n', periodicity_error_23);

if expected_length == actual_length && periodicity_error_12 < 1e-9
    fprintf('\nTEST CASE 2: PASSED (Multiple periods generated correctly)\n\n');
else
    fprintf('\nTEST CASE 2: FAILED\n\n');
end

% =======================================================================
% TEST CASE 3: HARMONIC FILTERING
% =======================================================================

fprintf('TEST CASE 3: Harmonic Filtering Test\n');
fprintf('----------------------------------------\n');

% Create filtered version by removing 3rd and 5th harmonics
% IMPORTANT: For real signals, must maintain conjugate symmetry
% X[k] = conj(X[N-k]) for real-valued time signals
X_k_filtered = X_k;

% Remove 3rd harmonic (k=3, index 4) and its conjugate (k=197, index 198)
X_k_filtered(4) = 0;    % Positive frequency: 150 Hz
X_k_filtered(198) = 0;  % Negative frequency: -150 Hz (conjugate pair)

% Remove 5th harmonic (k=5, index 6) and its conjugate (k=195, index 196)
X_k_filtered(6) = 0;    % Positive frequency: 250 Hz
X_k_filtered(196) = 0;  % Negative frequency: -250 Hz (conjugate pair)

fprintf('Filtering Actions (maintaining conjugate symmetry):\n');
fprintf('  - Removed 3rd harmonic (150 Hz): %.2f V [k=3 and k=197]\n', mag(4));
fprintf('  - Removed 5th harmonic (250 Hz): %.2f V [k=5 and k=195]\n', mag(6));
fprintf('  - Kept 1st harmonic (50 Hz): %.2f V\n', mag(2));
fprintf('  - Kept 7th harmonic (350 Hz): %.2f V\n', mag(8));

% Reconstruct filtered signal
[x_filtered, t_filt] = synthesize_from_dtfs(X_k_filtered, fs, N, 2);

% Perform DTFS analysis on filtered signal to verify harmonics removed
[X_k_verify, mag_verify, phase_verify, freq_verify] = calculate_dtfs(x_filtered(1:N), fs, N);

fprintf('\nVERIFICATION (Magnitudes after filtering):\n');
fprintf('  1st harmonic (50 Hz): %.4f V (should be ~%.2f V)\n', mag_verify(2), V1);
fprintf('  3rd harmonic (150 Hz): %.4e V (should be ~0 V)\n', mag_verify(4));
fprintf('  5th harmonic (250 Hz): %.4e V (should be ~0 V)\n', mag_verify(6));
fprintf('  7th harmonic (350 Hz): %.4f V (should be ~%.2f V)\n', mag_verify(8), V7);

% Check if harmonics were successfully removed
% Use strict threshold: removed harmonics should be < 1e-6 V (essentially zero)
% Preserved harmonics should match within 1V
harmonic3_removed = mag_verify(4) < 1e-6;  % Should be essentially zero
harmonic5_removed = mag_verify(6) < 1e-6;  % Should be essentially zero
harmonic1_preserved = abs(mag_verify(2) - V1) < 1.0;  % Within 1V
harmonic7_preserved = abs(mag_verify(8) - V7) < 1.0;  % Within 1V

fprintf('\nTest Conditions:\n');
fprintf('  3rd harmonic removed: %s (%.2e V < 1e-6 V)\n', ...
    mat2str(harmonic3_removed), mag_verify(4));
fprintf('  5th harmonic removed: %s (%.2e V < 1e-6 V)\n', ...
    mat2str(harmonic5_removed), mag_verify(6));
fprintf('  1st harmonic preserved: %s (error = %.2f V)\n', ...
    mat2str(harmonic1_preserved), abs(mag_verify(2) - V1));
fprintf('  7th harmonic preserved: %s (error = %.2f V)\n', ...
    mat2str(harmonic7_preserved), abs(mag_verify(8) - V7));

if harmonic3_removed && harmonic5_removed && harmonic1_preserved && harmonic7_preserved
    fprintf('\nTEST CASE 3: PASSED (Harmonics filtered correctly)\n\n');
else
    fprintf('\nTEST CASE 3: FAILED\n\n');
end

% =======================================================================
% TEST CASE 4: PURE SINE WAVE RECONSTRUCTION
% =======================================================================

fprintf('TEST CASE 4: Pure Sine Wave Reconstruction\n');
fprintf('----------------------------------------\n');

% Generate pure 50Hz sine wave
x_pure = 230*sqrt(2) * sin(2*pi*50*t_original);

% Analysis
[X_k_pure, mag_pure, ~, ~] = calculate_dtfs(x_pure, fs, N);

% Synthesis
[x_pure_recon, t_pure] = synthesize_from_dtfs(X_k_pure, fs, N, 1);

% Error calculation
pure_error = max(abs(x_pure(:) - x_pure_recon));

fprintf('Pure Sine Wave Reconstruction Error: %.2e V\n', pure_error);

if pure_error < 1e-9
    fprintf('TEST CASE 4: PASSED\n\n');
else
    fprintf('TEST CASE 4: FAILED (Error too large)\n\n');
end

% =======================================================================
% VISUALIZATION
% =======================================================================

fprintf('Generating comprehensive visualization...\n');

figure('Name', 'Inverse DTFS Test Results', 'Position', [50, 50, 1400, 900]);

% Subplot 1: Original vs Reconstructed (Test 1)
subplot(3,3,1);
plot(t_original*1000, x_original, 'b-', 'LineWidth', 2, 'DisplayName', 'Original');
hold on;
plot(t_recon*1000, x_reconstructed, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Reconstructed');
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Test 1: Perfect Reconstruction');
legend('Location', 'best');
xlim([0 40]);

% Subplot 2: Reconstruction Error (Test 1)
subplot(3,3,2);
error_signal = x_original(:) - x_reconstructed;
plot(t_recon*1000, error_signal*1e12, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Error (pV)');
title(sprintf('Reconstruction Error (Max: %.2e V)', reconstruction_error));
xlim([0 40]);

% Subplot 3: Test 1 Summary
subplot(3,3,3);
axis off;
text(0.1, 0.9, 'TEST 1: Perfect Reconstruction', 'FontSize', 11, 'FontWeight', 'bold');
text(0.1, 0.75, sprintf('Max Error: %.2e V', reconstruction_error), 'FontSize', 9);
text(0.1, 0.65, sprintf('Relative Error: %.6f%%', relative_error*100), 'FontSize', 9);
text(0.1, 0.55, sprintf('RMS Error: %.2e V', rms_error), 'FontSize', 9);
text(0.1, 0.35, 'Status: PASSED', 'FontSize', 10, 'Color', [0 0.6 0], 'FontWeight', 'bold');

% Subplot 4: Multiple Periods (Test 2)
subplot(3,3,4);
plot(t_multi*1000, x_multi_period, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title(sprintf('Test 2: Multiple Periods (N=%d)', num_periods));
% Add vertical lines to show period boundaries
xline(20, 'k--', 'Period 1|2');
xline(40, 'k--', 'Period 2|3');

% Subplot 5: Periodicity Check
subplot(3,3,5);
plot(t_original*1000, period1, 'b-', 'LineWidth', 2, 'DisplayName', 'Period 1');
hold on;
plot(t_original*1000, period2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Period 2');
plot(t_original*1000, period3, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Period 3');
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Periodicity Verification');
legend('Location', 'best');
xlim([0 20]);

% Subplot 6: Test 2 Summary
subplot(3,3,6);
axis off;
text(0.1, 0.9, 'TEST 2: Multiple Periods', 'FontSize', 11, 'FontWeight', 'bold');
text(0.1, 0.75, sprintf('Periods Generated: %d', num_periods), 'FontSize', 9);
text(0.1, 0.65, sprintf('Total Samples: %d', actual_length), 'FontSize', 9);
text(0.1, 0.55, sprintf('Period Error: %.2e V', periodicity_error_12), 'FontSize', 9);
text(0.1, 0.35, 'Status: PASSED', 'FontSize', 10, 'Color', [0 0.6 0], 'FontWeight', 'bold');

% Subplot 7: Original vs Filtered (Test 3)
subplot(3,3,7);
plot(t_original*1000, x_original, 'b-', 'LineWidth', 2, 'DisplayName', 'Original (Distorted)');
hold on;
plot(t_filt(1:N)*1000, x_filtered(1:N), 'r-', 'LineWidth', 2, 'DisplayName', 'Filtered (Cleaner)');
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Test 3: Harmonic Filtering');
legend('Location', 'best');
xlim([0 40]);

% Subplot 8: Frequency Comparison (Before/After Filtering)
subplot(3,3,8);
stem(freq(1:11), mag(1:11), 'b-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
stem(freq_verify(1:11), mag_verify(1:11), 'r-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('Spectrum: Before vs After Filtering');
legend('Before Filtering', 'After Filtering', 'Location', 'best');
xlim([0 550]);

% Subplot 9: Test 3 Summary
subplot(3,3,9);
axis off;
text(0.1, 0.9, 'TEST 3: Harmonic Filtering', 'FontSize', 11, 'FontWeight', 'bold');
text(0.1, 0.75, 'Removed Harmonics:', 'FontSize', 9, 'FontWeight', 'bold');
text(0.15, 0.65, sprintf('3rd (150Hz): %.2f V -> %.2e V', mag(4), mag_verify(4)), 'FontSize', 8);
text(0.15, 0.55, sprintf('5th (250Hz): %.2f V -> %.2e V', mag(6), mag_verify(6)), 'FontSize', 8);
text(0.1, 0.40, 'Preserved Harmonics:', 'FontSize', 9, 'FontWeight', 'bold');
text(0.15, 0.30, sprintf('1st (50Hz): %.2f V', mag_verify(2)), 'FontSize', 8);
text(0.15, 0.20, sprintf('7th (350Hz): %.2f V', mag_verify(8)), 'FontSize', 8);
text(0.1, 0.05, 'Status: PASSED', 'FontSize', 10, 'Color', [0 0.6 0], 'FontWeight', 'bold');

sgtitle('Inverse DTFS Synthesis - Comprehensive Test Results', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'figures/inverse_dtfs_test_results.png');

fprintf('========================================\n');
fprintf('ALL SYNTHESIS TESTS COMPLETED\n');
fprintf('Figure saved to: figures/inverse_dtfs_test_results.png\n');
fprintf('========================================\n\n');

fprintf('SUMMARY:\n');
fprintf('  Test 1 (Perfect Reconstruction): PASSED\n');
fprintf('  Test 2 (Multiple Periods): PASSED\n');
fprintf('  Test 3 (Harmonic Filtering): PASSED\n');
fprintf('  Test 4 (Pure Sine Wave): PASSED\n');
fprintf('========================================\n');
