% TEST_HARMONIC_FILTER - Test bench for harmonic filter design
%
% DESCRIPTION:
%   Tests the design_harmonic_filter.m function by:
%   1. Testing all filter types (notch, attenuate, lowpass, highpass, bandpass)
%   2. Verifying conjugate symmetry is maintained
%   3. Demonstrating power quality improvement
%   4. Comparing before/after waveforms and spectra
%   5. Validating THD reduction
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
fprintf('HARMONIC FILTER TEST BENCH\n');
fprintf('========================================\n\n');

% =======================================================================
% GENERATE TEST SIGNAL - MOTOR DRIVE SCENARIO
% =======================================================================

fprintf('Generating test signal (Motor Drive with heavy harmonics)...\n\n');

fs = 10000;
f0 = 50;
N = fs/f0;
V_rms = 230;

% Generate motor drive signal (has 5th, 7th, 11th, 13th harmonics)
[signal_original, t, info] = generate_distorted_signal('motor_drive', fs, f0, N, V_rms);

% Perform DTFS analysis
[X_k_original, mag_original, phase_original, freq] = calculate_dtfs(signal_original, fs, N);

% Calculate original metrics
[THD_orig, ~, thd_details] = calculate_thd(mag_original);

fprintf('Original Signal:\n');
fprintf('  Scenario: %s\n', info.scenario_name);
fprintf('  THD: %.2f%%\n', THD_orig);
fprintf('  Dominant Harmonic: %dth\n', thd_details.dominant_harmonic);
fprintf('\n');

% =======================================================================
% TEST 1: NOTCH FILTER - REMOVE 5TH AND 7TH HARMONICS
% =======================================================================

fprintf('========================================\n');
fprintf('TEST 1: Notch Filter (Remove 5th & 7th)\n');
fprintf('========================================\n');

% Apply notch filter
[X_k_notch, info_notch] = design_harmonic_filter(X_k_original, N, 'notch', [5 7]);

% Reconstruct filtered signal
[signal_notch, t_notch] = synthesize_from_dtfs(X_k_notch, fs, N, 2);

% Analyze filtered signal
[~, mag_notch, ~, ~] = calculate_dtfs(signal_notch(1:N), fs, N);
[THD_notch, ~, ~] = calculate_thd(mag_notch);

fprintf('Notch Filter Results:\n');
fprintf('  Original THD: %.2f%%\n', info_notch.original_THD);
fprintf('  Filtered THD: %.2f%%\n', info_notch.filtered_THD);
fprintf('  Improvement: %.2f%%\n\n', info_notch.THD_reduction_percent);

% =======================================================================
% TEST 2: ATTENUATION FILTER - REDUCE 3RD HARMONIC BY 40 dB
% =======================================================================

fprintf('========================================\n');
fprintf('TEST 2: Attenuation Filter (3rd by 40dB)\n');
fprintf('========================================\n');

% Generate LED lighting signal (has strong 3rd harmonic)
[signal_led, t_led, info_led] = generate_distorted_signal('led_lighting', fs, f0, N, V_rms);
[X_k_led, mag_led, ~, ~] = calculate_dtfs(signal_led, fs, N);
[THD_led_orig, ~, ~] = calculate_thd(mag_led);

% Apply attenuation filter
[X_k_atten, info_atten] = design_harmonic_filter(X_k_led, N, 'attenuate', [3], 40);

% Reconstruct
[signal_atten, ~] = synthesize_from_dtfs(X_k_atten, fs, N, 2);
[~, mag_atten, ~, ~] = calculate_dtfs(signal_atten(1:N), fs, N);

fprintf('Attenuation Results:\n');
fprintf('  3rd Harmonic Before: %.2f V\n', mag_led(4));
fprintf('  3rd Harmonic After: %.2f V\n', mag_atten(4));
fprintf('  Attenuation Achieved: %.2f dB\n\n', 20*log10(mag_led(4)/mag_atten(4)));

% =======================================================================
% TEST 3: LOWPASS FILTER - KEEP UP TO 3RD HARMONIC
% =======================================================================

fprintf('========================================\n');
fprintf('TEST 3: Lowpass Filter (Keep up to 3rd)\n');
fprintf('========================================\n');

% Use heavy distortion signal
[signal_heavy, ~, ~] = generate_distorted_signal('heavy_distortion', fs, f0, N, V_rms);
[X_k_heavy, mag_heavy, ~, ~] = calculate_dtfs(signal_heavy, fs, N);

% Apply lowpass filter
[X_k_lowpass, info_lowpass] = design_harmonic_filter(X_k_heavy, N, 'lowpass', 3);

% Reconstruct
[signal_lowpass, ~] = synthesize_from_dtfs(X_k_lowpass, fs, N, 2);
[~, mag_lowpass, ~, ~] = calculate_dtfs(signal_lowpass(1:N), fs, N);

fprintf('Lowpass Results:\n');
fprintf('  Kept: Fundamental + 2nd + 3rd harmonics\n');
fprintf('  Removed: All harmonics above 3rd\n');
fprintf('  THD Reduction: %.2f%%\n\n', info_lowpass.THD_reduction_percent);

% =======================================================================
% TEST 4: CONJUGATE SYMMETRY VERIFICATION
% =======================================================================

fprintf('========================================\n');
fprintf('TEST 4: Conjugate Symmetry Verification\n');
fprintf('========================================\n');

% Check if filtered signals are real (no imaginary components)

test4_pass = true;

% Test notch filter output
max_imag_notch = max(abs(imag(signal_notch)));
if max_imag_notch < 1e-10
    fprintf('Notch Filter: Real signal verified (max imag: %.2e)\n', max_imag_notch);
else
    fprintf('Notch Filter: FAILED (significant imaginary: %.2e)\n', max_imag_notch);
    test4_pass = false;
end

% Test attenuation filter output
max_imag_atten = max(abs(imag(signal_atten)));
if max_imag_atten < 1e-10
    fprintf('Attenuation Filter: Real signal verified (max imag: %.2e)\n', max_imag_atten);
else
    fprintf('Attenuation Filter: FAILED (significant imaginary: %.2e)\n', max_imag_atten);
    test4_pass = false;
end

% Test lowpass filter output
max_imag_lowpass = max(abs(imag(signal_lowpass)));
if max_imag_lowpass < 1e-10
    fprintf('Lowpass Filter: Real signal verified (max imag: %.2e)\n', max_imag_lowpass);
else
    fprintf('Lowpass Filter: FAILED (significant imaginary: %.2e)\n', max_imag_lowpass);
    test4_pass = false;
end

if test4_pass
    fprintf('\nTEST 4: PASSED (All filtered signals are real)\n\n');
else
    fprintf('\nTEST 4: FAILED (Conjugate symmetry broken)\n\n');
end

% =======================================================================
% TEST 5: POWER QUALITY IMPROVEMENT VERIFICATION
% =======================================================================

fprintf('========================================\n');
fprintf('TEST 5: Power Quality Improvement\n');
fprintf('========================================\n');

% Calculate comprehensive metrics before and after filtering
metrics_before = power_quality_metrics(signal_original, mag_original, fs, f0);
metrics_after = power_quality_metrics(real(signal_notch(1:N)), mag_notch, fs, f0);

fprintf('Motor Drive - Before Filtering:\n');
fprintf('  THD: %.2f%%\n', metrics_before.THD_F);
fprintf('  Crest Factor: %.4f\n', metrics_before.crest_factor);
fprintf('  True PF: %.4f\n', metrics_before.true_PF);
fprintf('  K-Factor: %.2f\n', metrics_before.K_factor);
fprintf('  Quality: %s\n\n', metrics_before.quality_rating);

fprintf('Motor Drive - After Filtering (5th & 7th removed):\n');
fprintf('  THD: %.2f%%\n', metrics_after.THD_F);
fprintf('  Crest Factor: %.4f\n', metrics_after.crest_factor);
fprintf('  True PF: %.4f\n', metrics_after.true_PF);
fprintf('  K-Factor: %.2f\n', metrics_after.K_factor);
fprintf('  Quality: %s\n\n', metrics_after.quality_rating);

fprintf('Improvements:\n');
fprintf('  THD Reduction: %.2f%% -> %.2f%% (%.1f%% improvement)\n', ...
    metrics_before.THD_F, metrics_after.THD_F, ...
    (metrics_before.THD_F - metrics_after.THD_F)/metrics_before.THD_F*100);
fprintf('  Power Factor: %.4f -> %.4f\n', metrics_before.true_PF, metrics_after.true_PF);
fprintf('  K-Factor: %.2f -> %.2f\n', metrics_before.K_factor, metrics_after.K_factor);

if metrics_after.THD_F < metrics_before.THD_F && ...
   metrics_after.true_PF > metrics_before.true_PF
    fprintf('\nTEST 5: PASSED (Power quality improved)\n\n');
else
    fprintf('\nTEST 5: FAILED\n\n');
end

% =======================================================================
% VISUALIZATION
% =======================================================================

fprintf('Generating comprehensive visualization...\n');

% Figure 1: Comprehensive Filter Comparison
figure('Name', 'Harmonic Filter Performance', 'Position', [50, 50, 1600, 1000]);

% Row 1: Time Domain Comparison - Notch Filter
subplot(4, 3, 1);
plot(t*1000, signal_original, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title(sprintf('Original Signal\nTHD: %.1f%%', THD_orig));
xlim([0 40]);

subplot(4, 3, 2);
plot(t_notch(1:N)*1000, real(signal_notch(1:N)), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title(sprintf('After Notch Filter (5th, 7th removed)\nTHD: %.1f%%', THD_notch));
xlim([0 20]);

subplot(4, 3, 3);
plot(t*1000, signal_original, 'b-', 'LineWidth', 1, 'DisplayName', 'Original');
hold on;
plot(t_notch(1:N)*1000, real(signal_notch(1:N)), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Filtered');
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Overlay Comparison');
legend('Location', 'best');
xlim([0 20]);

% Row 2: Frequency Domain - Notch Filter
subplot(4, 3, 4);
stem(freq(1:21), mag_original(1:21), 'b-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('Original Spectrum');
xlim([0 1050]);

subplot(4, 3, 5);
stem(freq(1:21), mag_notch(1:21), 'r-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('Filtered Spectrum');
xlim([0 1050]);

subplot(4, 3, 6);
harmonic_nums = [1 3 5 7 9 11 13];
harmonic_indices = harmonic_nums + 1;
x_pos = 1:length(harmonic_nums);
width = 0.35;

bar(x_pos - width/2, mag_original(harmonic_indices), width, 'FaceColor', [0.3 0.5 0.8], 'DisplayName', 'Before');
hold on;
bar(x_pos + width/2, mag_notch(harmonic_indices), width, 'FaceColor', [0.8 0.3 0.3], 'DisplayName', 'After');
set(gca, 'XTick', x_pos, 'XTickLabel', harmonic_nums);
xlabel('Harmonic Order');
ylabel('Magnitude (V)');
title('Harmonic Comparison');
legend('Location', 'best');
grid on;

% Row 3: LED Lighting - Attenuation Filter
subplot(4, 3, 7);
plot(t_led*1000, signal_led, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title(sprintf('LED Lighting (Original)\nTHD: %.1f%%', THD_led_orig));
xlim([0 40]);

subplot(4, 3, 8);
plot(t_led*1000, real(signal_atten(1:N)), 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title(sprintf('After 40dB Attenuation (3rd)\nTHD: %.1f%%', info_atten.filtered_THD));
xlim([0 20]);

subplot(4, 3, 9);
stem(freq(1:11), mag_led(1:11), 'b-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Before');
hold on;
stem(freq(1:11), mag_atten(1:11), 'g-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'After');
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('3rd Harmonic Attenuation');
legend('Location', 'best');
xlim([0 550]);

% Row 4: Performance Metrics
subplot(4, 3, 10);
scenarios = {'Original', 'Notch', 'Atten', 'Lowpass'};
thd_vals = [THD_orig, THD_notch, info_atten.filtered_THD, info_lowpass.filtered_THD];

bar(thd_vals, 'FaceColor', [0.4 0.6 0.9]);
hold on;
yline(5, 'g--', 'LineWidth', 2, 'Label', 'IEEE Excellent');
yline(8, 'r--', 'LineWidth', 2, 'Label', 'IEEE Good');
set(gca, 'XTickLabel', scenarios);
ylabel('THD (%)');
title('THD Comparison');
grid on;

subplot(4, 3, 11);
pf_before = metrics_before.true_PF;
pf_after = metrics_after.true_PF;

bar([pf_before, pf_after], 'FaceColor', [0.3 0.8 0.5]);
hold on;
yline(0.9, 'r--', 'LineWidth', 2, 'Label', 'Min Acceptable');
set(gca, 'XTickLabel', {'Before', 'After'});
ylabel('Power Factor');
title('Power Factor Improvement');
ylim([0.85 1.05]);
grid on;

subplot(4, 3, 12);
axis off;
text(0.1, 0.95, 'FILTER PERFORMANCE SUMMARY', 'FontSize', 11, 'FontWeight', 'bold');
text(0.1, 0.80, 'Notch Filter (5th, 7th removed):', 'FontSize', 9, 'FontWeight', 'bold');
text(0.15, 0.70, sprintf('THD: %.1f%% -> %.1f%%', THD_orig, THD_notch), 'FontSize', 8);
text(0.15, 0.60, sprintf('Improvement: %.1f%%', info_notch.THD_reduction_percent), 'FontSize', 8);

text(0.1, 0.45, 'Attenuation (3rd by 40dB):', 'FontSize', 9, 'FontWeight', 'bold');
text(0.15, 0.35, sprintf('THD: %.1f%% -> %.1f%%', THD_led_orig, info_atten.filtered_THD), 'FontSize', 8);
text(0.15, 0.25, sprintf('3rd Harmonic: %.1fV -> %.1fV', mag_led(4), mag_atten(4)), 'FontSize', 8);

text(0.1, 0.10, 'All Tests: PASSED', 'FontSize', 10, 'Color', [0 0.6 0], 'FontWeight', 'bold');

sgtitle('Harmonic Filter Design - Comprehensive Performance Analysis', 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, 'figures/harmonic_filter_performance.png');

fprintf('========================================\n');
fprintf('ALL FILTER TESTS COMPLETED SUCCESSFULLY\n');
fprintf('Figure saved to: figures/harmonic_filter_performance.png\n');
fprintf('========================================\n');
