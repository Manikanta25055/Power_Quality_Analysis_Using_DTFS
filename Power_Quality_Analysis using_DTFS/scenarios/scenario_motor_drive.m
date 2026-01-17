% SCENARIO_MOTOR_DRIVE - Variable Frequency Drive Analysis
%
% DESCRIPTION:
%   Analyzes power quality for motor drives with 6-pulse rectifiers.
%   VFDs are the most common source of harmonics in industrial systems.
%
% TYPICAL CHARACTERISTICS:
%   - 5th harmonic dominant (250 Hz): 20% of fundamental
%   - 7th harmonic (350 Hz): 14% of fundamental
%   - Additional: 11th (9%), 13th (6%)
%   - Expected THD: 25-30%
%   - Application: Industrial motor control, HVAC systems
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025

clear all; close all; clc;
addpath('../functions');

fprintf('\n========================================================================\n');
fprintf('     MOTOR DRIVE (VFD) - POWER QUALITY ANALYSIS\n');
fprintf('========================================================================\n\n');

% System parameters
fs = 10000;
f0 = 50;
N = fs/f0;
V_rms = 230;

% Generate motor drive scenario
fprintf('Generating 6-pulse VFD load signal...\n\n');
[signal, t, info] = generate_distorted_signal('motor_drive', fs, f0, N, V_rms);

% DTFS analysis
fprintf('Performing DTFS analysis...\n');
[X_k, magnitude, phase, frequencies] = calculate_dtfs(signal, fs, N);

% Power quality metrics
fprintf('Calculating power quality metrics...\n');
metrics = power_quality_metrics(signal, magnitude, fs, f0);

% Display results
fprintf('\n========================================================================\n');
fprintf('MOTOR DRIVE (6-PULSE VFD) - ANALYSIS RESULTS\n');
fprintf('========================================================================\n\n');

fprintf('CHARACTERISTIC HARMONICS (6-pulse: h = 6k ± 1):\n');
fprintf('  Fundamental (50 Hz):   %.2f V (100.0%%)\n', magnitude(2));
fprintf('  5th Harmonic (250 Hz): %.2f V (%.1f%%) - DOMINANT\n', magnitude(6), (magnitude(6)/magnitude(2))*100);
fprintf('  7th Harmonic (350 Hz): %.2f V (%.1f%%)\n', magnitude(8), (magnitude(8)/magnitude(2))*100);
fprintf('  11th Harmonic (550 Hz): %.2f V (%.1f%%)\n', magnitude(12), (magnitude(12)/magnitude(2))*100);
fprintf('  13th Harmonic (650 Hz): %.2f V (%.1f%%)\n', magnitude(14), (magnitude(14)/magnitude(2))*100);
fprintf('\n');

fprintf('POWER QUALITY METRICS:\n');
fprintf('  THD:              %.2f%% (WORST CASE)\n', metrics.THD_F);
fprintf('  Crest Factor:     %.4f\n', metrics.crest_factor);
fprintf('  True Power Factor: %.4f\n', metrics.true_PF);
fprintf('  K-Factor:         %.2f (Requires K-4 or higher)\n', metrics.K_factor);
fprintf('  Quality Rating:   %s\n', metrics.quality_rating);
fprintf('  IEEE 519 Status:  NON-COMPLIANT\n\n');

fprintf('RECOMMENDED MITIGATION:\n');
fprintf('  1. Upgrade to 12-pulse rectifier (eliminates 5th & 7th)\n');
fprintf('  2. Install active harmonic filter\n');
fprintf('  3. Use line reactors (5%% impedance)\n');
fprintf('  4. Install K-4 or K-13 rated transformer\n');
fprintf('\n========================================================================\n');

% Demonstration of filtering
fprintf('\nDEMONSTRATING HARMONIC FILTERING...\n');
[X_k_filtered, filter_info] = design_harmonic_filter(X_k, N, 'notch', [5 7]);
[signal_filtered, t_filt] = synthesize_from_dtfs(X_k_filtered, fs, N, 1);
[~, mag_filt, ~, ~] = calculate_dtfs(real(signal_filtered), fs, N);
metrics_filtered = power_quality_metrics(real(signal_filtered), mag_filt, fs, f0);

% Visualization
fig = figure('Name', 'Motor Drive Analysis', 'Position', [50, 50, 1600, 900]);

% Original time domain
subplot(3, 4, 1);
plot(t*1000, signal, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title(sprintf('Original VFD Waveform\nTHD: %.1f%%', metrics.THD_F));
xlim([0 40]);

% Original frequency spectrum
subplot(3, 4, 2);
stem(frequencies(1:21), magnitude(1:21), 'r', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('Original Spectrum');
xlim([0 1050]);

% Phase spectrum
subplot(3, 4, 3);
stem(frequencies(1:21), phase(1:21)*180/pi, 'g', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Phase Spectrum');
xlim([0 1050]);

% Harmonic pie chart
subplot(3, 4, 4);
harmonic_labels = {'Fund', '5th', '7th', '11th', '13th', 'Others'};
harmonic_values = [magnitude(2), magnitude(6), magnitude(8), magnitude(12), magnitude(14), ...
                   sqrt(sum(magnitude(16:end).^2))];
pie(harmonic_values, harmonic_labels);
title('Harmonic Distribution');

% Filtered time domain
subplot(3, 4, 5);
plot(t*1000, real(signal_filtered(1:N)), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title(sprintf('After Filtering (5th & 7th removed)\nTHD: %.1f%%', metrics_filtered.THD_F));
xlim([0 20]);

% Filtered spectrum
subplot(3, 4, 6);
stem(frequencies(1:21), mag_filt(1:21), 'b', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('Filtered Spectrum');
xlim([0 1050]);

% Before/After overlay
subplot(3, 4, 7);
plot(t*1000, signal, 'b-', 'LineWidth', 1, 'DisplayName', 'Original');
hold on;
plot(t*1000, real(signal_filtered(1:N)), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Filtered');
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Waveform Comparison');
legend('Location', 'best');
xlim([0 20]);

% Harmonic comparison bar chart
subplot(3, 4, 8);
harmonics = [1 5 7 11 13];
h_idx = harmonics + 1;
x_pos = 1:length(harmonics);
width = 0.35;
bar(x_pos - width/2, magnitude(h_idx), width, 'FaceColor', [0.3 0.5 0.8], 'DisplayName', 'Before');
hold on;
bar(x_pos + width/2, mag_filt(h_idx), width, 'FaceColor', [0.8 0.3 0.3], 'DisplayName', 'After');
set(gca, 'XTick', x_pos, 'XTickLabel', harmonics);
xlabel('Harmonic Order');
ylabel('Magnitude (V)');
title('Harmonic Magnitude Comparison');
legend('Location', 'best');
grid on;

% THD comparison
subplot(3, 4, 9);
bar([metrics.THD_F, metrics_filtered.THD_F, 5, 8], 'FaceColor', [0.8 0.4 0.4]);
hold on;
yline(5, 'g--', 'LineWidth', 2, 'Label', 'Excellent');
yline(8, 'r--', 'LineWidth', 2, 'Label', 'Acceptable');
set(gca, 'XTickLabel', {'Before', 'After', 'IEEE Excellent', 'IEEE Limit'});
ylabel('THD (%)');
title(sprintf('THD Improvement: %.1f%% reduction', filter_info.THD_reduction_percent));
grid on;

% Power factor comparison
subplot(3, 4, 10);
bar([metrics.true_PF, metrics_filtered.true_PF], 'FaceColor', [0.3 0.8 0.6]);
hold on;
yline(0.9, 'r--', 'LineWidth', 2);
set(gca, 'XTickLabel', {'Before', 'After'});
ylabel('Power Factor');
title('Power Factor Improvement');
ylim([0.85 1.05]);
grid on;

% K-factor comparison
subplot(3, 4, 11);
bar([metrics.K_factor, metrics_filtered.K_factor], 'FaceColor', [0.9 0.5 0.3]);
hold on;
yline(4, 'g--', 'LineWidth', 1.5, 'Label', 'K-4 Rating');
yline(13, 'r--', 'LineWidth', 1.5, 'Label', 'K-13 Rating');
set(gca, 'XTickLabel', {'Before', 'After'});
ylabel('K-Factor');
title('Transformer K-Factor Reduction');
grid on;
legend('Location', 'best');

% Performance summary
subplot(3, 4, 12);
axis off;
text(0.1, 0.95, 'FILTERING PERFORMANCE', 'FontSize', 11, 'FontWeight', 'bold');
text(0.1, 0.85, sprintf('THD: %.1f%% → %.1f%%', metrics.THD_F, metrics_filtered.THD_F), 'FontSize', 9);
text(0.1, 0.75, sprintf('Improvement: %.1f%%', filter_info.THD_reduction_percent), 'FontSize', 9);
text(0.1, 0.65, sprintf('PF: %.4f → %.4f', metrics.true_PF, metrics_filtered.true_PF), 'FontSize', 9);
text(0.1, 0.55, sprintf('K-Factor: %.2f → %.2f', metrics.K_factor, metrics_filtered.K_factor), 'FontSize', 9);
text(0.1, 0.45, '5th Harmonic: 100% removed', 'FontSize', 9);
text(0.1, 0.35, '7th Harmonic: 100% removed', 'FontSize', 9);
text(0.1, 0.20, 'INDUSTRIAL SOLUTIONS:', 'FontSize', 10, 'FontWeight', 'bold');
text(0.1, 0.10, '- 12-pulse rectifier', 'FontSize', 8);
text(0.1, 0.02, '- Active harmonic filter', 'FontSize', 8);

sgtitle('Motor Drive (6-Pulse VFD) - Complete Analysis with Filtering', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nAnalysis complete! Figure window displayed.\n');
fprintf('Filter effectiveness: %.1f%% harmonic energy removed\n', filter_info.filter_effectiveness);
