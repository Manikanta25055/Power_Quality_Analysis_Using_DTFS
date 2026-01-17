% SCENARIO_DATA_CENTER - Data Center / SMPS Load Analysis
%
% DESCRIPTION:
%   Analyzes power quality for data centers with switch-mode power supplies.
%   SMPS loads produce multiple odd harmonics with relatively uniform distribution.
%
% TYPICAL CHARACTERISTICS:
%   - Multiple odd harmonics: 3rd (12%), 5th (10%), 7th (8%), 9th (6%), 11th (4%)
%   - Expected THD: 18-22%
%   - Application: Data centers, server farms, IT equipment
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025

clear all; close all; clc;
addpath('../functions');

fprintf('\n========================================================================\n');
fprintf('         DATA CENTER LOAD - POWER QUALITY ANALYSIS\n');
fprintf('========================================================================\n\n');

% System parameters
fs = 10000;
f0 = 50;
N = fs/f0;
V_rms = 230;

% Generate data center scenario
fprintf('Generating data center / SMPS load signal...\n\n');
[signal, t, info] = generate_distorted_signal('data_center', fs, f0, N, V_rms);

% DTFS analysis
fprintf('Performing DTFS analysis...\n');
[X_k, magnitude, phase, frequencies] = calculate_dtfs(signal, fs, N);

% Power quality metrics
fprintf('Calculating power quality metrics...\n');
metrics = power_quality_metrics(signal, magnitude, fs, f0);

% Display results
fprintf('\n========================================================================\n');
fprintf('DATA CENTER / SMPS LOAD - ANALYSIS RESULTS\n');
fprintf('========================================================================\n\n');

fprintf('HARMONIC CONTENT (Uniformly distributed odd harmonics):\n');
fprintf('  Fundamental (50 Hz):   %.2f V (100.0%%)\n', magnitude(2));
fprintf('  3rd Harmonic (150 Hz): %.2f V (%.1f%%)\n', magnitude(4), (magnitude(4)/magnitude(2))*100);
fprintf('  5th Harmonic (250 Hz): %.2f V (%.1f%%)\n', magnitude(6), (magnitude(6)/magnitude(2))*100);
fprintf('  7th Harmonic (350 Hz): %.2f V (%.1f%%)\n', magnitude(8), (magnitude(8)/magnitude(2))*100);
fprintf('  9th Harmonic (450 Hz): %.2f V (%.1f%%)\n', magnitude(10), (magnitude(10)/magnitude(2))*100);
fprintf('  11th Harmonic (550 Hz): %.2f V (%.1f%%)\n', magnitude(12), (magnitude(12)/magnitude(2))*100);
fprintf('  13th Harmonic (650 Hz): %.2f V (%.1f%%)\n', magnitude(14), (magnitude(14)/magnitude(2))*100);
fprintf('\n');

fprintf('POWER QUALITY METRICS:\n');
fprintf('  THD:              %.2f%%\n', metrics.THD_F);
fprintf('  Crest Factor:     %.4f\n', metrics.crest_factor);
fprintf('  True Power Factor: %.4f\n', metrics.true_PF);
fprintf('  K-Factor:         %.2f (Use K-4 or K-13 transformer)\n', metrics.K_factor);
fprintf('  Quality Rating:   %s\n', metrics.quality_rating);
fprintf('  IEEE 519 Status:  %s\n\n', mat2str(metrics.compliance_IEEE519));

fprintf('DATA CENTER CHARACTERISTICS:\n');
fprintf('  - High-frequency switching converters (50-200 kHz)\n');
fprintf('  - Multiple odd harmonics (no dominant harmonic)\n');
fprintf('  - Variable load depending on server utilization\n');
fprintf('  - Power factor correction often included\n\n');

fprintf('RECOMMENDED MITIGATION:\n');
fprintf('  1. Active harmonic filter (adaptive to load changes)\n');
fprintf('  2. K-13 rated transformers for heavy IT loads\n');
fprintf('  3. Hybrid passive-active filtering system\n');
fprintf('  4. Distributed power factor correction\n');
fprintf('\n========================================================================\n');

% Visualization
fig = figure('Name', 'Data Center Analysis', 'Position', [100, 100, 1400, 800]);

% Time domain waveform
subplot(2, 3, 1);
plot(t*1000, signal, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title(sprintf('Data Center Waveform\nTHD: %.1f%%', metrics.THD_F));
xlim([0 40]);

% Frequency spectrum
subplot(2, 3, 2);
stem(frequencies(1:21), magnitude(1:21), 'r', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('Harmonic Spectrum (Uniform Distribution)');
xlim([0 1050]);

% Phase spectrum
subplot(2, 3, 3);
stem(frequencies(1:21), phase(1:21)*180/pi, 'g', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Phase Spectrum');
xlim([0 1050]);

% Harmonic contribution pie chart
subplot(2, 3, 4);
harmonic_labels = {'Fund', '3rd', '5th', '7th', '9th', '11th', '13th', 'Others'};
harmonic_values = [magnitude(2), magnitude(4), magnitude(6), magnitude(8), ...
                   magnitude(10), magnitude(12), magnitude(14), ...
                   sqrt(sum(magnitude(16:end).^2))];
pie(harmonic_values, harmonic_labels);
title('Harmonic Energy Distribution');

% Individual harmonic distortion
subplot(2, 3, 5);
harmonics = [3 5 7 9 11 13];
ihd_values = (magnitude(harmonics+1) ./ magnitude(2)) * 100;
bar(harmonics, ihd_values, 'FaceColor', [0.4 0.6 0.9]);
grid on;
xlabel('Harmonic Order');
ylabel('IHD (%)');
title('Individual Harmonic Distortion');
hold on;
yline(5, 'r--', 'LineWidth', 2, 'Label', 'IEEE Limit (5%)');

% Power quality metrics summary
subplot(2, 3, 6);
axis off;
text(0.1, 0.95, 'POWER QUALITY SUMMARY', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.85, sprintf('THD: %.2f%% (Limit: 8%%)', metrics.THD_F), 'FontSize', 10);
text(0.1, 0.75, sprintf('Crest Factor: %.4f', metrics.crest_factor), 'FontSize', 10);
text(0.1, 0.65, sprintf('Power Factor: %.4f', metrics.true_PF), 'FontSize', 10);
text(0.1, 0.55, sprintf('K-Factor: %.2f', metrics.K_factor), 'FontSize', 10);
text(0.1, 0.45, sprintf('Dominant: 3rd harmonic (%.1f%%)', (magnitude(4)/magnitude(2))*100), 'FontSize', 10);
text(0.1, 0.30, 'CHARACTERISTICS:', 'FontSize', 10, 'FontWeight', 'bold');
text(0.1, 0.20, '- SMPS with PFC (high freq switching)', 'FontSize', 9);
text(0.1, 0.12, '- Multiple odd harmonics (uniform)', 'FontSize', 9);
text(0.1, 0.04, '- Variable load profile', 'FontSize', 9);

sgtitle('Data Center / SMPS Load - Complete Power Quality Analysis', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Analysis complete! Figure window displayed.\n');
