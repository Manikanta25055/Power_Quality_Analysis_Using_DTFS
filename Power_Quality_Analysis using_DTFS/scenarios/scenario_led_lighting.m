% SCENARIO_LED_LIGHTING - LED Lighting Load Analysis
%
% DESCRIPTION:
%   Analyzes power quality for LED lighting systems. LED drivers use
%   single-phase rectifiers that produce strong 3rd harmonic distortion.
%
% TYPICAL CHARACTERISTICS:
%   - 3rd harmonic dominant (150 Hz): 18% of fundamental
%   - Additional odd harmonics: 5th (8%), 7th (4%), 9th (2%)
%   - Expected THD: 18-22%
%   - Application: Residential and commercial lighting
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025

clear all; close all; clc;
addpath('../functions');

fprintf('\n========================================================================\n');
fprintf('         LED LIGHTING LOAD - POWER QUALITY ANALYSIS\n');
fprintf('========================================================================\n\n');

% System parameters
fs = 10000;           % Sampling frequency (Hz)
f0 = 50;              % Fundamental frequency (Hz)
N = fs/f0;            % Period length (200 samples)
V_rms = 230;          % Nominal RMS voltage (V)

% Generate LED lighting scenario
fprintf('Generating LED lighting load signal...\n\n');
[signal, t, info] = generate_distorted_signal('led_lighting', fs, f0, N, V_rms);

% Perform DTFS analysis
fprintf('Performing DTFS analysis...\n');
[X_k, magnitude, phase, frequencies] = calculate_dtfs(signal, fs, N);

% Calculate power quality metrics
fprintf('Calculating power quality metrics...\n');
metrics = power_quality_metrics(signal, magnitude, fs, f0);

% Display results
fprintf('\n========================================================================\n');
fprintf('LED LIGHTING LOAD - ANALYSIS RESULTS\n');
fprintf('========================================================================\n\n');

fprintf('HARMONIC CONTENT:\n');
fprintf('  Fundamental (50 Hz):  %.2f V (100.0%%)\n', magnitude(2));
fprintf('  3rd Harmonic (150 Hz): %.2f V (%.1f%%)\n', magnitude(4), (magnitude(4)/magnitude(2))*100);
fprintf('  5th Harmonic (250 Hz): %.2f V (%.1f%%)\n', magnitude(6), (magnitude(6)/magnitude(2))*100);
fprintf('  7th Harmonic (350 Hz): %.2f V (%.1f%%)\n', magnitude(8), (magnitude(8)/magnitude(2))*100);
fprintf('  9th Harmonic (450 Hz): %.2f V (%.1f%%)\n', magnitude(10), (magnitude(10)/magnitude(2))*100);
fprintf('\n');

fprintf('POWER QUALITY METRICS:\n');
fprintf('  THD:              %.2f%%\n', metrics.THD_F);
fprintf('  Crest Factor:     %.4f\n', metrics.crest_factor);
fprintf('  True Power Factor: %.4f\n', metrics.true_PF);
fprintf('  K-Factor:         %.2f\n', metrics.K_factor);
fprintf('  Quality Rating:   %s\n', metrics.quality_rating);
fprintf('  IEEE 519 Status:  %s\n\n', mat2str(metrics.compliance_IEEE519));

fprintf('RECOMMENDED ACTIONS:\n');
fprintf('  %s\n', metrics.recommended_action);
fprintf('\n========================================================================\n');

% Visualization
fig = figure('Name', 'LED Lighting Analysis', 'Position', [100, 100, 1400, 800]);

% Time domain waveform
subplot(2, 3, 1);
plot(t*1000, signal, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Voltage (V)');
title(sprintf('LED Lighting Waveform\nTHD: %.1f%%', metrics.THD_F));
xlim([0 40]);

% Frequency spectrum (magnitude)
subplot(2, 3, 2);
stem(frequencies(1:21), magnitude(1:21), 'r', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('Harmonic Spectrum');
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
harmonic_labels = {'Fund (50Hz)', '3rd (150Hz)', '5th (250Hz)', '7th (350Hz)', '9th (450Hz)', 'Others'};
harmonic_values = [magnitude(2), magnitude(4), magnitude(6), magnitude(8), magnitude(10), ...
                   sqrt(sum(magnitude(12:end).^2))];
pie(harmonic_values, harmonic_labels);
title('Harmonic Energy Distribution');

% THD comparison with IEEE 519 limits
subplot(2, 3, 5);
bar([metrics.THD_F, 5, 8], 'FaceColor', [0.8 0.4 0.4]);
hold on;
yline(5, 'g--', 'LineWidth', 2, 'Label', 'Excellent (<5%)');
yline(8, 'r--', 'LineWidth', 2, 'Label', 'Acceptable (<8%)');
set(gca, 'XTickLabel', {'Measured', 'IEEE Excellent', 'IEEE Limit'});
ylabel('THD (%)');
title('THD vs IEEE 519 Standards');
grid on;

% Power quality metrics summary
subplot(2, 3, 6);
axis off;
text(0.1, 0.95, 'POWER QUALITY SUMMARY', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.85, sprintf('THD: %.2f%% (IEEE Limit: 8%%)', metrics.THD_F), 'FontSize', 10);
text(0.1, 0.75, sprintf('Crest Factor: %.4f (Ideal: 1.414)', metrics.crest_factor), 'FontSize', 10);
text(0.1, 0.65, sprintf('Power Factor: %.4f', metrics.true_PF), 'FontSize', 10);
text(0.1, 0.55, sprintf('K-Factor: %.2f (Use K-4 transformer)', metrics.K_factor), 'FontSize', 10);
text(0.1, 0.45, sprintf('Dominant: 3rd harmonic (%.1f%%)', (magnitude(4)/magnitude(2))*100), 'FontSize', 10);
text(0.1, 0.30, 'CHARACTERISTICS:', 'FontSize', 10, 'FontWeight', 'bold');
text(0.1, 0.20, '- Single-phase rectifiers in LED drivers', 'FontSize', 9);
text(0.1, 0.12, '- Strong triplen harmonics (3rd, 9th)', 'FontSize', 9);
text(0.1, 0.04, '- Mitigation: Delta-wye transformer', 'FontSize', 9);

sgtitle('LED Lighting Load - Complete Power Quality Analysis', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Analysis complete! Figure window displayed.\n');
