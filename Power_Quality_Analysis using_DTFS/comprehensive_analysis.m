% COMPREHENSIVE_ANALYSIS - Complete Power Quality Analysis with All Visualizations
%
% DESCRIPTION:
%   This script provides the most comprehensive power quality analysis with:
%   1. Time domain comparison
%   2. DTFS magnitude spectrum
%   3. Phase spectrum (ADDED)
%   4. Harmonic contribution pie charts (ADDED)
%   5. THD analysis charts
%   6. Power quality dashboard
%   7. Filtering results with before/after comparison
%   8. Detailed metrics tables
%
% INCLUDES ALL REQUIREMENTS FROM PROJECT SPECIFICATION
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025
% COURSE: Digital Signal Processing (FISAC Assessment)

clear all;
close all;
clc;

% Add functions folder to path
addpath('functions');

fprintf('\n');
fprintf('================================================================================\n');
fprintf('     COMPREHENSIVE POWER QUALITY ANALYSIS USING DTFS\n');
fprintf('     Complete Analysis with All Visualizations\n');
fprintf('================================================================================\n');
fprintf('Project:  Power Quality Assessment and Harmonic Mitigation\n');
fprintf('Author:   Manikanta Gonugondla (230906450)\n');
fprintf('Course:   Digital Signal Processing (FISAC Assessment)\n');
fprintf('Date:     October 2025\n');
fprintf('================================================================================\n\n');

%% SYSTEM PARAMETERS
fprintf('Initializing system parameters...\n\n');

fs = 10000;           % Sampling frequency (Hz)
f0 = 50;              % Fundamental frequency (Hz) - Indian grid
N = fs/f0;            % Samples per period (200 samples)
V_rms_nominal = 230;  % Nominal RMS voltage (V)

fprintf('SYSTEM CONFIGURATION:\n');
fprintf('  Sampling Frequency:    %d Hz\n', fs);
fprintf('  Fundamental Frequency: %d Hz\n', f0);
fprintf('  Period Length:         %d samples\n', N);
fprintf('  Nominal Voltage:       %d V RMS\n\n', V_rms_nominal);

%% PART 1: GENERATE AND ANALYZE MULTIPLE SCENARIOS
fprintf('================================================================================\n');
fprintf('PART 1: MULTI-SCENARIO ANALYSIS\n');
fprintf('================================================================================\n\n');

scenarios = {'ideal', 'led_lighting', 'motor_drive', 'data_center'};
scenario_names = {'Ideal Power', 'LED Lighting', 'Motor Drive (VFD)', 'Data Center'};
num_scenarios = length(scenarios);

% Storage
signals = cell(num_scenarios, 1);
dtfs_results = cell(num_scenarios, 1);
metrics_results = cell(num_scenarios, 1);

% Analyze each scenario
for i = 1:num_scenarios
    fprintf('--- Analyzing Scenario %d/%d: %s ---\n', i, num_scenarios, scenario_names{i});

    [sig, t, info] = generate_distorted_signal(scenarios{i}, fs, f0, N, V_rms_nominal);
    signals{i}.waveform = sig;
    signals{i}.time = t;
    signals{i}.info = info;

    [X_k, mag, phase, freq] = calculate_dtfs(sig, fs, N);
    dtfs_results{i}.X_k = X_k;
    dtfs_results{i}.magnitude = mag;
    dtfs_results{i}.phase = phase;
    dtfs_results{i}.frequencies = freq;

    metrics = power_quality_metrics(sig, mag, fs, f0);
    metrics_results{i} = metrics;

    fprintf('\n');
end

%% PART 2: COMPARATIVE METRICS TABLE
fprintf('================================================================================\n');
fprintf('PART 2: COMPARATIVE METRICS TABLE\n');
fprintf('================================================================================\n\n');

fprintf('┌──────────────────────┬─────────┬─────────┬────────┬──────────┬─────────────┐\n');
fprintf('│ Scenario             │ THD (%%) │ Cr.Fact │   PF   │ K-Factor │ IEEE Status │\n');
fprintf('├──────────────────────┼─────────┼─────────┼────────┼──────────┼─────────────┤\n');

for i = 1:num_scenarios
    m = metrics_results{i};
    status_str = '✓ PASS';
    if ~m.compliance_IEEE519
        status_str = '✗ FAIL';
    end
    fprintf('│ %-20s │ %7.2f │ %7.4f │ %6.4f │ %8.2f │ %-11s │\n', ...
        scenario_names{i}, m.THD_F, m.crest_factor, m.true_PF, m.K_factor, status_str);
end

fprintf('└──────────────────────┴─────────┴─────────┴────────┴──────────┴─────────────┘\n\n');

%% VISUALIZATION 1: TIME DOMAIN COMPARISON
fprintf('Creating Figure 1: Time Domain Comparison...\n');

fig1 = figure('Name', 'Time Domain Comparison', 'Position', [50, 50, 1600, 900]);

for i = 1:num_scenarios
    subplot(2, 2, i);
    plot(signals{i}.time*1000, signals{i}.waveform, 'LineWidth', 1.5);
    grid on;
    xlabel('Time (ms)', 'FontSize', 11);
    ylabel('Voltage (V)', 'FontSize', 11);
    title(sprintf('%s\nTHD: %.1f%%, CF: %.3f, PF: %.3f', ...
        scenario_names{i}, metrics_results{i}.THD_F, ...
        metrics_results{i}.crest_factor, metrics_results{i}.true_PF), ...
        'FontSize', 12);
    xlim([0 40]);
end

sgtitle('Power Quality Scenarios - Time Domain Waveforms', 'FontSize', 14, 'FontWeight', 'bold');

%% VISUALIZATION 2: MAGNITUDE SPECTRUM
fprintf('Creating Figure 2: DTFS Magnitude Spectrum...\n');

fig2 = figure('Name', 'Magnitude Spectrum', 'Position', [100, 100, 1600, 900]);

for i = 1:num_scenarios
    subplot(2, 2, i);
    mag = dtfs_results{i}.magnitude;
    freq = dtfs_results{i}.frequencies;

    max_idx = find(freq <= 1000, 1, 'last');
    stem(freq(1:max_idx), mag(1:max_idx), 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    xlabel('Frequency (Hz)', 'FontSize', 11);
    ylabel('Magnitude (V)', 'FontSize', 11);
    title(sprintf('%s - Harmonic Spectrum\nDominant: %dth harmonic (%.0f Hz)', ...
        scenario_names{i}, metrics_results{i}.dominant_harmonic, ...
        metrics_results{i}.dominant_harmonic*f0), 'FontSize', 12);
    xlim([0 1000]);
end

sgtitle('Power Quality Scenarios - DTFS Magnitude Spectrum', 'FontSize', 14, 'FontWeight', 'bold');

%% VISUALIZATION 3: PHASE SPECTRUM (NEW)
fprintf('Creating Figure 3: Phase Spectrum...\n');

fig3 = figure('Name', 'Phase Spectrum', 'Position', [150, 150, 1600, 900]);

for i = 1:num_scenarios
    subplot(2, 2, i);
    phase = dtfs_results{i}.phase;
    freq = dtfs_results{i}.frequencies;

    max_idx = find(freq <= 1000, 1, 'last');
    stem(freq(1:max_idx), phase(1:max_idx)*180/pi, 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'Color', [0.2 0.7 0.3]);
    grid on;
    xlabel('Frequency (Hz)', 'FontSize', 11);
    ylabel('Phase (degrees)', 'FontSize', 11);
    title(sprintf('%s - Phase Spectrum', scenario_names{i}), 'FontSize', 12);
    xlim([0 1000]);
    ylim([-200 200]);
end

sgtitle('Power Quality Scenarios - Phase Spectrum Analysis', 'FontSize', 14, 'FontWeight', 'bold');

%% VISUALIZATION 4: HARMONIC CONTRIBUTION PIE CHARTS (NEW)
fprintf('Creating Figure 4: Harmonic Contribution Pie Charts...\n');

fig4 = figure('Name', 'Harmonic Contribution', 'Position', [200, 200, 1600, 900]);

for i = 1:num_scenarios
    subplot(2, 2, i);
    mag = dtfs_results{i}.magnitude;

    % Define harmonic contributions
    if i == 1  % Ideal
        labels = {'Fundamental'};
        values = [mag(2)];
    elseif i == 2  % LED
        labels = {'Fund (50Hz)', '3rd (150Hz)', '5th (250Hz)', '7th (350Hz)', '9th (450Hz)', 'Others'};
        values = [mag(2), mag(4), mag(6), mag(8), mag(10), sqrt(sum(mag(12:end).^2))];
    elseif i == 3  % Motor
        labels = {'Fund (50Hz)', '5th (250Hz)', '7th (350Hz)', '11th (550Hz)', '13th (650Hz)', 'Others'};
        values = [mag(2), mag(6), mag(8), mag(12), mag(14), sqrt(sum(mag([4,5,7,9,10,11,13,15:end]).^2))];
    else  % Data Center
        labels = {'Fund', '3rd', '5th', '7th', '9th', '11th', '13th', 'Others'};
        values = [mag(2), mag(4), mag(6), mag(8), mag(10), mag(12), mag(14), sqrt(sum(mag(16:end).^2))];
    end

    % Create pie chart with percentage labels
    pie(values.^2, labels);  % Use energy (squared) for more visible differences
    title(sprintf('%s\nHarmonic Energy Distribution', scenario_names{i}), ...
        'FontSize', 12);
end

sgtitle('Harmonic Contribution Analysis (Energy-based)', 'FontSize', 14, 'FontWeight', 'bold');

%% VISUALIZATION 5: THD COMPARISON WITH IEEE 519 LIMITS
fprintf('Creating Figure 5: THD Analysis Chart...\n');

fig5 = figure('Name', 'THD Analysis', 'Position', [250, 250, 1400, 700]);

% THD bar chart
subplot(1, 2, 1);
thd_values = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    thd_values(i) = metrics_results{i}.THD_F;
end

colors = [0.4 0.8 0.4;  % Green for ideal
          0.9 0.6 0.3;  % Orange for LED
          0.9 0.3 0.3;  % Red for motor
          0.7 0.5 0.8]; % Purple for data center

for i = 1:num_scenarios
    bar(i, thd_values(i), 'FaceColor', colors(i, :));
    hold on;
end

yline(5, 'g--', 'LineWidth', 3, 'Label', 'Excellent (<5%)', 'FontSize', 11);
yline(8, 'r--', 'LineWidth', 3, 'Label', 'IEEE 519 Limit (<8%)', 'FontSize', 11);
set(gca, 'XTickLabel', scenario_names, 'FontSize', 11);
ylabel('THD (%)', 'FontSize', 12);
title('Total Harmonic Distortion Comparison', 'FontSize', 13);
grid on;
ylim([0 max(thd_values)*1.2]);

% Individual scenarios with harmonic breakdown
subplot(1, 2, 2);
harmonics = [3, 5, 7, 9, 11, 13];
scenario_idx_for_detail = [2, 3, 4];  % LED, Motor, Data Center

for s = 1:length(scenario_idx_for_detail)
    idx = scenario_idx_for_detail(s);
    mag = dtfs_results{idx}.magnitude;
    ihd_values = (mag(harmonics+1) ./ mag(2)) * 100;

    x_offset = (s-1)*0.25;
    bar(harmonics + x_offset, ihd_values, 0.2, 'DisplayName', scenario_names{idx});
    hold on;
end

xlabel('Harmonic Order', 'FontSize', 12);
ylabel('Individual Harmonic Distortion (%)', 'FontSize', 12);
title('Individual Harmonic Distortion (IHD) Breakdown', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 10);
grid on;

sgtitle('THD Analysis with IEEE 519-2014 Standards', 'FontSize', 14, 'FontWeight', 'bold');

%% PART 3: HARMONIC FILTERING DEMONSTRATION
fprintf('\n');
fprintf('================================================================================\n');
fprintf('PART 3: HARMONIC FILTERING (Motor Drive Scenario)\n');
fprintf('================================================================================\n\n');

motor_idx = 3;
sig_motor_orig = signals{motor_idx}.waveform;
X_k_motor_orig = dtfs_results{motor_idx}.X_k;
mag_motor_orig = dtfs_results{motor_idx}.magnitude;
metrics_motor_orig = metrics_results{motor_idx};

fprintf('BEFORE FILTERING:\n');
fprintf('  THD:          %.2f%%\n', metrics_motor_orig.THD_F);
fprintf('  Crest Factor: %.4f\n', metrics_motor_orig.crest_factor);
fprintf('  Power Factor: %.4f\n', metrics_motor_orig.true_PF);
fprintf('  K-Factor:     %.2f\n\n', metrics_motor_orig.K_factor);

% Apply notch filter
fprintf('Applying notch filter to remove 5th and 7th harmonics...\n');
[X_k_motor_filt, filter_info] = design_harmonic_filter(X_k_motor_orig, N, 'notch', [5 7]);
[sig_motor_filt, t_filt] = synthesize_from_dtfs(X_k_motor_filt, fs, N, 1);
[~, mag_motor_filt, ~, ~] = calculate_dtfs(real(sig_motor_filt), fs, N);
metrics_motor_filt = power_quality_metrics(real(sig_motor_filt), mag_motor_filt, fs, f0);

fprintf('AFTER FILTERING:\n');
fprintf('  THD:          %.2f%%\n', metrics_motor_filt.THD_F);
fprintf('  Crest Factor: %.4f\n', metrics_motor_filt.crest_factor);
fprintf('  Power Factor: %.4f\n', metrics_motor_filt.true_PF);
fprintf('  K-Factor:     %.2f\n\n', metrics_motor_filt.K_factor);

%% VISUALIZATION 6: FILTERING RESULTS
fprintf('Creating Figure 6: Filtering Results...\n');

fig6 = figure('Name', 'Filtering Results', 'Position', [300, 300, 1600, 900]);

% Before time domain
subplot(2, 4, 1);
plot(signals{motor_idx}.time*1000, sig_motor_orig, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)', 'FontSize', 10);
ylabel('Voltage (V)', 'FontSize', 10);
title(sprintf('Before Filtering\nTHD: %.1f%%', metrics_motor_orig.THD_F), 'FontSize', 11);
xlim([0 40]);

% After time domain
subplot(2, 4, 2);
plot(t_filt*1000, real(sig_motor_filt), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)', 'FontSize', 10);
ylabel('Voltage (V)', 'FontSize', 10);
title(sprintf('After Filtering\nTHD: %.1f%%', metrics_motor_filt.THD_F), 'FontSize', 11);
xlim([0 20]);

% Overlay
subplot(2, 4, 3);
plot(signals{motor_idx}.time*1000, sig_motor_orig, 'b-', 'LineWidth', 1, 'DisplayName', 'Original');
hold on;
plot(t_filt(1:N)*1000, real(sig_motor_filt(1:N)), 'r-', 'LineWidth', 1.8, 'DisplayName', 'Filtered');
grid on;
xlabel('Time (ms)', 'FontSize', 10);
ylabel('Voltage (V)', 'FontSize', 10);
title('Waveform Overlay', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 9);
xlim([0 20]);

% Improvement metrics
subplot(2, 4, 4);
axis off;
text(0.1, 0.95, 'FILTERING IMPROVEMENTS', 'FontSize', 11, 'FontWeight', 'bold');
text(0.1, 0.80, sprintf('THD: %.1f%% → %.1f%%', metrics_motor_orig.THD_F, metrics_motor_filt.THD_F), 'FontSize', 10);
text(0.1, 0.65, sprintf('Reduction: %.1f%%', filter_info.THD_reduction_percent), 'FontSize', 10);
text(0.1, 0.50, sprintf('PF: %.4f → %.4f', metrics_motor_orig.true_PF, metrics_motor_filt.true_PF), 'FontSize', 10);
text(0.1, 0.35, sprintf('K-Factor: %.2f → %.2f', metrics_motor_orig.K_factor, metrics_motor_filt.K_factor), 'FontSize', 10);
text(0.1, 0.20, sprintf('CF: %.4f → %.4f', metrics_motor_orig.crest_factor, metrics_motor_filt.crest_factor), 'FontSize', 10);
text(0.1, 0.05, sprintf('Energy Removed: %.1f%%', filter_info.filter_effectiveness), 'FontSize', 10);

% Before spectrum
subplot(2, 4, 5);
stem(dtfs_results{motor_idx}.frequencies(1:21), mag_motor_orig(1:21), 'b', 'LineWidth', 1.5, 'MarkerSize', 7);
grid on;
xlabel('Frequency (Hz)', 'FontSize', 10);
ylabel('Magnitude (V)', 'FontSize', 10);
title('Original Spectrum', 'FontSize', 11);
xlim([0 1050]);

% After spectrum
subplot(2, 4, 6);
stem(dtfs_results{motor_idx}.frequencies(1:21), mag_motor_filt(1:21), 'r', 'LineWidth', 1.5, 'MarkerSize', 7);
grid on;
xlabel('Frequency (Hz)', 'FontSize', 10);
ylabel('Magnitude (V)', 'FontSize', 10);
title('Filtered Spectrum (5th & 7th removed)', 'FontSize', 11);
xlim([0 1050]);

% Harmonic comparison bars
subplot(2, 4, 7);
harmonics = [1 5 7 11 13];
h_idx = harmonics + 1;
x_pos = 1:length(harmonics);
width = 0.35;

bar(x_pos - width/2, mag_motor_orig(h_idx), width, 'FaceColor', [0.3 0.5 0.8], 'DisplayName', 'Before');
hold on;
bar(x_pos + width/2, mag_motor_filt(h_idx), width, 'FaceColor', [0.8 0.3 0.3], 'DisplayName', 'After');
set(gca, 'XTick', x_pos, 'XTickLabel', harmonics, 'FontSize', 10);
xlabel('Harmonic Order', 'FontSize', 10);
ylabel('Magnitude (V)', 'FontSize', 10);
title('Harmonic Magnitude Comparison', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 9);
grid on;

% Pie chart comparison
subplot(2, 4, 8);
pie_labels = {'Fund', '5th', '7th', 'Others'};
pie_before = [mag_motor_orig(2)^2, mag_motor_orig(6)^2, mag_motor_orig(8)^2, ...
              sum(mag_motor_orig([4,5,7,9,10,11,12,13,14,15:end]).^2)];
pie(pie_before, pie_labels);
title('Energy Distribution (Before)', 'FontSize', 11);

sgtitle('Motor Drive Filtering Results - Complete Analysis', 'FontSize', 14, 'FontWeight', 'bold');

%% VISUALIZATION 7: POWER QUALITY DASHBOARD
fprintf('Creating Figure 7: Power Quality Dashboard...\n');

fig7 = figure('Name', 'Power Quality Dashboard', 'Position', [350, 350, 1600, 900]);

% THD comparison
subplot(2, 3, 1);
thd_vals = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    thd_vals(i) = metrics_results{i}.THD_F;
end
bar(thd_vals, 'FaceColor', [0.4 0.6 0.9]);
hold on;
yline(5, 'g--', 'LineWidth', 2, 'Label', 'Excellent');
yline(8, 'r--', 'LineWidth', 2, 'Label', 'Limit');
set(gca, 'XTickLabel', scenario_names, 'FontSize', 9);
ylabel('THD (%)', 'FontSize', 11);
title('Total Harmonic Distortion', 'FontSize', 12);
grid on;

% Power Factor
subplot(2, 3, 2);
pf_vals = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    pf_vals(i) = metrics_results{i}.true_PF;
end
bar(pf_vals, 'FaceColor', [0.3 0.9 0.6]);
hold on;
yline(0.9, 'r--', 'LineWidth', 2);
set(gca, 'XTickLabel', scenario_names, 'FontSize', 9);
ylabel('Power Factor', 'FontSize', 11);
title('True Power Factor', 'FontSize', 12);
ylim([0.85 1.05]);
grid on;

% K-Factor
subplot(2, 3, 3);
k_vals = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    k_vals(i) = metrics_results{i}.K_factor;
end
bar(k_vals, 'FaceColor', [0.9 0.3 0.6]);
hold on;
yline(4, 'y--', 'LineWidth', 1.5, 'Label', 'K-4');
yline(13, 'r--', 'LineWidth', 1.5, 'Label', 'K-13');
set(gca, 'XTickLabel', scenario_names, 'FontSize', 9);
ylabel('K-Factor', 'FontSize', 11);
title('Transformer K-Factor', 'FontSize', 12);
grid on;
legend('Location', 'best', 'FontSize', 9);

% Crest Factor
subplot(2, 3, 4);
cf_vals = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    cf_vals(i) = metrics_results{i}.crest_factor;
end
bar(cf_vals, 'FaceColor', [0.9 0.6 0.3]);
hold on;
yline(1.414, 'g--', 'LineWidth', 2, 'Label', 'Ideal');
set(gca, 'XTickLabel', scenario_names, 'FontSize', 9);
ylabel('Crest Factor', 'FontSize', 11);
title('Crest Factor (Peak Stress)', 'FontSize', 12);
grid on;
legend('Location', 'best', 'FontSize', 9);

% IEEE 519 Compliance
subplot(2, 3, 5);
compliance = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    compliance(i) = metrics_results{i}.compliance_IEEE519;
end
bar(compliance, 'FaceColor', [0.5 0.8 0.5]);
set(gca, 'XTickLabel', scenario_names, 'YTick', [0 1], 'YTickLabel', {'Fail', 'Pass'}, 'FontSize', 9);
ylabel('Compliance Status', 'FontSize', 11);
title('IEEE 519 Compliance (THD < 5%)', 'FontSize', 12);
ylim([0 1.2]);
grid on;

% Summary
subplot(2, 3, 6);
axis off;
text(0.1, 0.95, 'PROJECT SUMMARY', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.82, sprintf('Scenarios: %d', num_scenarios), 'FontSize', 10);
text(0.1, 0.70, sprintf('Worst THD: %.1f%% (%s)', max(thd_vals), scenario_names{thd_vals == max(thd_vals)}), 'FontSize', 10);
text(0.1, 0.58, sprintf('Best PF: %.4f (%s)', max(pf_vals), scenario_names{pf_vals == max(pf_vals)}), 'FontSize', 10);
text(0.1, 0.46, sprintf('Highest K: %.2f (%s)', max(k_vals), scenario_names{k_vals == max(k_vals)}), 'FontSize', 10);
text(0.1, 0.30, 'FILTERING PERFORMANCE:', 'FontSize', 11, 'FontWeight', 'bold');
text(0.1, 0.18, sprintf('THD: %.1f%% → %.1f%%', metrics_motor_orig.THD_F, metrics_motor_filt.THD_F), 'FontSize', 10);
text(0.1, 0.06, sprintf('Improvement: %.1f%%', filter_info.THD_reduction_percent), 'FontSize', 10);

sgtitle('Comprehensive Power Quality Metrics Dashboard', 'FontSize', 14, 'FontWeight', 'bold');

%% FINAL SUMMARY
fprintf('\n');
fprintf('================================================================================\n');
fprintf('COMPREHENSIVE ANALYSIS COMPLETE\n');
fprintf('================================================================================\n\n');

fprintf('FIGURES GENERATED:\n');
fprintf('  Figure 1: Time Domain Comparison\n');
fprintf('  Figure 2: DTFS Magnitude Spectrum\n');
fprintf('  Figure 3: Phase Spectrum Analysis\n');
fprintf('  Figure 4: Harmonic Contribution Pie Charts\n');
fprintf('  Figure 5: THD Analysis Chart\n');
fprintf('  Figure 6: Filtering Results\n');
fprintf('  Figure 7: Power Quality Dashboard\n\n');

fprintf('KEY FINDINGS:\n');
fprintf('  1. Motor drives produce worst THD (%.1f%%)\n', max(thd_vals));
fprintf('  2. All non-linear loads violate IEEE 519 (>5%%)\n');
fprintf('  3. Filtering achieved %.1f%% THD reduction\n', filter_info.THD_reduction_percent);
fprintf('  4. Power factor improved from %.4f to %.4f\n', metrics_motor_orig.true_PF, metrics_motor_filt.true_PF);
fprintf('  5. K-factor reduced from %.2f to %.2f\n\n', metrics_motor_orig.K_factor, metrics_motor_filt.K_factor);

fprintf('================================================================================\n');
fprintf('All visualizations complete! Review the 7 figure windows.\n');
fprintf('================================================================================\n\n');
