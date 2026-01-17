% TEST_GENERATE_SIGNALS - Test bench for signal generation function
%
% DESCRIPTION:
%   Tests the generate_distorted_signal.m function by:
%   1. Generating all predefined power quality scenarios
%   2. Verifying signal characteristics (RMS, peak, harmonics)
%   3. Performing DTFS analysis on each scenario
%   4. Visualizing time and frequency domain representations
%   5. Comparing expected vs actual THD
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
fprintf('SIGNAL GENERATION TEST BENCH\n');
fprintf('========================================\n\n');

% =======================================================================
% DEFINE TEST SCENARIOS
% =======================================================================

scenarios = {'ideal', 'led_lighting', 'motor_drive', 'data_center', 'heavy_distortion'};
num_scenarios = length(scenarios);

% Common parameters
fs = 10000;   % 10 kHz sampling
f0 = 50;      % 50 Hz fundamental
N = fs/f0;    % 200 samples per period
V_rms = 230;  % 230V RMS

% Storage for all scenarios
all_signals = cell(num_scenarios, 1);
all_time = cell(num_scenarios, 1);
all_info = cell(num_scenarios, 1);
all_dtfs = cell(num_scenarios, 1);

% =======================================================================
% GENERATE ALL SCENARIOS
% =======================================================================

fprintf('Generating all power quality scenarios...\n\n');

for i = 1:num_scenarios
    fprintf('--- Scenario %d/%d: %s ---\n', i, num_scenarios, upper(scenarios{i}));

    % Generate signal
    [signal, t, info] = generate_distorted_signal(scenarios{i}, fs, f0, N, V_rms);

    % Store results
    all_signals{i} = signal;
    all_time{i} = t;
    all_info{i} = info;

    % Perform DTFS analysis
    [X_k, mag, phase, freq] = calculate_dtfs(signal, fs, N);

    % Store DTFS results
    all_dtfs{i}.X_k = X_k;
    all_dtfs{i}.magnitude = mag;
    all_dtfs{i}.phase = phase;
    all_dtfs{i}.frequencies = freq;

    fprintf('\n');
end

% =======================================================================
% CALCULATE ACTUAL THD FOR ALL SCENARIOS
% =======================================================================

fprintf('========================================\n');
fprintf('THD VERIFICATION\n');
fprintf('========================================\n\n');

fprintf('Scenario             | Expected THD | Actual THD | Error\n');
fprintf('---------------------|--------------|------------|-------\n');

for i = 1:num_scenarios
    % Get DTFS magnitude spectrum
    mag = all_dtfs{i}.magnitude;

    % Calculate THD
    % THD = sqrt(V2^2 + V3^2 + ... + Vn^2) / V1 * 100%
    V1 = mag(2);  % Fundamental (k=1, index 2)

    % Sum of squares of harmonics (indices 3 onwards represent harmonics)
    % We'll consider harmonics up to 20th (index 22)
    max_harmonic_idx = min(22, length(mag));
    harmonic_sum_sq = sum(mag(3:max_harmonic_idx).^2);

    actual_THD = sqrt(harmonic_sum_sq) / V1 * 100;
    expected_THD = all_info{i}.expected_THD;

    % Calculate error
    thd_error = abs(actual_THD - expected_THD);

    fprintf('%-20s | %8.2f%%   | %7.2f%%  | %.2f%%\n', ...
        all_info{i}.scenario_name(1:min(20,end)), expected_THD, actual_THD, thd_error);

    % Store actual THD
    all_info{i}.actual_THD = actual_THD;
end

fprintf('========================================\n\n');

% =======================================================================
% SIGNAL CHARACTERISTICS COMPARISON
% =======================================================================

fprintf('========================================\n');
fprintf('SIGNAL CHARACTERISTICS\n');
fprintf('========================================\n\n');

fprintf('Scenario             | RMS (V) | Peak (V) | Crest Factor\n');
fprintf('---------------------|---------|----------|--------------\n');

for i = 1:num_scenarios
    name = all_info{i}.scenario_name(1:min(20,end));
    rms_val = all_info{i}.signal_rms;
    peak_val = all_info{i}.signal_peak;
    cf = all_info{i}.crest_factor;

    fprintf('%-20s | %7.2f | %8.2f | %12.4f\n', name, rms_val, peak_val, cf);
end

fprintf('========================================\n\n');

% =======================================================================
% HARMONIC SPECTRUM ANALYSIS
% =======================================================================

fprintf('========================================\n');
fprintf('HARMONIC SPECTRUM COMPARISON\n');
fprintf('========================================\n\n');

for i = 1:num_scenarios
    fprintf('--- %s ---\n', all_info{i}.scenario_name);
    fprintf('Harmonic | Expected Mag | Actual Mag | Error\n');
    fprintf('---------|--------------|------------|---------\n');

    mag = all_dtfs{i}.magnitude;
    info = all_info{i};

    for j = 1:length(info.harmonics)
        h = info.harmonics(j);
        idx = h + 1;  % MATLAB indexing

        expected_mag = info.magnitudes(j);
        actual_mag = mag(idx);
        error_pct = abs(actual_mag - expected_mag) / expected_mag * 100;

        fprintf('   %2d    | %10.2f V | %9.2f V | %6.3f%%\n', ...
            h, expected_mag, actual_mag, error_pct);
    end

    fprintf('\n');
end

fprintf('========================================\n\n');

% =======================================================================
% VISUALIZATION
% =======================================================================

fprintf('Generating comprehensive visualization...\n');

% Figure 1: Time Domain Comparison (All Scenarios)
figure('Name', 'Power Quality Scenarios - Time Domain', 'Position', [50, 50, 1400, 900]);

for i = 1:num_scenarios
    subplot(3, 2, i);
    plot(all_time{i}*1000, all_signals{i}, 'LineWidth', 1.5);
    grid on;
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title(sprintf('%s\nTHD: %.1f%%, CF: %.3f', ...
        all_info{i}.scenario_name, all_info{i}.actual_THD, all_info{i}.crest_factor));
    xlim([0 40]);  % Show 2 cycles (40 ms)
end

sgtitle('Power Quality Scenarios - Time Domain Comparison', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'figures/scenarios_time_domain.png');

% Figure 2: Frequency Domain Comparison (All Scenarios)
figure('Name', 'Power Quality Scenarios - Frequency Domain', 'Position', [100, 100, 1400, 900]);

for i = 1:num_scenarios
    subplot(3, 2, i);
    mag = all_dtfs{i}.magnitude;
    freq = all_dtfs{i}.frequencies;

    % Plot harmonics up to 1 kHz
    max_idx = find(freq <= 1000, 1, 'last');
    stem(freq(1:max_idx), mag(1:max_idx), 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (V)');
    title(sprintf('%s\nTHD: %.1f%%', ...
        all_info{i}.scenario_name, all_info{i}.actual_THD));
    xlim([0 1000]);
end

sgtitle('Power Quality Scenarios - Harmonic Spectrum', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'figures/scenarios_frequency_domain.png');

% Figure 3: Detailed Comparison - First 4 Scenarios
figure('Name', 'Detailed Scenario Comparison', 'Position', [150, 150, 1600, 1000]);

for i = 1:min(4, num_scenarios)
    % Time domain
    subplot(4, 3, (i-1)*3 + 1);
    plot(all_time{i}*1000, all_signals{i}, 'LineWidth', 1.8, 'Color', [0 0.4 0.8]);
    grid on;
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title(sprintf('%s - Waveform', all_info{i}.scenario_name));
    xlim([0 20]);  % One cycle

    % Frequency spectrum (bar plot)
    subplot(4, 3, (i-1)*3 + 2);
    mag = all_dtfs{i}.magnitude;
    harmonic_indices = [2, 4, 6, 8, 10, 12, 14, 16];  % 1st, 3rd, 5th... 15th harmonics
    harmonic_numbers = [1, 3, 5, 7, 9, 11, 13, 15];

    bar(harmonic_numbers, mag(harmonic_indices), 'FaceColor', [0.8 0.2 0.2]);
    grid on;
    xlabel('Harmonic Order');
    ylabel('Magnitude (V)');
    title('Harmonic Content');
    xlim([0 16]);

    % THD pie chart (fundamental vs harmonics)
    subplot(4, 3, (i-1)*3 + 3);
    V1 = mag(2);
    V_harmonics = sqrt(sum(mag(3:22).^2));
    pie_data = [V1, V_harmonics];
    pie_labels = {sprintf('Fundamental\n%.1f V', V1), sprintf('Harmonics\n%.1f V', V_harmonics)};

    pie(pie_data, pie_labels);
    title(sprintf('THD: %.1f%%', all_info{i}.actual_THD));
    colormap([0.2 0.8 0.2; 0.8 0.2 0.2]);
end

sgtitle('Detailed Power Quality Analysis - Four Main Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'figures/scenarios_detailed_comparison.png');

% Figure 4: THD Comparison Bar Chart
figure('Name', 'THD Comparison', 'Position', [200, 200, 1000, 600]);

scenario_names_short = cell(num_scenarios, 1);
thd_values = zeros(num_scenarios, 1);

for i = 1:num_scenarios
    % Extract first word or abbreviate
    full_name = all_info{i}.scenario_name;
    if contains(full_name, 'Ideal')
        scenario_names_short{i} = 'Ideal';
    elseif contains(full_name, 'LED')
        scenario_names_short{i} = 'LED Lighting';
    elseif contains(full_name, 'Motor')
        scenario_names_short{i} = 'Motor Drive';
    elseif contains(full_name, 'Data')
        scenario_names_short{i} = 'Data Center';
    elseif contains(full_name, 'Heavy')
        scenario_names_short{i} = 'Heavy Distortion';
    else
        scenario_names_short{i} = full_name(1:min(15,end));
    end

    thd_values(i) = all_info{i}.actual_THD;
end

bar(thd_values, 'FaceColor', [0.3 0.6 0.9]);
hold on;

% Add IEEE standard reference lines
yline(5, 'g--', 'LineWidth', 2, 'Label', 'IEEE Excellent (<5%)');
yline(8, 'y--', 'LineWidth', 2, 'Label', 'IEEE Good (<8%)');

set(gca, 'XTickLabel', scenario_names_short);
ylabel('Total Harmonic Distortion (%)');
xlabel('Power Quality Scenario');
title('THD Comparison Across Scenarios (IEEE Standard Reference)');
grid on;
ylim([0 max(thd_values)*1.2]);

% Add value labels on bars
for i = 1:num_scenarios
    text(i, thd_values(i) + 1, sprintf('%.1f%%', thd_values(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

saveas(gcf, 'figures/thd_comparison.png');

fprintf('========================================\n');
fprintf('ALL TESTS COMPLETED SUCCESSFULLY\n');
fprintf('Figures saved to figures/ folder:\n');
fprintf('  - scenarios_time_domain.png\n');
fprintf('  - scenarios_frequency_domain.png\n');
fprintf('  - scenarios_detailed_comparison.png\n');
fprintf('  - thd_comparison.png\n');
fprintf('========================================\n');
