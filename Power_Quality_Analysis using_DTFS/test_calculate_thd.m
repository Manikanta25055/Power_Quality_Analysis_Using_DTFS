% TEST_CALCULATE_THD - Test bench for THD calculation function
%
% DESCRIPTION:
%   Tests the calculate_thd.m function by:
%   1. Generating signals with known harmonic content
%   2. Calculating THD using both methods (fundamental and RMS)
%   3. Comparing calculated THD with theoretical values
%   4. Testing with all power quality scenarios
%   5. Validating IEEE standard compliance checks
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
fprintf('THD CALCULATION TEST BENCH\n');
fprintf('========================================\n\n');

% =======================================================================
% TEST CASE 1: THEORETICAL VERIFICATION
% =======================================================================

fprintf('TEST CASE 1: Theoretical THD Verification\n');
fprintf('----------------------------------------\n');

% Create a signal with known harmonics:
% Fundamental: 100V (peak)
% 3rd Harmonic: 30V (30% of fundamental)
% 5th Harmonic: 20V (20% of fundamental)
%
% Theoretical THD = sqrt(0.3^2 + 0.2^2) * 100% = 36.06%

V1 = 100;  % Fundamental
V3 = 30;   % 3rd harmonic
V5 = 20;   % 5th harmonic

% Theoretical THD calculation
theoretical_THD = sqrt((V3/V1)^2 + (V5/V1)^2) * 100;

% Create magnitude spectrum (simulating DTFS output)
mag_spectrum = zeros(12, 1);
mag_spectrum(1) = 0;     % DC component
mag_spectrum(2) = V1;    % Fundamental (k=1)
mag_spectrum(4) = V3;    % 3rd harmonic (k=3)
mag_spectrum(6) = V5;    % 5th harmonic (k=5)

% Calculate THD
[THD_pct, THD_ratio, details] = calculate_thd(mag_spectrum, 'fundamental');

% Verify results
fprintf('THEORETICAL THD: %.4f%%\n', theoretical_THD);
fprintf('CALCULATED THD:  %.4f%%\n', THD_pct);
error_pct = abs(THD_pct - theoretical_THD);
fprintf('ERROR: %.6f%%\n', error_pct);

if error_pct < 0.01
    fprintf('\nTEST CASE 1: PASSED\n\n');
else
    fprintf('\nTEST CASE 1: FAILED (Error = %.4f%%)\n\n', error_pct);
end

% =======================================================================
% TEST CASE 2: ALL POWER QUALITY SCENARIOS
% =======================================================================

fprintf('TEST CASE 2: Power Quality Scenarios THD Analysis\n');
fprintf('----------------------------------------\n\n');

scenarios = {'ideal', 'led_lighting', 'motor_drive', 'data_center', 'heavy_distortion'};
num_scenarios = length(scenarios);

% Parameters
fs = 10000;
f0 = 50;
N = fs/f0;
V_rms = 230;

% Storage
thd_results = zeros(num_scenarios, 1);
quality_ratings = cell(num_scenarios, 1);

fprintf('Scenario             | THD (%%)|   Quality Rating\n');
fprintf('---------------------|--------|---------------------------\n');

for i = 1:num_scenarios
    % Generate signal
    [signal, t, info] = generate_distorted_signal(scenarios{i}, fs, f0, N, V_rms);

    % Perform DTFS analysis
    [X_k, mag, phase, freq] = calculate_dtfs(signal, fs, N);

    % Calculate THD
    [THD_pct, THD_ratio, details] = calculate_thd(mag, 'fundamental');

    % Store results
    thd_results(i) = THD_pct;
    quality_ratings{i} = details.quality_rating;

    % Display results
    scenario_name = info.scenario_name(1:min(20, end));
    fprintf('%-20s | %6.2f | %s\n', scenario_name, THD_pct, details.quality_rating);
end

fprintf('\nTEST CASE 2: COMPLETED\n\n');

% =======================================================================
% TEST CASE 3: IEEE STANDARD COMPLIANCE
% =======================================================================

fprintf('TEST CASE 3: IEEE Standard Compliance Check\n');
fprintf('----------------------------------------\n');

% Test boundary conditions
test_thd_values = [3.0, 4.9, 5.0, 6.5, 7.9, 8.0, 10.0, 15.0];
expected_ratings = {'EXCELLENT', 'EXCELLENT', 'GOOD', 'GOOD', 'GOOD', ...
                    'POOR (Action Required)', 'POOR (Action Required)', 'POOR (Action Required)'};

fprintf('\nTHD Value | Expected Rating          | Actual Rating            | Status\n');
fprintf('----------|--------------------------|--------------------------|--------\n');

all_passed = true;

for i = 1:length(test_thd_values)
    % Create magnitude spectrum with specific THD
    target_thd = test_thd_values(i) / 100;  % Convert to ratio

    V1 = 100;  % Fundamental
    V_harmonics = target_thd * V1;  % Total harmonic content

    % Create spectrum with single harmonic for simplicity
    mag_spectrum = zeros(6, 1);
    mag_spectrum(2) = V1;
    mag_spectrum(4) = V_harmonics;  % Put all harmonic content in 3rd harmonic

    % Calculate THD
    [THD_pct, ~, details] = calculate_thd(mag_spectrum, 'fundamental');

    % Check rating
    passed = strcmp(details.quality_rating, expected_ratings{i});
    all_passed = all_passed && passed;

    status = 'PASS';
    if ~passed
        status = 'FAIL';
    end

    fprintf('%6.1f%%   | %-24s | %-24s | %s\n', ...
        test_thd_values(i), expected_ratings{i}, details.quality_rating, status);
end

if all_passed
    fprintf('\nTEST CASE 3: PASSED (All IEEE ratings correct)\n\n');
else
    fprintf('\nTEST CASE 3: FAILED\n\n');
end

% =======================================================================
% TEST CASE 4: METHOD COMPARISON (THD-F vs THD-R)
% =======================================================================

fprintf('TEST CASE 4: THD Method Comparison\n');
fprintf('----------------------------------------\n');

% Generate a distorted signal
[signal, t, info] = generate_distorted_signal('motor_drive', fs, f0, N, V_rms);

% Perform DTFS
[X_k, mag, phase, freq] = calculate_dtfs(signal, fs, N);

% Calculate using both methods
[THD_F, ratio_F, details_F] = calculate_thd(mag, 'fundamental');
[THD_R, ratio_R, details_R] = calculate_thd(mag, 'rms');

fprintf('Method Comparison for Motor Drive Scenario:\n');
fprintf('  THD-F (Fundamental): %.4f%%\n', THD_F);
fprintf('  THD-R (RMS):         %.4f%%\n', THD_R);
fprintf('  Difference:          %.4f%%\n', abs(THD_F - THD_R));
fprintf('\nNote: THD-F is always >= THD-R\n');
fprintf('IEEE 519 standard uses THD-F (Fundamental reference)\n');

if THD_F >= THD_R
    fprintf('\nTEST CASE 4: PASSED (THD-F >= THD-R as expected)\n\n');
else
    fprintf('\nTEST CASE 4: FAILED (THD-F should be >= THD-R)\n\n');
end

% =======================================================================
% TEST CASE 5: INDIVIDUAL HARMONIC DISTORTION (IHD)
% =======================================================================

fprintf('TEST CASE 5: Individual Harmonic Distortion Analysis\n');
fprintf('----------------------------------------\n');

% Use LED lighting scenario (3rd harmonic dominant)
[signal, t, info] = generate_distorted_signal('led_lighting', fs, f0, N, V_rms);
[X_k, mag, phase, freq] = calculate_dtfs(signal, fs, N);
[THD_pct, THD_ratio, details] = calculate_thd(mag, 'fundamental');

fprintf('LED Lighting Scenario - Top 5 Harmonics:\n');
fprintf('  Harmonic | IHD (%%)  | Expected (~)\n');
fprintf('  ---------|----------|-------------\n');

% Expected values from generate_distorted_signal LED scenario
expected_ihd = [18.0, 8.0, 4.0, 2.0];  % 3rd, 5th, 7th, 9th

% Check 3rd, 5th, 7th, 9th harmonics
harmonics_to_check = [3, 5, 7, 9];

for i = 1:length(harmonics_to_check)
    h = harmonics_to_check(i);
    idx = h;  % Index in details.individual_thd (starts from 2nd harmonic)

    if idx <= length(details.individual_thd)
        actual_ihd = details.individual_thd(idx);
        expected = expected_ihd(i);

        fprintf('    %2d     | %7.3f  | %7.1f\n', h+1, actual_ihd, expected);
    end
end

fprintf('\nTEST CASE 5: COMPLETED\n\n');

% =======================================================================
% VISUALIZATION
% =======================================================================

fprintf('Generating visualization...\n');

figure('Name', 'THD Analysis Results', 'Position', [100, 100, 1400, 800]);

% Subplot 1: THD Comparison Across Scenarios
subplot(2, 3, 1);
scenario_names_short = {'Ideal', 'LED', 'Motor', 'Data Ctr', 'Heavy'};
colors = [0.2 0.8 0.2; 0.9 0.9 0.2; 0.9 0.6 0.2; 0.9 0.6 0.2; 0.9 0.2 0.2];

bar_handle = bar(thd_results);
bar_handle.FaceColor = 'flat';

for i = 1:num_scenarios
    if thd_results(i) < 5
        bar_handle.CData(i,:) = [0.2 0.8 0.2];  % Green
    elseif thd_results(i) < 8
        bar_handle.CData(i,:) = [0.9 0.9 0.2];  % Yellow
    else
        bar_handle.CData(i,:) = [0.9 0.2 0.2];  % Red
    end
end

hold on;
yline(5, 'g--', 'LineWidth', 2, 'Label', '5% (Excellent)');
yline(8, 'r--', 'LineWidth', 2, 'Label', '8% (Good/Poor)');

set(gca, 'XTickLabel', scenario_names_short);
ylabel('THD (%)');
title('THD Comparison - All Scenarios');
grid on;

% Subplot 2: THD-F vs THD-R Comparison
subplot(2, 3, 2);

% Calculate both for all scenarios
thd_f_all = zeros(num_scenarios, 1);
thd_r_all = zeros(num_scenarios, 1);

for i = 1:num_scenarios
    [signal, ~, ~] = generate_distorted_signal(scenarios{i}, fs, f0, N, V_rms);
    [X_k, mag, ~, ~] = calculate_dtfs(signal, fs, N);

    [thd_f_all(i), ~, ~] = calculate_thd(mag, 'fundamental');
    [thd_r_all(i), ~, ~] = calculate_thd(mag, 'rms');
end

x_pos = 1:num_scenarios;
bar_width = 0.35;

bar(x_pos - bar_width/2, thd_f_all, bar_width, 'FaceColor', [0.3 0.5 0.8], 'DisplayName', 'THD-F');
hold on;
bar(x_pos + bar_width/2, thd_r_all, bar_width, 'FaceColor', [0.8 0.5 0.3], 'DisplayName', 'THD-R');

set(gca, 'XTickLabel', scenario_names_short);
ylabel('THD (%)');
title('THD-F vs THD-R Comparison');
legend('Location', 'northwest');
grid on;

% Subplot 3: Dominant Harmonics
subplot(2, 3, 3);
axis off;
text(0.1, 0.95, 'Dominant Harmonics by Scenario', 'FontSize', 11, 'FontWeight', 'bold');
text(0.1, 0.85, 'LED Lighting: 3rd harmonic', 'FontSize', 9);
text(0.1, 0.75, 'Motor Drive: 5th harmonic', 'FontSize', 9);
text(0.1, 0.65, 'Data Center: 3rd harmonic', 'FontSize', 9);
text(0.1, 0.55, 'Heavy Distortion: 3rd harmonic', 'FontSize', 9);
text(0.1, 0.40, 'IEEE 519 Standard Limits:', 'FontSize', 10, 'FontWeight', 'bold');
text(0.1, 0.30, 'THD < 5%: Excellent', 'FontSize', 9, 'Color', [0 0.6 0]);
text(0.1, 0.22, 'THD 5-8%: Good', 'FontSize', 9, 'Color', [0.7 0.7 0]);
text(0.1, 0.14, 'THD > 8%: Poor', 'FontSize', 9, 'Color', [0.8 0 0]);

% Subplots 4-6: Harmonic spectrum for selected scenarios
scenarios_to_plot = {'led_lighting', 'motor_drive', 'heavy_distortion'};
scenario_titles = {'LED Lighting', 'Motor Drive', 'Heavy Distortion'};

for i = 1:3
    subplot(2, 3, 3+i);

    [signal, ~, ~] = generate_distorted_signal(scenarios_to_plot{i}, fs, f0, N, V_rms);
    [X_k, mag, ~, freq] = calculate_dtfs(signal, fs, N);
    [THD_pct, ~, details] = calculate_thd(mag, 'fundamental');

    % Plot first 15 harmonics
    harmonic_numbers = 1:15;
    harmonic_mags = mag(harmonic_numbers + 1);  % +1 for MATLAB indexing

    bar(harmonic_numbers, harmonic_mags, 'FaceColor', [0.4 0.6 0.9]);
    xlabel('Harmonic Order');
    ylabel('Magnitude (V)');
    title(sprintf('%s\nTHD: %.2f%% (%s)', scenario_titles{i}, THD_pct, details.quality_rating));
    grid on;
    xlim([0 16]);
end

sgtitle('Total Harmonic Distortion - Comprehensive Analysis', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'figures/thd_analysis_results.png');

fprintf('========================================\n');
fprintf('ALL THD TESTS COMPLETED\n');
fprintf('Figure saved to: figures/thd_analysis_results.png\n');
fprintf('========================================\n');
