% TEST_POWER_QUALITY_METRICS - Test bench for comprehensive power quality metrics
%
% DESCRIPTION:
%   Tests the power_quality_metrics.m function by:
%   1. Analyzing all power quality scenarios
%   2. Verifying metric calculations
%   3. Comparing metrics across scenarios
%   4. Generating comprehensive reports
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
fprintf('POWER QUALITY METRICS TEST BENCH\n');
fprintf('========================================\n\n');

% =======================================================================
% TEST ALL SCENARIOS
% =======================================================================

scenarios = {'ideal', 'led_lighting', 'motor_drive', 'data_center', 'heavy_distortion'};
num_scenarios = length(scenarios);

% Parameters
fs = 10000;
f0 = 50;
N = fs/f0;
V_rms = 230;

% Storage for all metrics
all_metrics = cell(num_scenarios, 1);

fprintf('Analyzing all power quality scenarios...\n\n');

for i = 1:num_scenarios
    fprintf('--- Scenario %d/%d: %s ---\n', i, num_scenarios, upper(scenarios{i}));

    % Generate signal
    [signal, t, info] = generate_distorted_signal(scenarios{i}, fs, f0, N, V_rms);

    % Perform DTFS analysis
    [X_k, mag, phase, freq] = calculate_dtfs(signal, fs, N);

    % Calculate comprehensive metrics
    metrics = power_quality_metrics(signal, mag, fs, f0);

    % Store results
    all_metrics{i} = metrics;
    all_metrics{i}.scenario_name = info.scenario_name;
end

% =======================================================================
% COMPARATIVE ANALYSIS
% =======================================================================

fprintf('\n========================================\n');
fprintf('COMPARATIVE METRICS TABLE\n');
fprintf('========================================\n\n');

% Table 1: Distortion Metrics
fprintf('HARMONIC DISTORTION COMPARISON:\n');
fprintf('%-20s | THD-F(%%) | THD-R(%%) | Even(%%) | Odd(%%)\n', 'Scenario');
fprintf('---------------------|---------|---------|---------|--------\n');

for i = 1:num_scenarios
    m = all_metrics{i};
    name = m.scenario_name(1:min(20,end));
    fprintf('%-20s | %7.2f | %7.2f | %7.2f | %6.2f\n', ...
        name, m.THD_F, m.THD_R, m.even_harmonic_distortion, m.odd_harmonic_distortion);
end
fprintf('\n');

% Table 2: Waveform Shape Metrics
fprintf('WAVEFORM CHARACTERISTICS:\n');
fprintf('%-20s | Crest F | Form F  | Dist F  \n', 'Scenario');
fprintf('---------------------|---------|---------|----------\n');

for i = 1:num_scenarios
    m = all_metrics{i};
    name = m.scenario_name(1:min(20,end));
    fprintf('%-20s | %7.4f | %7.4f | %8.4f\n', ...
        name, m.crest_factor, m.form_factor, m.distortion_factor);
end
fprintf('\n');

% Table 3: Power Factor Metrics
fprintf('POWER FACTOR ANALYSIS:\n');
fprintf('%-20s | Disp PF | Dist PF | True PF\n', 'Scenario');
fprintf('---------------------|---------|---------|--------\n');

for i = 1:num_scenarios
    m = all_metrics{i};
    name = m.scenario_name(1:min(20,end));
    fprintf('%-20s | %7.4f | %7.4f | %7.4f\n', ...
        name, m.displacement_PF, m.distortion_PF, m.true_PF);
end
fprintf('\n');

% Table 4: Equipment Stress Factors
fprintf('EQUIPMENT STRESS FACTORS:\n');
fprintf('%-20s | K-Factor |   TIF   | I²t\n', 'Scenario');
fprintf('---------------------|----------|---------|----------\n');

for i = 1:num_scenarios
    m = all_metrics{i};
    name = m.scenario_name(1:min(20,end));
    fprintf('%-20s | %8.2f | %7.4f | %8.2f\n', ...
        name, m.K_factor, m.TIF, m.IT_product);
end
fprintf('\n');

% Table 5: Quality Ratings
fprintf('IEEE 519 QUALITY ASSESSMENT:\n');
fprintf('%-20s | Rating         | Compliant | Dominant\n', 'Scenario');
fprintf('---------------------|----------------|-----------|----------\n');

for i = 1:num_scenarios
    m = all_metrics{i};
    name = m.scenario_name(1:min(20,end));
    compliance = 'Yes';
    if ~m.compliance_IEEE519
        compliance = 'No ';
    end
    fprintf('%-20s | %-14s | %9s | %dth\n', ...
        name, m.quality_rating, compliance, m.dominant_harmonic);
end
fprintf('\n');

% =======================================================================
% DETAILED HARMONIC CONTENT FOR SELECTED SCENARIOS
% =======================================================================

fprintf('========================================\n');
fprintf('DETAILED HARMONIC CONTENT\n');
fprintf('========================================\n\n');

% Show detailed harmonics for Motor Drive scenario
fprintf('Motor Drive Scenario - Top 10 Harmonics:\n');
fprintf('Order | Freq (Hz) | Magnitude (V) | %% of Fund | IHD (%%)\n');
fprintf('------|-----------|---------------|-----------|--------\n');

m = all_metrics{3};  % Motor drive
if ~isempty(m.harmonic_content)
    for h = 1:min(10, size(m.harmonic_content, 1))
        row = m.harmonic_content(h, :);
        fprintf('%5.0f | %9.0f | %13.4f | %9.2f | %6.2f\n', ...
            row(1), row(2), row(3), row(4), row(5));
    end
end
fprintf('\n');

% =======================================================================
% VERIFICATION TESTS
% =======================================================================

fprintf('========================================\n');
fprintf('METRIC VERIFICATION TESTS\n');
fprintf('========================================\n\n');

% Test 1: Ideal signal should have perfect metrics
fprintf('TEST 1: Ideal Signal Verification\n');
fprintf('----------------------------------------\n');
m_ideal = all_metrics{1};

test1_pass = true;

% Check THD is zero
if abs(m_ideal.THD_F) < 0.01
    fprintf('THD = %.4f%% : PASS (expected ~0%%)\n', m_ideal.THD_F);
else
    fprintf('THD = %.4f%% : FAIL\n', m_ideal.THD_F);
    test1_pass = false;
end

% Check Crest Factor is sqrt(2)
expected_cf = sqrt(2);
if abs(m_ideal.crest_factor - expected_cf) < 0.01
    fprintf('Crest Factor = %.4f : PASS (expected %.4f)\n', m_ideal.crest_factor, expected_cf);
else
    fprintf('Crest Factor = %.4f : FAIL\n', m_ideal.crest_factor);
    test1_pass = false;
end

% Check Power Factor is 1.0
if abs(m_ideal.true_PF - 1.0) < 0.01
    fprintf('True PF = %.4f : PASS (expected 1.0)\n', m_ideal.true_PF);
else
    fprintf('True PF = %.4f : FAIL\n', m_ideal.true_PF);
    test1_pass = false;
end

% Check K-Factor is 1.0
if abs(m_ideal.K_factor - 1.0) < 0.1
    fprintf('K-Factor = %.2f : PASS (expected 1.0)\n', m_ideal.K_factor);
else
    fprintf('K-Factor = %.2f : FAIL\n', m_ideal.K_factor);
    test1_pass = false;
end

if test1_pass
    fprintf('\nTEST 1: PASSED\n\n');
else
    fprintf('\nTEST 1: FAILED\n\n');
end

% Test 2: Higher THD should decrease power factor
fprintf('TEST 2: Power Factor Correlation Test\n');
fprintf('----------------------------------------\n');

test2_pass = true;

% Compare ideal vs heavy distortion
pf_ideal = all_metrics{1}.true_PF;
pf_heavy = all_metrics{5}.true_PF;

fprintf('Ideal PF:      %.4f (THD = %.2f%%)\n', pf_ideal, all_metrics{1}.THD_F);
fprintf('Heavy Dist PF: %.4f (THD = %.2f%%)\n', pf_heavy, all_metrics{5}.THD_F);

if pf_ideal > pf_heavy
    fprintf('PF decreases with THD: PASS\n');
else
    fprintf('PF should decrease with THD: FAIL\n');
    test2_pass = false;
end

if test2_pass
    fprintf('\nTEST 2: PASSED\n\n');
else
    fprintf('\nTEST 2: FAILED\n\n');
end

% Test 3: K-Factor increases with harmonics
fprintf('TEST 3: K-Factor Scaling Test\n');
fprintf('----------------------------------------\n');

test3_pass = true;

k_ideal = all_metrics{1}.K_factor;
k_led = all_metrics{2}.K_factor;
k_motor = all_metrics{3}.K_factor;
k_heavy = all_metrics{5}.K_factor;

fprintf('Ideal:      K = %.2f\n', k_ideal);
fprintf('LED:        K = %.2f\n', k_led);
fprintf('Motor:      K = %.2f\n', k_motor);
fprintf('Heavy Dist: K = %.2f\n', k_heavy);

if k_ideal < k_led && k_led < k_heavy && k_motor < k_heavy
    fprintf('K-Factor increases with distortion: PASS\n');
else
    fprintf('K-Factor should increase with distortion: FAIL\n');
    test3_pass = false;
end

if test3_pass
    fprintf('\nTEST 3: PASSED\n\n');
else
    fprintf('\nTEST 3: FAILED\n\n');
end

% =======================================================================
% VISUALIZATION
% =======================================================================

fprintf('Generating comprehensive visualization...\n');

figure('Name', 'Power Quality Metrics Comparison', 'Position', [50, 50, 1600, 1000]);

% Subplot 1: THD Comparison
subplot(3, 3, 1);
scenario_names = {'Ideal', 'LED', 'Motor', 'Data Ctr', 'Heavy'};
thd_values = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    thd_values(i) = all_metrics{i}.THD_F;
end

bar(thd_values, 'FaceColor', [0.3 0.6 0.9]);
hold on;
yline(5, 'g--', 'LineWidth', 2);
yline(8, 'r--', 'LineWidth', 2);
set(gca, 'XTickLabel', scenario_names);
ylabel('THD-F (%)');
title('Total Harmonic Distortion');
grid on;

% Subplot 2: Crest Factor Comparison
subplot(3, 3, 2);
cf_values = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    cf_values(i) = all_metrics{i}.crest_factor;
end

bar(cf_values, 'FaceColor', [0.9 0.6 0.3]);
hold on;
yline(1.414, 'g--', 'LineWidth', 2, 'Label', 'Ideal');
set(gca, 'XTickLabel', scenario_names);
ylabel('Crest Factor');
title('Crest Factor (Peak Stress)');
grid on;

% Subplot 3: Power Factor Comparison
subplot(3, 3, 3);
pf_values = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    pf_values(i) = all_metrics{i}.true_PF;
end

bar(pf_values, 'FaceColor', [0.3 0.9 0.6]);
hold on;
yline(0.9, 'r--', 'LineWidth', 2, 'Label', 'Min Acceptable');
set(gca, 'XTickLabel', scenario_names);
ylabel('Power Factor');
title('True Power Factor');
ylim([0.85 1.05]);
grid on;

% Subplot 4: K-Factor Comparison
subplot(3, 3, 4);
k_values = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    k_values(i) = all_metrics{i}.K_factor;
end

bar(k_values, 'FaceColor', [0.9 0.3 0.6]);
hold on;
yline(4, 'y--', 'LineWidth', 1.5, 'Label', 'K-4');
yline(13, 'r--', 'LineWidth', 1.5, 'Label', 'K-13');
set(gca, 'XTickLabel', scenario_names);
ylabel('K-Factor');
title('Transformer K-Factor');
grid on;

% Subplot 5: Even vs Odd Harmonics
subplot(3, 3, 5);
even_values = zeros(num_scenarios, 1);
odd_values = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    even_values(i) = all_metrics{i}.even_harmonic_distortion;
    odd_values(i) = all_metrics{i}.odd_harmonic_distortion;
end

x = 1:num_scenarios;
width = 0.35;
bar(x - width/2, even_values, width, 'FaceColor', [0.5 0.5 0.9], 'DisplayName', 'Even');
hold on;
bar(x + width/2, odd_values, width, 'FaceColor', [0.9 0.5 0.5], 'DisplayName', 'Odd');
set(gca, 'XTickLabel', scenario_names);
ylabel('THD (%)');
title('Even vs Odd Harmonic Distortion');
legend('Location', 'northwest');
grid on;

% Subplot 6: Distortion Factor
subplot(3, 3, 6);
df_values = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    df_values(i) = all_metrics{i}.distortion_factor;
end

bar(df_values, 'FaceColor', [0.6 0.3 0.9]);
set(gca, 'XTickLabel', scenario_names);
ylabel('Distortion Factor');
title('Distortion Factor (Harmonic/Total)');
grid on;

% Subplot 7: RMS Voltage Comparison
subplot(3, 3, 7);
rms_values = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    rms_values(i) = all_metrics{i}.signal_rms;
end

bar(rms_values, 'FaceColor', [0.3 0.9 0.9]);
hold on;
yline(230, 'r--', 'LineWidth', 2, 'Label', 'Nominal');
set(gca, 'XTickLabel', scenario_names);
ylabel('RMS Voltage (V)');
title('Total RMS Voltage');
grid on;

% Subplot 8: TIF Comparison
subplot(3, 3, 8);
tif_values = zeros(num_scenarios, 1);
for i = 1:num_scenarios
    tif_values(i) = all_metrics{i}.TIF;
end

bar(tif_values, 'FaceColor', [0.9 0.9 0.3]);
set(gca, 'XTickLabel', scenario_names);
ylabel('TIF');
title('Telephone Influence Factor');
grid on;

% Subplot 9: Compliance Summary
subplot(3, 3, 9);
axis off;

y_pos = 0.95;
text(0.1, y_pos, 'IEEE 519 COMPLIANCE SUMMARY', 'FontSize', 11, 'FontWeight', 'bold');
y_pos = y_pos - 0.12;

for i = 1:num_scenarios
    m = all_metrics{i};

    if m.compliance_IEEE519
        color = [0 0.6 0];
        status = 'PASS';
    else
        color = [0.8 0 0];
        status = 'FAIL';
    end

    text(0.1, y_pos, sprintf('%s:', scenario_names{i}), 'FontSize', 9);
    text(0.5, y_pos, sprintf('%s (THD=%.1f%%)', status, m.THD_F), ...
        'FontSize', 9, 'Color', color, 'FontWeight', 'bold');

    y_pos = y_pos - 0.12;
end

sgtitle('Comprehensive Power Quality Metrics - Multi-Scenario Analysis', ...
    'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, 'figures/power_quality_metrics_comparison.png');

fprintf('========================================\n');
fprintf('ALL TESTS COMPLETED\n');
fprintf('Figure saved to: figures/power_quality_metrics_comparison.png\n');
fprintf('========================================\n');
