function [signal, t, harmonic_info] = generate_distorted_signal(scenario, fs, f0, N, V_rms_fundamental)
% GENERATE_DISTORTED_SIGNAL - Creates realistic power signals with harmonic distortion
%
% DESCRIPTION:
%   Generates power signals with various types of harmonic distortion
%   representing real-world power quality scenarios. Each scenario mimics
%   actual non-linear loads found in electrical systems.
%
% SYNTAX:
%   [signal, t, harmonic_info] = generate_distorted_signal(scenario)
%   [signal, t, harmonic_info] = generate_distorted_signal(scenario, fs, f0, N, V_rms_fundamental)
%
% INPUTS:
%   scenario           - String specifying distortion scenario:
%                        'ideal'        - Pure sinusoid (no distortion)
%                        'led_lighting' - 3rd harmonic dominant (residential)
%                        'motor_drive'  - 5th & 7th harmonics (industrial)
%                        'data_center'  - Multiple odd harmonics (IT loads)
%                        'heavy_distortion' - Severe distortion (worst case)
%                        'custom'       - User-defined (prompts for inputs)
%
%   fs                 - (Optional) Sampling frequency in Hz (default: 10000 Hz)
%   f0                 - (Optional) Fundamental frequency in Hz (default: 50 Hz)
%   N                  - (Optional) Period length in samples (default: fs/f0)
%   V_rms_fundamental  - (Optional) RMS voltage of fundamental (default: 230 V)
%
% OUTPUTS:
%   signal        - Generated time-domain signal (voltage waveform)
%   t             - Time vector in seconds
%   harmonic_info - Structure containing harmonic information:
%                   .harmonics     - Harmonic numbers [1, 3, 5, 7, ...]
%                   .magnitudes    - Peak magnitudes in Volts
%                   .phases        - Phase angles in radians
%                   .percentages   - Harmonic as % of fundamental
%                   .scenario_name - Scenario description
%                   .expected_THD  - Approximate expected THD
%
% POWER QUALITY SCENARIOS:
%
%   1. IDEAL - Pure 50 Hz sinusoid
%      - THD: 0%
%      - Use case: Reference/baseline comparison
%
%   2. LED_LIGHTING - Residential/commercial LED loads
%      - 3rd Harmonic: 18% (dominant)
%      - 5th Harmonic: 8%
%      - 7th Harmonic: 4%
%      - THD: ~20%
%      - Cause: Single-phase rectifiers in LED drivers
%
%   3. MOTOR_DRIVE - Variable frequency drives (VFDs)
%      - 5th Harmonic: 20% (dominant)
%      - 7th Harmonic: 14%
%      - 11th Harmonic: 9%
%      - 13th Harmonic: 6%
%      - THD: ~28%
%      - Cause: 6-pulse rectifiers in motor drives
%
%   4. DATA_CENTER - Switch-mode power supplies (SMPS)
%      - 3rd Harmonic: 12%
%      - 5th Harmonic: 10%
%      - 7th Harmonic: 8%
%      - 9th Harmonic: 6%
%      - 11th Harmonic: 4%
%      - THD: ~20%
%      - Cause: High-frequency switching converters
%
%   5. HEAVY_DISTORTION - Severely polluted grid
%      - Multiple harmonics up to 15th
%      - THD: ~45%
%      - Use case: Worst-case analysis
%
% EXAMPLE USAGE:
%   % Generate LED lighting scenario
%   [signal, t, info] = generate_distorted_signal('led_lighting');
%   plot(t*1000, signal);
%   fprintf('Expected THD: %.1f%%\n', info.expected_THD);
%
%   % Generate motor drive with custom parameters
%   [signal, t, info] = generate_distorted_signal('motor_drive', 10000, 50, 200, 230);
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025
% COURSE: Digital Signal Processing (FISAC Assessment)
% PROJECT: Power Quality Analysis using DTFS

% =======================================================================
% INPUT VALIDATION AND DEFAULT VALUES
% =======================================================================

% Check if scenario is provided
if nargin < 1
    error('generate_distorted_signal:NoScenario', ...
        'Scenario must be specified. Options: ideal, led_lighting, motor_drive, data_center, heavy_distortion');
end

% Set default parameters
if nargin < 2 || isempty(fs)
    fs = 10000;  % 10 kHz sampling (standard for power analysis)
end

if nargin < 3 || isempty(f0)
    f0 = 50;     % 50 Hz fundamental (use 60 Hz for US systems)
end

if nargin < 4 || isempty(N)
    N = fs/f0;   % Samples per period
end

if nargin < 5 || isempty(V_rms_fundamental)
    V_rms_fundamental = 230;  % 230V RMS (standard in India/Europe)
end

% Validate scenario string
valid_scenarios = {'ideal', 'led_lighting', 'motor_drive', 'data_center', ...
                   'heavy_distortion', 'custom'};
if ~ismember(lower(scenario), valid_scenarios)
    error('generate_distorted_signal:InvalidScenario', ...
        'Invalid scenario. Choose: ideal, led_lighting, motor_drive, data_center, heavy_distortion, custom');
end

% =======================================================================
% CALCULATE FUNDAMENTAL PEAK VOLTAGE
% =======================================================================

V_peak_fundamental = V_rms_fundamental * sqrt(2);  % Convert RMS to peak

% =======================================================================
% GENERATE TIME VECTOR
% =======================================================================

t = (0:N-1)' / fs;  % Time vector for one period

% =======================================================================
% DEFINE HARMONIC CONTENT FOR EACH SCENARIO
% =======================================================================

scenario = lower(scenario);

switch scenario

    case 'ideal'
        % Pure sinusoid - no distortion
        harmonic_numbers = [1];
        harmonic_percentages = [100];
        harmonic_phases = [0];
        scenario_name = 'Ideal Power Signal (No Distortion)';
        expected_THD = 0.0;

    case 'led_lighting'
        % LED lighting with 3rd harmonic dominant
        % Typical of single-phase rectifiers in LED drivers
        harmonic_numbers = [1, 3, 5, 7, 9];
        harmonic_percentages = [100, 18, 8, 4, 2];
        harmonic_phases = [0, pi/6, -pi/4, pi/3, 0];
        scenario_name = 'LED Lighting Load (3rd Harmonic Dominant)';
        expected_THD = 20.5;

    case 'motor_drive'
        % Variable Frequency Drive (VFD) - 6-pulse rectifier
        % 5th and 7th harmonics dominant
        harmonic_numbers = [1, 5, 7, 11, 13];
        harmonic_percentages = [100, 20, 14, 9, 6];
        harmonic_phases = [0, -pi/4, pi/4, -pi/6, pi/6];
        scenario_name = 'Motor Drive / VFD (6-Pulse Rectifier)';
        expected_THD = 28.2;

    case 'data_center'
        % Switch-mode power supplies (SMPS)
        % Multiple odd harmonics
        harmonic_numbers = [1, 3, 5, 7, 9, 11, 13];
        harmonic_percentages = [100, 12, 10, 8, 6, 4, 3];
        harmonic_phases = [0, pi/8, -pi/6, pi/4, -pi/5, pi/7, -pi/8];
        scenario_name = 'Data Center / SMPS Loads';
        expected_THD = 20.4;

    case 'heavy_distortion'
        % Severely distorted grid - worst case
        % Multiple harmonics up to 15th
        harmonic_numbers = [1, 3, 5, 7, 9, 11, 13, 15];
        harmonic_percentages = [100, 25, 22, 18, 15, 12, 10, 8];
        harmonic_phases = [0, pi/6, -pi/4, pi/3, -pi/5, pi/4, -pi/6, pi/5];
        scenario_name = 'Heavy Distortion (Worst Case)';
        expected_THD = 45.6;

    case 'custom'
        % Custom scenario - prompt user for inputs
        fprintf('\n=== CUSTOM HARMONIC SCENARIO ===\n');
        harmonic_numbers = input('Enter harmonic numbers (e.g., [1 3 5 7]): ');
        harmonic_percentages = input('Enter percentages of fundamental (e.g., [100 15 10 5]): ');
        harmonic_phases = input('Enter phase angles in radians (e.g., [0 pi/6 -pi/4 pi/3]): ');
        scenario_name = 'Custom User-Defined Scenario';

        % Calculate expected THD
        THD_squared = sum((harmonic_percentages(2:end)/100).^2);
        expected_THD = sqrt(THD_squared) * 100;

end

% =======================================================================
% VALIDATE HARMONIC DATA CONSISTENCY
% =======================================================================

if length(harmonic_numbers) ~= length(harmonic_percentages) || ...
   length(harmonic_numbers) ~= length(harmonic_phases)
    error('generate_distorted_signal:DataMismatch', ...
        'Harmonic numbers, percentages, and phases must have same length');
end

% Ensure fundamental is first harmonic
if harmonic_numbers(1) ~= 1
    error('generate_distorted_signal:NoFundamental', ...
        'First harmonic must be fundamental (harmonic number 1)');
end

% =======================================================================
% CALCULATE HARMONIC MAGNITUDES (PEAK VALUES)
% =======================================================================

num_harmonics = length(harmonic_numbers);
harmonic_magnitudes = zeros(num_harmonics, 1);

for i = 1:num_harmonics
    harmonic_magnitudes(i) = V_peak_fundamental * (harmonic_percentages(i) / 100);
end

% =======================================================================
% GENERATE COMPOSITE SIGNAL
% =======================================================================

signal = zeros(N, 1);

fprintf('\n========================================\n');
fprintf('GENERATING DISTORTED POWER SIGNAL\n');
fprintf('========================================\n');
fprintf('Scenario: %s\n', scenario_name);
fprintf('Fundamental: %.2f Hz, %.2f V RMS (%.2f V peak)\n', ...
    f0, V_rms_fundamental, V_peak_fundamental);
fprintf('Sampling Frequency: %.0f Hz\n', fs);
fprintf('Period Length: %d samples\n', N);
fprintf('========================================\n');
fprintf('Harmonic Composition:\n');
fprintf('  Order | Frequency | Magnitude | Percentage | Phase\n');
fprintf('  ------|-----------|-----------|------------|-------\n');

for i = 1:num_harmonics
    h = harmonic_numbers(i);
    freq = h * f0;
    mag = harmonic_magnitudes(i);
    pct = harmonic_percentages(i);
    ph = harmonic_phases(i);

    % Add this harmonic component to signal
    signal = signal + mag * sin(2*pi*freq*t + ph);

    fprintf('   %2d   | %6.1f Hz | %6.2f V  |   %5.1f%%   | %+.3f\n', ...
        h, freq, mag, pct, ph);
end

fprintf('========================================\n');
fprintf('Expected THD: %.2f%%\n', expected_THD);
fprintf('Signal RMS: %.2f V\n', sqrt(mean(signal.^2)));
fprintf('Signal Peak: %.2f V\n', max(abs(signal)));
fprintf('Crest Factor: %.3f\n', max(abs(signal))/sqrt(mean(signal.^2)));
fprintf('========================================\n\n');

% =======================================================================
% PREPARE OUTPUT STRUCTURE
% =======================================================================

harmonic_info.harmonics = harmonic_numbers;
harmonic_info.magnitudes = harmonic_magnitudes;
harmonic_info.phases = harmonic_phases;
harmonic_info.percentages = harmonic_percentages;
harmonic_info.frequencies = harmonic_numbers * f0;
harmonic_info.scenario_name = scenario_name;
harmonic_info.expected_THD = expected_THD;
harmonic_info.fundamental_rms = V_rms_fundamental;
harmonic_info.fundamental_peak = V_peak_fundamental;
harmonic_info.signal_rms = sqrt(mean(signal.^2));
harmonic_info.signal_peak = max(abs(signal));
harmonic_info.crest_factor = max(abs(signal)) / sqrt(mean(signal.^2));

% =======================================================================
% END OF FUNCTION
% =======================================================================

end
