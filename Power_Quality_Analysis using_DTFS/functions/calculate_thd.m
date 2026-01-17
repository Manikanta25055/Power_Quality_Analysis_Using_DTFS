function [THD_percent, THD_ratio, harmonic_details] = calculate_thd(magnitude_spectrum, method)
% CALCULATE_THD - Calculates Total Harmonic Distortion from DTFS magnitude spectrum
%
% DESCRIPTION:
%   Computes Total Harmonic Distortion (THD) according to IEEE standards.
%   THD is the primary metric for assessing power quality and quantifying
%   the deviation of a power signal from an ideal sinusoid.
%
% SYNTAX:
%   [THD_percent, THD_ratio, harmonic_details] = calculate_thd(magnitude_spectrum)
%   [THD_percent, THD_ratio, harmonic_details] = calculate_thd(magnitude_spectrum, method)
%
% INPUTS:
%   magnitude_spectrum - Magnitude spectrum from DTFS analysis [N x 1 vector]
%                        magnitude_spectrum(1) = DC component
%                        magnitude_spectrum(2) = Fundamental (50 Hz)
%                        magnitude_spectrum(3+) = Harmonics
%
%   method             - (Optional) THD calculation method:
%                        'fundamental' - THD_F (default, IEEE standard)
%                                       THD = sqrt(sum(Vh^2)) / V1 * 100%
%                        'rms'         - THD_R (alternative definition)
%                                       THD = sqrt(sum(Vh^2)) / Vrms * 100%
%
% OUTPUTS:
%   THD_percent       - THD as percentage (%)
%   THD_ratio         - THD as ratio (0-1 range)
%   harmonic_details  - Structure with detailed harmonic analysis:
%                       .fundamental_mag     - Fundamental magnitude (V)
%                       .fundamental_rms     - Fundamental RMS (V)
%                       .total_rms          - Total signal RMS (V)
%                       .harmonic_rms       - RMS of all harmonics (V)
%                       .individual_THD     - Individual harmonic distortion (%)
%                       .harmonic_magnitudes - Individual harmonic magnitudes (V)
%                       .harmonic_percentages - Each harmonic as % of fundamental
%                       .dominant_harmonic   - Harmonic order with max magnitude
%                       .num_harmonics      - Number of harmonics analyzed
%
% IEEE STANDARDS:
%   According to IEEE 519-2014 (Harmonic Control in Electrical Power Systems):
%
%   THD Limits for Voltage:
%   - THD < 5%    : Excellent power quality
%   - THD 5-8%    : Good power quality (acceptable)
%   - THD > 8%    : Poor power quality (corrective action needed)
%
%   THD_F (Fundamental) is the standard definition:
%   THD_F = sqrt(V2^2 + V3^2 + V4^2 + ... + Vn^2) / V1 * 100%
%
%   where:
%   - V1 = Fundamental component magnitude
%   - V2, V3, ... Vn = Harmonic component magnitudes
%
% EXAMPLE USAGE:
%   % Perform DTFS analysis
%   [X_k, mag, phase, freq] = calculate_dtfs(signal, fs, N);
%
%   % Calculate THD
%   [THD_pct, THD_ratio, details] = calculate_thd(mag);
%
%   fprintf('Total Harmonic Distortion: %.2f%%\n', THD_pct);
%   fprintf('Dominant Harmonic: %d (%.2f V)\n', ...
%           details.dominant_harmonic, details.harmonic_magnitudes(details.dominant_harmonic));
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025
% COURSE: Digital Signal Processing (FISAC Assessment)
% PROJECT: Power Quality Analysis using DTFS

% =======================================================================
% INPUT VALIDATION
% =======================================================================

% Check if magnitude spectrum is provided
if nargin < 1
    error('calculate_thd:NoInput', 'Magnitude spectrum is required');
end

% Validate magnitude spectrum
if ~isnumeric(magnitude_spectrum) || ~isvector(magnitude_spectrum)
    error('calculate_thd:InvalidInput', 'Magnitude spectrum must be a numeric vector');
end

% Ensure magnitude spectrum is a column vector
magnitude_spectrum = magnitude_spectrum(:);

% Set default method
if nargin < 2 || isempty(method)
    method = 'fundamental';  % IEEE standard method
end

% Validate method
valid_methods = {'fundamental', 'rms'};
if ~ismember(lower(method), valid_methods)
    error('calculate_thd:InvalidMethod', ...
        'Method must be ''fundamental'' or ''rms''');
end

method = lower(method);

% Check minimum length
if length(magnitude_spectrum) < 2
    error('calculate_thd:InsufficientData', ...
        'Magnitude spectrum must have at least 2 elements (DC and fundamental)');
end

% =======================================================================
% EXTRACT FUNDAMENTAL AND HARMONICS
% =======================================================================

% DC component (index 1, k=0)
V_dc = magnitude_spectrum(1);

% Fundamental component (index 2, k=1)
V_fundamental_peak = magnitude_spectrum(2);
V_fundamental_rms = V_fundamental_peak / sqrt(2);

% Check if fundamental is non-zero
if V_fundamental_peak < 1e-10
    warning('calculate_thd:ZeroFundamental', ...
        'Fundamental component is nearly zero (%.2e V). THD calculation may be inaccurate.', ...
        V_fundamental_peak);

    % Return NaN if fundamental is essentially zero
    THD_percent = NaN;
    THD_ratio = NaN;
    harmonic_details = struct();
    return;
end

% Harmonic components (indices 3 onwards, k=2, 3, 4, ...)
% We analyze up to 50th harmonic or end of spectrum, whichever is smaller
max_harmonic_index = min(52, length(magnitude_spectrum));  % Index 52 = 50th harmonic
num_harmonics = max_harmonic_index - 2;  % Exclude DC and fundamental

if num_harmonics < 1
    warning('calculate_thd:NoHarmonics', ...
        'No harmonic components found. THD = 0%%');
    THD_percent = 0;
    THD_ratio = 0;
    harmonic_details = struct();
    harmonic_details.fundamental_mag = V_fundamental_peak;
    harmonic_details.fundamental_rms = V_fundamental_rms;
    return;
end

% Extract harmonic magnitudes (peak values)
harmonic_magnitudes_peak = magnitude_spectrum(3:max_harmonic_index);
harmonic_magnitudes_rms = harmonic_magnitudes_peak / sqrt(2);

% =======================================================================
% CALCULATE THD
% =======================================================================

% Sum of squares of all harmonic RMS values
harmonic_sum_of_squares = sum(harmonic_magnitudes_rms.^2);
harmonic_rms_total = sqrt(harmonic_sum_of_squares);

% Calculate total signal RMS
% V_total_rms = sqrt(V_dc^2 + V_fundamental_rms^2 + sum(V_harmonic_rms^2))
V_total_rms = sqrt(V_dc^2 + V_fundamental_rms^2 + harmonic_sum_of_squares);

% Calculate THD based on selected method
switch method
    case 'fundamental'
        % THD_F = sqrt(sum(Vh_rms^2)) / V1_rms * 100%
        % This is the IEEE 519 standard definition
        THD_ratio = harmonic_rms_total / V_fundamental_rms;
        THD_percent = THD_ratio * 100;
        thd_method_name = 'THD-F (Fundamental Reference)';

    case 'rms'
        % THD_R = sqrt(sum(Vh_rms^2)) / Vtotal_rms * 100%
        % Alternative definition, less common
        THD_ratio = harmonic_rms_total / V_total_rms;
        THD_percent = THD_ratio * 100;
        thd_method_name = 'THD-R (RMS Reference)';
end

% =======================================================================
% CALCULATE INDIVIDUAL HARMONIC DISTORTION (IHD)
% =======================================================================

% Individual Harmonic Distortion for each harmonic
% IHD_n = V_n_rms / V1_rms * 100%
individual_thd = (harmonic_magnitudes_rms / V_fundamental_rms) * 100;

% Harmonic percentages (same as IHD)
harmonic_percentages = individual_thd;

% =======================================================================
% FIND DOMINANT HARMONIC
% =======================================================================

% Find harmonic with maximum magnitude
[max_harmonic_mag, max_harmonic_idx] = max(harmonic_magnitudes_peak);
dominant_harmonic_number = max_harmonic_idx + 1;  % +1 because harmonics start from k=2

% =======================================================================
% ASSESS POWER QUALITY (IEEE 519 STANDARD)
% =======================================================================

if THD_percent < 5.0
    quality_rating = 'EXCELLENT';
    quality_color = 'green';
elseif THD_percent < 8.0
    quality_rating = 'GOOD';
    quality_color = 'yellow';
else
    quality_rating = 'POOR (Action Required)';
    quality_color = 'red';
end

% =======================================================================
% POPULATE HARMONIC DETAILS STRUCTURE
% =======================================================================

harmonic_details.fundamental_mag = V_fundamental_peak;
harmonic_details.fundamental_rms = V_fundamental_rms;
harmonic_details.total_rms = V_total_rms;
harmonic_details.harmonic_rms = harmonic_rms_total;
harmonic_details.individual_thd = individual_thd;
harmonic_details.harmonic_magnitudes = harmonic_magnitudes_peak;
harmonic_details.harmonic_percentages = harmonic_percentages;
harmonic_details.dominant_harmonic = dominant_harmonic_number;
harmonic_details.dominant_harmonic_magnitude = max_harmonic_mag;
harmonic_details.num_harmonics = num_harmonics;
harmonic_details.method = thd_method_name;
harmonic_details.quality_rating = quality_rating;
harmonic_details.dc_component = V_dc;

% =======================================================================
% DISPLAY RESULTS
% =======================================================================

fprintf('\n========================================\n');
fprintf('TOTAL HARMONIC DISTORTION ANALYSIS\n');
fprintf('========================================\n');
fprintf('Method: %s\n', thd_method_name);
fprintf('----------------------------------------\n');
fprintf('Fundamental (50 Hz):\n');
fprintf('  Peak Voltage: %.4f V\n', V_fundamental_peak);
fprintf('  RMS Voltage:  %.4f V\n', V_fundamental_rms);
fprintf('----------------------------------------\n');
fprintf('Harmonic Content:\n');
fprintf('  Total Harmonic RMS: %.4f V\n', harmonic_rms_total);
fprintf('  Number of Harmonics: %d (2nd to %dth)\n', num_harmonics, num_harmonics+1);
fprintf('  Dominant Harmonic: %dth (%.4f V peak)\n', ...
    dominant_harmonic_number, max_harmonic_mag);
fprintf('----------------------------------------\n');
fprintf('THD Results:\n');
fprintf('  THD: %.4f%% (%.6f ratio)\n', THD_percent, THD_ratio);
fprintf('  Quality Rating: %s\n', quality_rating);
fprintf('----------------------------------------\n');

% Display IEEE standards reference
fprintf('IEEE 519 Standards:\n');
fprintf('  THD < 5%%  : Excellent\n');
fprintf('  THD 5-8%% : Good\n');
fprintf('  THD > 8%%  : Poor (requires correction)\n');
fprintf('========================================\n');

% Display top 5 harmonics
fprintf('\nTop 5 Harmonics by Magnitude:\n');
fprintf('  Order | Frequency | Magnitude | IHD (%%) \n');
fprintf('  ------|-----------|-----------|----------\n');

% Sort harmonics by magnitude
[sorted_mags, sorted_indices] = sort(harmonic_magnitudes_peak, 'descend');
num_to_display = min(5, length(sorted_mags));

for i = 1:num_to_display
    harmonic_order = sorted_indices(i) + 1;  % +1 because harmonics start from 2nd
    harmonic_freq = harmonic_order * 50;     % Assuming 50 Hz fundamental
    mag = sorted_mags(i);
    ihd = individual_thd(sorted_indices(i));

    fprintf('   %2d   | %6d Hz | %7.4f V | %7.3f%%\n', ...
        harmonic_order, harmonic_freq, mag, ihd);
end

fprintf('========================================\n\n');

% =======================================================================
% END OF FUNCTION
% =======================================================================

end
