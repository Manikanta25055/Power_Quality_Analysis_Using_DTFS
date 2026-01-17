function [metrics] = power_quality_metrics(signal, magnitude_spectrum, fs, f0)
% POWER_QUALITY_METRICS - Comprehensive power quality analysis
%
% DESCRIPTION:
%   Calculates all important power quality metrics according to IEEE standards.
%   This function provides a complete assessment of electrical power quality
%   including distortion, waveform characteristics, and equipment stress factors.
%
% SYNTAX:
%   metrics = power_quality_metrics(signal, magnitude_spectrum, fs, f0)
%
% INPUTS:
%   signal             - Time-domain voltage/current signal [N x 1 vector]
%   magnitude_spectrum - DTFS magnitude spectrum from calculate_dtfs [N x 1 vector]
%   fs                 - Sampling frequency (Hz)
%   f0                 - Fundamental frequency (Hz), typically 50 or 60 Hz
%
% OUTPUTS:
%   metrics - Structure containing all power quality metrics:
%
%   TIME DOMAIN METRICS:
%   .signal_rms        - RMS value of signal (V or A)
%   .signal_peak       - Peak value (V or A)
%   .signal_avg        - Average (DC) value (V or A)
%   .signal_pp         - Peak-to-peak value (V or A)
%
%   WAVEFORM SHAPE METRICS:
%   .crest_factor      - Peak/RMS ratio (dimensionless)
%   .form_factor       - RMS/Average ratio (dimensionless)
%   .distortion_factor - Ratio of harmonic content to total
%
%   HARMONIC DISTORTION METRICS:
%   .THD_F             - Total Harmonic Distortion (Fundamental reference) %
%   .THD_R             - Total Harmonic Distortion (RMS reference) %
%   .fundamental_mag   - Fundamental component magnitude (V or A)
%   .fundamental_rms   - Fundamental component RMS (V or A)
%   .harmonic_rms      - RMS of all harmonics (V or A)
%
%   POWER SYSTEM METRICS:
%   .displacement_PF   - Displacement Power Factor (fundamental only)
%   .distortion_PF     - Distortion Power Factor (due to harmonics)
%   .true_PF           - True Power Factor (combined)
%
%   EQUIPMENT STRESS METRICS:
%   .K_factor          - Transformer K-factor (derating factor)
%   .TIF               - Telephone Influence Factor
%   .IT_product        - Current-Time product for fuses
%
%   FREQUENCY DOMAIN METRICS:
%   .harmonic_content  - Table of individual harmonics with:
%                        [order, frequency, magnitude, percentage, IHD]
%   .dominant_harmonic - Order of largest harmonic
%   .even_harmonic_distortion - THD from even harmonics only %
%   .odd_harmonic_distortion  - THD from odd harmonics only %
%
%   QUALITY ASSESSMENT:
%   .quality_rating    - IEEE 519 rating: EXCELLENT/GOOD/POOR
%   .compliance_IEEE519- Boolean: true if THD < 5%
%   .recommended_action- String with recommendations
%
% IEEE STANDARDS REFERENCE:
%   - IEEE 519-2014: Harmonic Control in Electrical Power Systems
%   - IEEE 1459-2010: Definitions for Measurement of Electric Power Quantities
%   - IEC 61000-4-7: Harmonics and interharmonics measurement
%
% EXAMPLE USAGE:
%   % Generate and analyze a signal
%   [signal, t, info] = generate_distorted_signal('motor_drive', 10000, 50, 200, 230);
%   [X_k, mag, phase, freq] = calculate_dtfs(signal, 10000, 200);
%   metrics = power_quality_metrics(signal, mag, 10000, 50);
%
%   fprintf('THD: %.2f%%\n', metrics.THD_F);
%   fprintf('Crest Factor: %.3f\n', metrics.crest_factor);
%   fprintf('K-Factor: %.2f\n', metrics.K_factor);
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025
% COURSE: Digital Signal Processing (FISAC Assessment)
% PROJECT: Power Quality Analysis using DTFS

% =======================================================================
% INPUT VALIDATION
% =======================================================================

if nargin < 4
    error('power_quality_metrics:InsufficientInputs', ...
        'All inputs required: signal, magnitude_spectrum, fs, f0');
end

% Ensure inputs are column vectors
signal = signal(:);
magnitude_spectrum = magnitude_spectrum(:);

% Validate signal
if ~isnumeric(signal) || length(signal) < 2
    error('power_quality_metrics:InvalidSignal', 'Signal must be a numeric vector');
end

% Validate magnitude spectrum
if ~isnumeric(magnitude_spectrum) || length(magnitude_spectrum) < 2
    error('power_quality_metrics:InvalidSpectrum', 'Magnitude spectrum must be a numeric vector');
end

% =======================================================================
% EXTRACT FUNDAMENTAL AND HARMONIC COMPONENTS
% =======================================================================

% DC component (k=0, index 1)
V_dc = magnitude_spectrum(1);

% Fundamental component (k=1, index 2)
V1_peak = magnitude_spectrum(2);
V1_rms = V1_peak / sqrt(2);

% All harmonics (k>=2, indices 3+)
max_harmonic_idx = min(52, length(magnitude_spectrum));  % Up to 50th harmonic
num_harmonics = max_harmonic_idx - 2;

if num_harmonics < 1
    harmonic_mags_peak = [];
    harmonic_mags_rms = [];
else
    harmonic_mags_peak = magnitude_spectrum(3:max_harmonic_idx);
    harmonic_mags_rms = harmonic_mags_peak / sqrt(2);
end

% =======================================================================
% TIME DOMAIN METRICS
% =======================================================================

% RMS value (Root Mean Square)
signal_rms = sqrt(mean(signal.^2));

% Peak value
signal_peak = max(abs(signal));

% Average (DC) value
signal_avg = mean(abs(signal));  % Rectified average for AC signals

% Peak-to-peak value
signal_pp = max(signal) - min(signal);

% =======================================================================
% WAVEFORM SHAPE METRICS
% =======================================================================

% Crest Factor (Peak/RMS)
% Ideal sine wave: 1.414 (sqrt(2))
% Higher values indicate sharp peaks (stress on insulation)
crest_factor = signal_peak / signal_rms;

% Form Factor (RMS/Average)
% Ideal sine wave: 1.111 (pi/(2*sqrt(2)))
% Indicates deviation from ideal waveform
if signal_avg > 1e-10
    form_factor = signal_rms / signal_avg;
else
    form_factor = NaN;
end

% Distortion Factor
% Ratio of harmonic content to total signal
if ~isempty(harmonic_mags_rms)
    harmonic_rms_total = sqrt(sum(harmonic_mags_rms.^2));
    distortion_factor = harmonic_rms_total / signal_rms;
else
    harmonic_rms_total = 0;
    distortion_factor = 0;
end

% =======================================================================
% HARMONIC DISTORTION METRICS
% =======================================================================

% Calculate THD using both methods
if V1_rms > 1e-10 && ~isempty(harmonic_mags_rms)
    % THD-F (Fundamental reference) - IEEE standard
    THD_F = (harmonic_rms_total / V1_rms) * 100;

    % THD-R (RMS reference) - Alternative
    THD_R = (harmonic_rms_total / signal_rms) * 100;
else
    THD_F = 0;
    THD_R = 0;
end

% =======================================================================
% SEPARATE EVEN AND ODD HARMONIC DISTORTION
% =======================================================================

if ~isempty(harmonic_mags_rms)
    % Even harmonics: 2nd, 4th, 6th, ... (indices 1, 3, 5, ...)
    even_indices = 1:2:length(harmonic_mags_rms);
    even_harmonic_rms = sqrt(sum(harmonic_mags_rms(even_indices).^2));
    even_THD = (even_harmonic_rms / V1_rms) * 100;

    % Odd harmonics: 3rd, 5th, 7th, ... (indices 2, 4, 6, ...)
    odd_indices = 2:2:length(harmonic_mags_rms);
    if ~isempty(odd_indices)
        odd_harmonic_rms = sqrt(sum(harmonic_mags_rms(odd_indices).^2));
        odd_THD = (odd_harmonic_rms / V1_rms) * 100;
    else
        odd_harmonic_rms = 0;
        odd_THD = 0;
    end
else
    even_THD = 0;
    odd_THD = 0;
end

% =======================================================================
% POWER FACTOR CALCULATIONS
% =======================================================================

% Displacement Power Factor (DPF)
% For fundamental frequency only (assumes resistive load at fundamental)
% In single-phase systems with harmonics but no phase shift at fundamental
displacement_PF = 1.0;  % Assumes in-phase fundamental (conservative)

% Distortion Power Factor
% Accounts for harmonic distortion
% DPF_distortion = 1 / sqrt(1 + THD^2)
if THD_F > 0
    distortion_PF = 1 / sqrt(1 + (THD_F/100)^2);
else
    distortion_PF = 1.0;
end

% True Power Factor (TPF)
% Combination of displacement and distortion
% TPF = DPF * Distortion_PF
true_PF = displacement_PF * distortion_PF;

% =======================================================================
% TRANSFORMER K-FACTOR
% =======================================================================

% K-Factor indicates additional heating in transformers due to harmonics
% K = sum(h^2 * I_h^2) / sum(I_h^2)
% where h is harmonic order, I_h is harmonic current
%
% For voltage analysis, we approximate using voltage harmonics
% Typical K-factors:
%   K-1: Linear loads (no harmonics)
%   K-4: Computers, office equipment
%   K-13: Data centers, heavy non-linear loads
%   K-20: Extreme harmonic environments

if ~isempty(harmonic_mags_rms) && V1_rms > 1e-10
    K_factor = 0;
    total_harmonic_power = 0;

    for h = 1:length(harmonic_mags_rms)
        harmonic_order = h + 1;  % 2nd, 3rd, 4th, ... harmonics
        I_h_normalized = harmonic_mags_rms(h) / V1_rms;  % Normalized to fundamental
        K_factor = K_factor + (harmonic_order^2) * (I_h_normalized^2);
        total_harmonic_power = total_harmonic_power + I_h_normalized^2;
    end

    % Normalize and add fundamental contribution
    K_factor = K_factor + 1.0;  % Fundamental (h=1, contributes 1.0)

else
    K_factor = 1.0;  % No harmonics = K-1 rating
end

% =======================================================================
% TELEPHONE INFLUENCE FACTOR (TIF)
% =======================================================================

% TIF measures interference with telephone systems
% TIF = sqrt(sum((h * V_h)^2)) / V1
% Higher frequencies have more influence on telephone circuits

if ~isempty(harmonic_mags_peak) && V1_peak > 1e-10
    TIF_sum = 0;

    for h = 1:min(15, length(harmonic_mags_peak))  % Usually calculated up to 15th
        harmonic_order = h + 1;
        % Weighting factor increases with frequency
        weight = harmonic_order;
        TIF_sum = TIF_sum + (weight * harmonic_mags_peak(h))^2;
    end

    TIF = sqrt(TIF_sum) / V1_peak;
else
    TIF = 0;
end

% =======================================================================
% I-T PRODUCT (FOR FUSE/BREAKER COORDINATION)
% =======================================================================

% Current-Time product affects fuse and breaker operation
% I^2*t product determines thermal stress
% For voltage signals, this is informational

IT_product = signal_rms^2;  % Normalized to 1 second

% =======================================================================
% HARMONIC CONTENT TABLE
% =======================================================================

if ~isempty(harmonic_mags_peak)
    % Create detailed harmonic table
    % Columns: [Order, Frequency(Hz), Magnitude, Percentage, IHD(%)]

    num_harmonics_table = min(20, length(harmonic_mags_peak));  % Show up to 20th
    harmonic_content = zeros(num_harmonics_table, 5);

    for h = 1:num_harmonics_table
        harmonic_order = h + 1;  % 2nd, 3rd, 4th, ...
        harmonic_freq = harmonic_order * f0;
        harmonic_mag = harmonic_mags_peak(h);
        harmonic_pct = (harmonic_mag / V1_peak) * 100;
        harmonic_IHD = (harmonic_mags_rms(h) / V1_rms) * 100;

        harmonic_content(h, :) = [harmonic_order, harmonic_freq, harmonic_mag, ...
                                   harmonic_pct, harmonic_IHD];
    end

    % Find dominant harmonic
    [~, max_idx] = max(harmonic_mags_peak);
    dominant_harmonic = max_idx + 1;  % Harmonic order

else
    harmonic_content = [];
    dominant_harmonic = 1;  % Only fundamental present
end

% =======================================================================
% QUALITY ASSESSMENT (IEEE 519)
% =======================================================================

if THD_F < 5.0
    quality_rating = 'EXCELLENT';
    compliance_IEEE519 = true;
    recommended_action = 'No action required. Power quality is excellent.';
elseif THD_F < 8.0
    quality_rating = 'GOOD';
    compliance_IEEE519 = false;
    recommended_action = 'Monitor power quality. Consider filtering if sensitive equipment present.';
else
    quality_rating = 'POOR (Action Required)';
    compliance_IEEE519 = false;

    % Provide specific recommendations based on dominant harmonics
    if dominant_harmonic == 3
        recommended_action = 'High 3rd harmonic detected. Install harmonic filter or use delta-wye transformers.';
    elseif dominant_harmonic == 5 || dominant_harmonic == 7
        recommended_action = 'Characteristic VFD harmonics detected. Install 6-pulse or 12-pulse filter.';
    else
        recommended_action = 'Multiple harmonics present. Conduct detailed analysis and install active harmonic filter.';
    end
end

% =======================================================================
% POPULATE OUTPUT STRUCTURE
% =======================================================================

% Time domain metrics
metrics.signal_rms = signal_rms;
metrics.signal_peak = signal_peak;
metrics.signal_avg = signal_avg;
metrics.signal_pp = signal_pp;

% Waveform shape metrics
metrics.crest_factor = crest_factor;
metrics.form_factor = form_factor;
metrics.distortion_factor = distortion_factor;

% Harmonic distortion metrics
metrics.THD_F = THD_F;
metrics.THD_R = THD_R;
metrics.fundamental_mag = V1_peak;
metrics.fundamental_rms = V1_rms;
metrics.harmonic_rms = harmonic_rms_total;

% Separate even/odd harmonics
metrics.even_harmonic_distortion = even_THD;
metrics.odd_harmonic_distortion = odd_THD;

% Power factor metrics
metrics.displacement_PF = displacement_PF;
metrics.distortion_PF = distortion_PF;
metrics.true_PF = true_PF;

% Equipment stress metrics
metrics.K_factor = K_factor;
metrics.TIF = TIF;
metrics.IT_product = IT_product;

% Frequency domain
metrics.harmonic_content = harmonic_content;
metrics.dominant_harmonic = dominant_harmonic;

% Quality assessment
metrics.quality_rating = quality_rating;
metrics.compliance_IEEE519 = compliance_IEEE519;
metrics.recommended_action = recommended_action;

% Additional info
metrics.fundamental_frequency = f0;
metrics.sampling_frequency = fs;

% =======================================================================
% DISPLAY SUMMARY
% =======================================================================

fprintf('\n========================================\n');
fprintf('COMPREHENSIVE POWER QUALITY METRICS\n');
fprintf('========================================\n\n');

fprintf('TIME DOMAIN METRICS:\n');
fprintf('  RMS Value:        %8.4f V/A\n', signal_rms);
fprintf('  Peak Value:       %8.4f V/A\n', signal_peak);
fprintf('  Average Value:    %8.4f V/A\n', signal_avg);
fprintf('  Peak-to-Peak:     %8.4f V/A\n', signal_pp);
fprintf('\n');

fprintf('WAVEFORM CHARACTERISTICS:\n');
fprintf('  Crest Factor:     %8.4f (ideal: 1.414)\n', crest_factor);
fprintf('  Form Factor:      %8.4f (ideal: 1.111)\n', form_factor);
fprintf('  Distortion Factor:%8.4f\n', distortion_factor);
fprintf('\n');

fprintf('HARMONIC DISTORTION:\n');
fprintf('  THD-F:            %8.4f %% (IEEE 519 standard)\n', THD_F);
fprintf('  THD-R:            %8.4f %%\n', THD_R);
fprintf('  Even Harmonics:   %8.4f %%\n', even_THD);
fprintf('  Odd Harmonics:    %8.4f %%\n', odd_THD);
fprintf('  Dominant:         %dth harmonic (%.0f Hz)\n', ...
    dominant_harmonic, dominant_harmonic * f0);
fprintf('\n');

fprintf('POWER SYSTEM METRICS:\n');
fprintf('  Displacement PF:  %8.4f\n', displacement_PF);
fprintf('  Distortion PF:    %8.4f\n', distortion_PF);
fprintf('  True PF:          %8.4f\n', true_PF);
fprintf('\n');

fprintf('EQUIPMENT STRESS FACTORS:\n');
fprintf('  K-Factor:         %8.2f', K_factor);
if K_factor < 4
    fprintf(' (Low - Linear load)\n');
elseif K_factor < 13
    fprintf(' (Medium - Office equipment)\n');
elseif K_factor < 20
    fprintf(' (High - Data center)\n');
else
    fprintf(' (Very High - Extreme)\n');
end
fprintf('  TIF:              %8.4f\n', TIF);
fprintf('  I²t Product:      %8.4f\n', IT_product);
fprintf('\n');

fprintf('QUALITY ASSESSMENT:\n');
fprintf('  IEEE 519 Rating:  %s\n', quality_rating);
fprintf('  Compliance:       %s\n', mat2str(compliance_IEEE519));
fprintf('  Recommendation:   %s\n', recommended_action);
fprintf('========================================\n\n');

% =======================================================================
% END OF FUNCTION
% =======================================================================

end
