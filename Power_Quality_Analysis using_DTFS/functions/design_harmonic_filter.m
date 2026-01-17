function [X_k_filtered, filter_info] = design_harmonic_filter(X_k, N, filter_type, harmonics_to_remove, attenuation_dB)
% DESIGN_HARMONIC_FILTER - Designs and applies harmonic filters in DTFS domain
%
% DESCRIPTION:
%   Removes or attenuates specific harmonics from power signals using
%   DTFS domain filtering. Maintains conjugate symmetry to ensure
%   real-valued output signals. This represents active harmonic filtering
%   used in industrial power quality improvement systems.
%
% SYNTAX:
%   [X_k_filtered, filter_info] = design_harmonic_filter(X_k, N, filter_type, harmonics_to_remove)
%   [X_k_filtered, filter_info] = design_harmonic_filter(X_k, N, filter_type, harmonics_to_remove, attenuation_dB)
%
% INPUTS:
%   X_k                 - DTFS coefficients from calculate_dtfs [N x 1 complex vector]
%   N                   - Period length (samples per fundamental cycle)
%   filter_type         - Filter strategy:
%                         'notch'    - Remove specific harmonics completely
%                         'attenuate'- Reduce specific harmonics by attenuation_dB
%                         'lowpass'  - Keep fundamental + harmonics up to order M
%                         'highpass' - Remove fundamental, keep only harmonics
%                         'bandpass' - Keep only harmonics in specified range
%
%   harmonics_to_remove - Vector of harmonic orders to filter
%                         For 'notch'/'attenuate': [3 5 7] removes 3rd, 5th, 7th
%                         For 'lowpass': scalar M (keep up to Mth harmonic)
%                         For 'bandpass': [M1 M2] (keep M1 to M2)
%
%   attenuation_dB      - (Optional) Attenuation in dB for 'attenuate' mode
%                         Default: Inf (complete removal for 'notch')
%                         Example: 40 dB reduces harmonic to 1% of original
%
% OUTPUTS:
%   X_k_filtered - Filtered DTFS coefficients [N x 1 complex vector]
%   filter_info  - Structure containing filter information:
%                  .filter_type        - Type of filter applied
%                  .harmonics_removed  - List of harmonic orders affected
%                  .attenuation_dB     - Attenuation applied
%                  .original_THD       - THD before filtering
%                  .filtered_THD       - THD after filtering
%                  .THD_reduction      - Percentage reduction in THD
%                  .harmonics_before   - Magnitudes before filtering
%                  .harmonics_after    - Magnitudes after filtering
%                  .filter_effectiveness - Percentage of harmonic energy removed
%
% CONJUGATE SYMMETRY:
%   For real-valued signals, DTFS coefficients satisfy:
%   X[k] = conj(X[N-k])
%
%   When filtering harmonic h, we must remove both:
%   - Positive frequency: X[h] at index h+1
%   - Negative frequency: X[N-h] at index N-h+1
%
% EXAMPLE USAGE:
%   % Remove 3rd and 5th harmonics (LED lighting mitigation)
%   [X_k, mag, ~, ~] = calculate_dtfs(signal, fs, N);
%   [X_filtered, info] = design_harmonic_filter(X_k, N, 'notch', [3 5]);
%   [signal_clean, t] = synthesize_from_dtfs(X_filtered, fs, N, 3);
%
%   % Attenuate 5th and 7th by 40 dB (VFD mitigation)
%   [X_filtered, info] = design_harmonic_filter(X_k, N, 'attenuate', [5 7], 40);
%
%   % Lowpass: Keep only fundamental + 3rd harmonic
%   [X_filtered, info] = design_harmonic_filter(X_k, N, 'lowpass', 3);
%
% INDUSTRIAL APPLICATIONS:
%   - Active Harmonic Filters (AHF): Inject compensation currents
%   - Passive LC Filters: Tuned to specific harmonic frequencies
%   - Delta-Wye Transformers: Cancel triplen (3rd, 9th, 15th) harmonics
%   - 12-Pulse Rectifiers: Eliminate 5th and 7th harmonics
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025
% COURSE: Digital Signal Processing (FISAC Assessment)
% PROJECT: Power Quality Analysis using DTFS

% =======================================================================
% INPUT VALIDATION
% =======================================================================

if nargin < 4
    error('design_harmonic_filter:InsufficientInputs', ...
        'Required: X_k, N, filter_type, harmonics_to_remove');
end

% Set default attenuation
if nargin < 5 || isempty(attenuation_dB)
    if strcmp(filter_type, 'attenuate')
        attenuation_dB = 40;  % Default 40 dB attenuation
    else
        attenuation_dB = Inf;  % Complete removal for other modes
    end
end

% Validate X_k
if ~isnumeric(X_k) || ~isvector(X_k)
    error('design_harmonic_filter:InvalidX_k', 'X_k must be a complex vector');
end

X_k = X_k(:);  % Ensure column vector

% Validate N
if length(X_k) ~= N
    warning('design_harmonic_filter:LengthMismatch', ...
        'X_k length (%d) does not match N (%d). Using length(X_k).', length(X_k), N);
    N = length(X_k);
end

% Validate filter type
valid_types = {'notch', 'attenuate', 'lowpass', 'highpass', 'bandpass'};
if ~ismember(lower(filter_type), valid_types)
    error('design_harmonic_filter:InvalidType', ...
        'Filter type must be: notch, attenuate, lowpass, highpass, or bandpass');
end

filter_type = lower(filter_type);

% =======================================================================
% CALCULATE ORIGINAL SPECTRUM AND THD
% =======================================================================

% Extract magnitudes
mag_original = abs(X_k);
mag_original(2:end) = 2 * mag_original(2:end);  % Scale for peak values

% Calculate original THD
V1_rms = mag_original(2) / sqrt(2);
if V1_rms > 1e-10
    harmonic_rms_original = sqrt(sum((mag_original(3:end) / sqrt(2)).^2));
    THD_original = (harmonic_rms_original / V1_rms) * 100;
else
    THD_original = 0;
end

% =======================================================================
% APPLY FILTER BASED ON TYPE
% =======================================================================

% Create a copy for filtering
X_k_filtered = X_k;

% Calculate attenuation factor from dB
% dB = 20*log10(factor) => factor = 10^(-dB/20)
attenuation_factor = 10^(-attenuation_dB/20);

% Track which harmonics are actually removed
harmonics_affected = [];

fprintf('\n========================================\n');
fprintf('HARMONIC FILTER DESIGN\n');
fprintf('========================================\n');
fprintf('Filter Type: %s\n', upper(filter_type));
fprintf('Period Length: %d samples\n', N);

switch filter_type

    case 'notch'
        % Remove specific harmonics completely
        fprintf('Harmonics to Remove: ');
        fprintf('%d ', harmonics_to_remove);
        fprintf('\n');
        fprintf('Attenuation: Complete removal (Inf dB)\n');
        fprintf('----------------------------------------\n');

        for i = 1:length(harmonics_to_remove)
            h = harmonics_to_remove(i);

            if h < 1 || h >= N/2
                warning('Harmonic %d out of valid range [1, %d). Skipping.', h, floor(N/2));
                continue;
            end

            % Positive frequency component: k=h, index=h+1
            idx_pos = h + 1;

            % Negative frequency component (conjugate): k=N-h, index=N-h+1
            idx_neg = N - h + 1;

            % Store original magnitude
            mag_before = mag_original(idx_pos);

            % Zero out both components (maintain conjugate symmetry)
            X_k_filtered(idx_pos) = 0;
            X_k_filtered(idx_neg) = 0;

            harmonics_affected = [harmonics_affected; h];

            fprintf('  Removed %dth harmonic (%d Hz): %.4f V -> 0 V\n', ...
                h, h*50, mag_before);
        end

    case 'attenuate'
        % Attenuate specific harmonics by specified dB
        fprintf('Harmonics to Attenuate: ');
        fprintf('%d ', harmonics_to_remove);
        fprintf('\n');
        fprintf('Attenuation: %.1f dB (factor: %.6f)\n', attenuation_dB, attenuation_factor);
        fprintf('----------------------------------------\n');

        for i = 1:length(harmonics_to_remove)
            h = harmonics_to_remove(i);

            if h < 1 || h >= N/2
                warning('Harmonic %d out of valid range [1, %d). Skipping.', h, floor(N/2));
                continue;
            end

            idx_pos = h + 1;
            idx_neg = N - h + 1;

            mag_before = mag_original(idx_pos);

            % Attenuate both components
            X_k_filtered(idx_pos) = X_k_filtered(idx_pos) * attenuation_factor;
            X_k_filtered(idx_neg) = X_k_filtered(idx_neg) * attenuation_factor;

            harmonics_affected = [harmonics_affected; h];

            mag_after = mag_before * attenuation_factor;
            fprintf('  Attenuated %dth harmonic (%d Hz): %.4f V -> %.4f V (%.1f dB)\n', ...
                h, h*50, mag_before, mag_after, attenuation_dB);
        end

    case 'lowpass'
        % Keep fundamental + harmonics up to order M
        M = harmonics_to_remove(1);  % Max harmonic order to keep

        fprintf('Lowpass Order: %d (keep up to %dth harmonic)\n', M, M);
        fprintf('Fundamental Frequency: 50 Hz\n');
        fprintf('Cutoff Frequency: %d Hz\n', M*50);
        fprintf('----------------------------------------\n');

        % Zero out all harmonics above M
        for h = M+1:floor(N/2)-1
            idx_pos = h + 1;
            idx_neg = N - h + 1;

            if mag_original(idx_pos) > 1e-6  % Only report significant harmonics
                fprintf('  Removed %dth harmonic (%d Hz): %.4f V -> 0 V\n', ...
                    h, h*50, mag_original(idx_pos));
                harmonics_affected = [harmonics_affected; h];
            end

            X_k_filtered(idx_pos) = 0;
            X_k_filtered(idx_neg) = 0;
        end

    case 'highpass'
        % Remove fundamental, keep only harmonics
        fprintf('Highpass Filter: Remove fundamental (50 Hz)\n');
        fprintf('Keep all harmonics (>= 100 Hz)\n');
        fprintf('----------------------------------------\n');

        % Zero out fundamental (k=1)
        mag_fund = mag_original(2);
        X_k_filtered(2) = 0;           % Positive frequency
        X_k_filtered(N) = 0;           % Negative frequency (k=N-1)

        fprintf('  Removed fundamental (50 Hz): %.4f V -> 0 V\n', mag_fund);
        harmonics_affected = [1];

    case 'bandpass'
        % Keep only harmonics in range [M1, M2]
        M1 = harmonics_to_remove(1);
        M2 = harmonics_to_remove(2);

        fprintf('Bandpass Filter: Keep harmonics %d to %d\n', M1, M2);
        fprintf('Passband: %d Hz to %d Hz\n', M1*50, M2*50);
        fprintf('----------------------------------------\n');

        % Remove fundamental if M1 > 1
        if M1 > 1
            mag_fund = mag_original(2);
            X_k_filtered(2) = 0;
            X_k_filtered(N) = 0;
            fprintf('  Removed fundamental (50 Hz): %.4f V -> 0 V\n', mag_fund);
            harmonics_affected = [1];
        end

        % Remove harmonics below M1
        for h = 2:M1-1
            idx_pos = h + 1;
            idx_neg = N - h + 1;

            if mag_original(idx_pos) > 1e-6
                fprintf('  Removed %dth harmonic (%d Hz): %.4f V -> 0 V\n', ...
                    h, h*50, mag_original(idx_pos));
                harmonics_affected = [harmonics_affected; h];
            end

            X_k_filtered(idx_pos) = 0;
            X_k_filtered(idx_neg) = 0;
        end

        % Remove harmonics above M2
        for h = M2+1:floor(N/2)-1
            idx_pos = h + 1;
            idx_neg = N - h + 1;

            if mag_original(idx_pos) > 1e-6
                fprintf('  Removed %dth harmonic (%d Hz): %.4f V -> 0 V\n', ...
                    h, h*50, mag_original(idx_pos));
                harmonics_affected = [harmonics_affected; h];
            end

            X_k_filtered(idx_pos) = 0;
            X_k_filtered(idx_neg) = 0;
        end
end

% =======================================================================
% CALCULATE FILTERED SPECTRUM AND THD
% =======================================================================

% Extract filtered magnitudes
mag_filtered = abs(X_k_filtered);
mag_filtered(2:end) = 2 * mag_filtered(2:end);

% Calculate filtered THD
V1_filtered_rms = mag_filtered(2) / sqrt(2);
if V1_filtered_rms > 1e-10
    harmonic_rms_filtered = sqrt(sum((mag_filtered(3:end) / sqrt(2)).^2));
    THD_filtered = (harmonic_rms_filtered / V1_filtered_rms) * 100;
else
    THD_filtered = 0;
end

% Calculate THD reduction
if THD_original > 1e-6
    THD_reduction = ((THD_original - THD_filtered) / THD_original) * 100;
else
    THD_reduction = 0;
end

% Calculate filter effectiveness (harmonic energy removed)
energy_original = sum(mag_original(3:end).^2);
energy_filtered = sum(mag_filtered(3:end).^2);

if energy_original > 1e-10
    filter_effectiveness = ((energy_original - energy_filtered) / energy_original) * 100;
else
    filter_effectiveness = 0;
end

% =======================================================================
% POPULATE OUTPUT STRUCTURE
% =======================================================================

filter_info.filter_type = filter_type;
filter_info.harmonics_affected = harmonics_affected;
filter_info.attenuation_dB = attenuation_dB;
filter_info.attenuation_factor = attenuation_factor;
filter_info.original_THD = THD_original;
filter_info.filtered_THD = THD_filtered;
filter_info.THD_reduction_percent = THD_reduction;
filter_info.harmonics_before = mag_original;
filter_info.harmonics_after = mag_filtered;
filter_info.filter_effectiveness = filter_effectiveness;
filter_info.fundamental_preserved = mag_filtered(2) > 1e-6;

% =======================================================================
% DISPLAY RESULTS
% =======================================================================

fprintf('========================================\n');
fprintf('FILTER PERFORMANCE\n');
fprintf('========================================\n');
fprintf('Original THD:    %.4f %%\n', THD_original);
fprintf('Filtered THD:    %.4f %%\n', THD_filtered);
fprintf('THD Reduction:   %.2f %%\n', THD_reduction);
fprintf('Filter Effectiveness: %.2f %% (harmonic energy removed)\n', filter_effectiveness);
fprintf('Fundamental:     %.4f V (preserved: %s)\n', ...
    mag_filtered(2), mat2str(filter_info.fundamental_preserved));
fprintf('========================================\n\n');

% =======================================================================
% END OF FUNCTION
% =======================================================================

end
