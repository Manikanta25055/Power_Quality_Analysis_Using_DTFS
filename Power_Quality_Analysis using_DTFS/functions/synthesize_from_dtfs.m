function [x_reconstructed, t] = synthesize_from_dtfs(X_k, fs, N, num_periods)
% SYNTHESIZE_FROM_DTFS - Reconstructs time-domain signal from DTFS coefficients
%
% DESCRIPTION:
%   Performs inverse DTFS to reconstruct the original periodic signal
%   from its DTFS coefficients. This is crucial for:
%   - Verifying DTFS analysis accuracy
%   - Reconstructing filtered signals after harmonic removal
%   - Understanding how individual harmonics contribute to waveform shape
%
% SYNTAX:
%   [x_reconstructed, t] = synthesize_from_dtfs(X_k, fs, N)
%   [x_reconstructed, t] = synthesize_from_dtfs(X_k, fs, N, num_periods)
%
% INPUTS:
%   X_k         - Complex DTFS coefficients [N x 1 vector]
%   fs          - Sampling frequency (Hz)
%   N           - Period length (samples per fundamental cycle)
%   num_periods - (Optional) Number of periods to generate (default: 1)
%
% OUTPUTS:
%   x_reconstructed - Reconstructed time-domain signal
%   t               - Time vector (seconds)
%
% MATHEMATICAL FOUNDATION:
%   DTFS Synthesis Equation:
%   x[n] = sum_{k=0}^{N-1} X[k] * exp(j*2*pi*k*n/N)
%
%   where:
%   - n = 0, 1, 2, ..., N-1 (time sample index)
%   - k = 0, 1, 2, ..., N-1 (harmonic index)
%   - X[k] = complex DTFS coefficient for kth harmonic
%
% POWER QUALITY APPLICATION:
%   After filtering out specific harmonics (e.g., setting X[3]=0 to remove
%   3rd harmonic), this function reconstructs the cleaned power signal.
%
% EXAMPLE USAGE:
%   % Perform DTFS analysis
%   [X_k, mag, phase, freq] = calculate_dtfs(x_distorted, fs, N);
%
%   % Remove 3rd harmonic by setting coefficient to zero
%   X_k_filtered = X_k;
%   X_k_filtered(4) = 0;  % Remove 3rd harmonic (k=3, index=4)
%
%   % Reconstruct filtered signal
%   [x_clean, t] = synthesize_from_dtfs(X_k_filtered, fs, N, 3);
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025
% COURSE: Digital Signal Processing (FISAC Assessment)
% PROJECT: Power Quality Analysis using DTFS

% =======================================================================
% INPUT VALIDATION
% =======================================================================

% Check if DTFS coefficients are provided
if nargin < 1
    error('synthesize_from_dtfs:InvalidInput', 'DTFS coefficients X_k are required');
end

% Set default sampling frequency if not provided
if nargin < 2 || isempty(fs)
    fs = 10000;  % Default 10 kHz sampling
    warning('synthesize_from_dtfs:DefaultFs', 'Using default fs = 10 kHz');
end

% Set default period length if not provided
if nargin < 3 || isempty(N)
    N = length(X_k);  % Use coefficient vector length as period
end

% Set default number of periods
if nargin < 4 || isempty(num_periods)
    num_periods = 1;  % Generate one period by default
end

% Validate DTFS coefficients
if ~isnumeric(X_k) || ~isvector(X_k)
    error('synthesize_from_dtfs:InvalidCoefficients', 'X_k must be a numeric vector');
end

% Ensure X_k is a column vector
X_k = X_k(:);

% Check coefficient length matches period
if length(X_k) ~= N
    warning('synthesize_from_dtfs:LengthMismatch', ...
        'Coefficient length (%d) does not match period (%d). Using length(X_k) as period.', ...
        length(X_k), N);
    N = length(X_k);
end

% Validate number of periods
if num_periods < 1 || num_periods > 100
    error('synthesize_from_dtfs:InvalidPeriods', ...
        'num_periods must be between 1 and 100');
end

% =======================================================================
% INVERSE DTFS CALCULATION
% =======================================================================

% Total number of samples to generate
total_samples = N * num_periods;

% Preallocate output signal
x_reconstructed = zeros(total_samples, 1);

% Generate time vector
t = (0:(total_samples-1))' / fs;

fprintf('\n========================================\n');
fprintf('INVERSE DTFS SYNTHESIS\n');
fprintf('========================================\n');
fprintf('Period Length: %d samples\n', N);
fprintf('Number of Periods: %d\n', num_periods);
fprintf('Total Samples: %d\n', total_samples);
fprintf('Sampling Frequency: %.0f Hz\n', fs);
fprintf('Synthesizing signal...\n');

% Inverse DTFS synthesis equation:
% x[n] = sum_{k=0}^{N-1} X[k] * exp(j*2*pi*k*n/N)

for n = 0:(total_samples-1)
    % Initialize accumulator for this time sample
    sum_value = 0;

    % Sum contributions from all harmonics
    for k = 0:(N-1)
        % Complex exponential: exp(j*2*pi*k*n/N)
        % Note: n is in range [0, total_samples-1], but we use modulo N
        % to handle multiple periods correctly
        n_mod = mod(n, N);
        complex_exponential = exp(1j * 2 * pi * k * n_mod / N);

        % Accumulate: X[k] * exp(j*2*pi*k*n/N)
        sum_value = sum_value + X_k(k+1) * complex_exponential;
    end

    % Store reconstructed sample
    x_reconstructed(n+1) = sum_value;
end

% The result should be real (imaginary part is numerical noise)
% Take real part to avoid any numerical issues
if max(abs(imag(x_reconstructed))) < 1e-10
    x_reconstructed = real(x_reconstructed);
else
    warning('synthesize_from_dtfs:ComplexResult', ...
        'Reconstructed signal has significant imaginary component (max = %.2e)', ...
        max(abs(imag(x_reconstructed))));
    % Still take real part, but warn user
    x_reconstructed = real(x_reconstructed);
end

% =======================================================================
% OUTPUT INFORMATION
% =======================================================================

% Calculate signal statistics
rms_value = sqrt(mean(x_reconstructed.^2));
peak_value = max(abs(x_reconstructed));
crest_factor = peak_value / rms_value;

fprintf('----------------------------------------\n');
fprintf('SYNTHESIS COMPLETE\n');
fprintf('----------------------------------------\n');
fprintf('Signal Statistics:\n');
fprintf('  RMS Value: %.4f V\n', rms_value);
fprintf('  Peak Value: %.4f V\n', peak_value);
fprintf('  Crest Factor: %.4f\n', crest_factor);
fprintf('  Duration: %.4f seconds\n', total_samples/fs);
fprintf('========================================\n\n');

% =======================================================================
% END OF FUNCTION
% =======================================================================

end
