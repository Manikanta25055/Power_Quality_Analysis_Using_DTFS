function [X_k, magnitude, phase, frequencies] = calculate_dtfs(x, fs, N)
% CALCULATE_DTFS - Calculates Discrete Time Fourier Series coefficients
%
% DESCRIPTION:
%   Performs DTFS analysis on a periodic signal to extract harmonic content.
%   This is the core function for power quality analysis, decomposing
%   distorted power waveforms into their constituent harmonics.
%
% SYNTAX:
%   [X_k, magnitude, phase, frequencies] = calculate_dtfs(x, fs, N)
%
% INPUTS:
%   x  - Input periodic signal (voltage or current waveform)
%   fs - Sampling frequency (Hz) - typically 10 kHz for power signals
%   N  - Period length (number of samples per fundamental cycle)
%
% OUTPUTS:
%   X_k        - Complex DTFS coefficients [N x 1 vector]
%   magnitude  - Magnitude spectrum |X[k]| in Volts [N x 1 vector]
%   phase      - Phase spectrum in radians [N x 1 vector]
%   frequencies- Corresponding frequency values in Hz [N x 1 vector]
%
% MATHEMATICAL FOUNDATION:
%   DTFS Analysis Equation:
%   X[k] = (1/N) * sum_{n=0}^{N-1} x[n] * exp(-j*2*pi*k*n/N)
%
%   where:
%   - k = 0, 1, 2, ..., N-1 (harmonic index)
%   - n = 0, 1, 2, ..., N-1 (time sample index)
%   - X[k] = complex DTFS coefficient for kth harmonic
%
% POWER QUALITY INTERPRETATION:
%   - X[0] = DC component (should be ~0 for AC power)
%   - X[1] = Fundamental component (50 Hz)
%   - X[k] for k>1 = Harmonic components (k*50 Hz)
%
% EXAMPLE USAGE:
%   fs = 10000;              % 10 kHz sampling
%   f0 = 50;                 % 50 Hz fundamental
%   N = fs/f0;               % 200 samples per cycle
%   t = (0:N-1)/fs;
%   x = 230*sqrt(2)*sin(2*pi*50*t);  % Pure 50Hz signal
%   [X_k, mag, phase, freq] = calculate_dtfs(x, fs, N);
%
% AUTHOR: Manikanta Gonugondla
% DATE: October 2025
% COURSE: Digital Signal Processing (FISAC Assessment)
% PROJECT: Power Quality Analysis using DTFS

% =======================================================================
% INPUT VALIDATION
% =======================================================================

% Check if input signal is provided
if nargin < 1
    error('calculate_dtfs:InvalidInput', 'Input signal x is required');
end

% Set default sampling frequency if not provided
if nargin < 2 || isempty(fs)
    fs = 10000;  % Default 10 kHz sampling for power signals
    warning('calculate_dtfs:DefaultFs', 'Using default fs = 10 kHz');
end

% Set default period length if not provided
if nargin < 3 || isempty(N)
    N = length(x);  % Use full signal length as one period
    warning('calculate_dtfs:DefaultN', 'Using signal length as period: N = %d', N);
end

% Validate input signal
if ~isnumeric(x) || ~isvector(x)
    error('calculate_dtfs:InvalidSignal', 'Input x must be a numeric vector');
end

% Ensure x is a column vector
x = x(:);

% Check if signal length matches period
if length(x) < N
    error('calculate_dtfs:SignalTooShort', ...
        'Signal length (%d) must be >= period length (%d)', length(x), N);
end

% Extract one period if signal is longer
if length(x) > N
    x = x(1:N);
    warning('calculate_dtfs:SignalTruncated', ...
        'Using only first %d samples (one period)', N);
end

% =======================================================================
% DTFS COEFFICIENT CALCULATION
% =======================================================================

% Preallocate output arrays for efficiency
X_k = zeros(N, 1);

% Calculate DTFS coefficients using the analysis equation
% X[k] = (1/N) * sum_{n=0}^{N-1} x[n] * exp(-j*2*pi*k*n/N)

for k = 0:(N-1)
    % Initialize accumulator for this harmonic
    sum_value = 0;

    % Sum over all time samples
    for n = 0:(N-1)
        % Complex exponential: exp(-j*2*pi*k*n/N)
        complex_exponential = exp(-1j * 2 * pi * k * n / N);

        % Accumulate: x[n] * exp(-j*2*pi*k*n/N)
        sum_value = sum_value + x(n+1) * complex_exponential;
    end

    % Normalize by period length
    X_k(k+1) = sum_value / N;
end

% =======================================================================
% MAGNITUDE AND PHASE EXTRACTION
% =======================================================================

% Calculate magnitude spectrum (in Volts or Amps)
magnitude = abs(X_k);

% For power analysis, scale DC and fundamental appropriately
% DC component (k=0) stays as is
% Fundamental and harmonics: multiply by 2 to get peak amplitude
% (This accounts for positive and negative frequency components)
magnitude(2:end) = 2 * magnitude(2:end);

% Calculate phase spectrum (in radians)
phase = angle(X_k);

% =======================================================================
% FREQUENCY AXIS GENERATION
% =======================================================================

% Generate frequency vector
% Frequency resolution: df = fs/N Hz
df = fs / N;

% Frequency values corresponding to each DTFS coefficient
% freq[k] = k * df = k * (fs/N) Hz
frequencies = (0:(N-1))' * df;

% =======================================================================
% OUTPUT INFORMATION (for debugging/verification)
% =======================================================================

% Display key information
fprintf('\n========================================\n');
fprintf('DTFS ANALYSIS COMPLETE\n');
fprintf('========================================\n');
fprintf('Sampling Frequency: %.0f Hz\n', fs);
fprintf('Period Length: %d samples\n', N);
fprintf('Frequency Resolution: %.2f Hz\n', df);
fprintf('DC Component: %.4f V\n', magnitude(1));

% Find fundamental frequency (assuming k=1 is fundamental)
if N > 1
    fprintf('Fundamental Frequency: %.2f Hz\n', frequencies(2));
    fprintf('Fundamental Magnitude: %.4f V (peak)\n', magnitude(2));
    fprintf('Fundamental RMS: %.4f V\n', magnitude(2)/sqrt(2));
end

fprintf('========================================\n\n');

% =======================================================================
% END OF FUNCTION
% =======================================================================

end
