%% functions

function plv = calculatePLV(signal1, signal2, fs, frequencyRange)
    % signal1: Time series data for channel 1
    % signal2: Time series data for channel 2
    % fs: Sampling frequency (in Hz)
    % frequencyRange: Frequency range of interest [lowerBound, upperBound] (in Hz)

    % Extract frequency range of interest using a bandpass filter
    filteredSignal1 = bandpass(signal1, frequencyRange, fs);
    filteredSignal2 = bandpass(signal2, frequencyRange, fs);

    % Compute the Hilbert transform to obtain instantaneous phases
    phase1 = angle(hilbert(filteredSignal1));
    phase2 = angle(hilbert(filteredSignal2));

    % Compute phase differences
    phaseDiff = phase2 - phase1;

    % Compute PLV
    plv = abs(mean(exp(1i * phaseDiff)));
end