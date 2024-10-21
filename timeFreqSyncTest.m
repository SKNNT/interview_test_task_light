% close all;
clear all;
clc;

%% Experiment's parameters
blockLength = 500;
dataBlockCount = 10;
sampleRate = 200e3;
toneFreq = 50e3;
freqError = 10e3;
snrDbRange = -20:2:20;
seedCount = 1000;

%% Generate transmitted signal
% Frequency correction sequence: single-tone with frequency toneFreq
freqSeq = exp(1j*2*pi*toneFreq/sampleRate*(1:blockLength));

% Synchronization sequence - some random complex signal normalized to 1
% Used for correlation
rng(1);
synchSeq = randn(1, blockLength) + 1j*randn(1, blockLength);
synchSeq = synchSeq / sqrt(mean(abs(synchSeq).^2));

% Data sequence - BPSK modulated signal: 1 = bit 1, -1 = bit 0
dataSeq = randi([0 1], 1, blockLength*dataBlockCount) * 2 - 1;

% Compose frame
txSignal = [freqSeq synchSeq dataSeq];

%% Run experiment
syncFailed = zeros(seedCount, length(snrDbRange));
% Iterate over all SNRs
for snrDbIdx = 1:length(snrDbRange)
    snrDb = snrDbRange(snrDbIdx);
    
    % Iterate over all seeds
    for seedIdx = 1:seedCount
        rng(seedIdx);
        
        %% Propagate through the channel
        % Apply frequency shift due to carrier mismatch
        rxSignal = txSignal .* exp(1j*2*pi*freqError/sampleRate*(1:length(txSignal)));
        
        % Sample AWGN noise from CN(0,1) distribution
        noise = (randn(size(rxSignal)) + 1j * randn(size(rxSignal))) * sqrt(1/2);
        
        % Add noise with given snr
        rxSignal = rxSignal + noise * db2mag(-snrDb);
        
        %% Synchronization
        % Phase Difference Calculation
        phaseDifferences = angle(rxSignal(2:length(rxSignal)) .* conj(rxSignal(1:length(rxSignal)-1)));
        
        % Calculate moving average of phaseDifferences using a window size equal to length of freqSeq
        windowSize = length(freqSeq);
        movingAvg = movmean([zeros(1, windowSize), phaseDifferences, zeros(1, windowSize)], windowSize);
        
        % Peak movingAvg finding
        [~,movingAvgMaxIdx] = max(abs(movingAvg));
        
        % Skip it if we went out of the array rxSignal
        if ( movingAvgMaxIdx + windowSize/2 > length(rxSignal))
            syncFailed(seedIdx, snrDbIdx) = 1;
            continue;
        end
        
        % FFT Analysis
        fftLength = movingAvgMaxIdx - windowSize/2;
        fftAmplitudeSignal = abs(fft(rxSignal(1:fftLength)));
        
        % Peak fft signal finding
        [~, fftAmplitudeSignalMaxIdx] = max(fftAmplitudeSignal);
        
        % Finding the shift
        fftFrequencies = (0:fftLength-1) * (sampleRate / fftLength);
        frequencyShift = fftFrequencies(fftAmplitudeSignalMaxIdx) - toneFreq;
        
        % Frequency compensation
        rxSignal = rxSignal .* exp(1j*2*pi*(-frequencyShift)/sampleRate*(1:length(rxSignal)));
        
        % Correlation
        correlation = conv(rxSignal, conj(synchSeq(end:-1:1)));
        
        % Peak correlation finding
        [~,correlationMaxIdx] = max(abs(correlation));
        
        % Finding sync start
        syncStart = correlationMaxIdx - length(synchSeq) + 1;
        
        % Calculate squared timing error compared to reference sync start
        refSyncStart = length(freqSeq) + 1;
        syncFailed(seedIdx, snrDbIdx) = abs(syncStart - refSyncStart) > 3;
    end
end

%% Plot results
fig = figure;semilogy(snrDbRange, mean(syncFailed,1),'Linewidth',2);
xlabel("SNR, dB"); ylabel("Probability of failed sync");
ylim([1e-3 1])
grid on; title("Initial frequency mismatch is "+string(round(freqError/1e3))+"KHz");
saveas(fig, "res_"+string(round(freqError/1e3))+".png");

