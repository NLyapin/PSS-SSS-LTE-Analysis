clear all; close all; clc;

cellIdentity = 1;                
resourceBlocks = 6;             
samplingFreq = 2e6;            
prefixType = 'Normal';         
duplexType = 'FDD';            

baseStation.NDLRB = resourceBlocks;
baseStation.CyclicPrefix = prefixType;
baseStation.CellRefP = 1;       
baseStation.NCellID = cellIdentity;
baseStation.DuplexMode = duplexType;

fprintf('Generating PSS signal for Cell ID: %d\n', cellIdentity);

baseStation.NSubframe = 0;
pssData = ltePSS(baseStation);      
pssLocations = ltePSSIndices(baseStation, 0); 

resourceGrid = lteDLResourceGrid(baseStation);
resourceGrid(pssLocations) = pssData;

pssSignal = lteOFDMModulate(baseStation, resourceGrid);

fprintf('Generating SSS signal for Cell ID: %d\n', cellIdentity);

sssData = lteSSS(baseStation);      
sssLocations = lteSSSIndices(baseStation, 0); 

sssGrid = lteDLResourceGrid(baseStation);
sssGrid(sssLocations) = sssData;

sssSignal = lteOFDMModulate(baseStation, sssGrid);

fprintf('Creating complete OFDM signal with PSS and SSS\n');

combinedGrid = lteDLResourceGrid(baseStation);
combinedGrid(pssLocations) = pssData;
combinedGrid(sssLocations) = sssData;

dataLocations = setdiff(1:numel(combinedGrid), [pssLocations; sssLocations]);
modulationData = randi([0 1], length(dataLocations), 1) * 2 - 1; 
combinedGrid(dataLocations) = modulationData;

ofdmWaveform = lteOFDMModulate(baseStation, combinedGrid);

function [timeShift, corrPeak] = syncPSS(inputSignal, refPSS, baseStation)
    refGrid = zeros(size(lteDLResourceGrid(baseStation)));
    refLocations = ltePSSIndices(baseStation, 0);
    refGrid(refLocations) = refPSS;
    
    refSignal = lteOFDMModulate(baseStation, refGrid);
    
    corrResult = xcorr(inputSignal, refSignal);
    [corrPeak, maxIndex] = max(abs(corrResult));
    
    timeShift = (maxIndex - length(refSignal)) / 3;
end

testShift = 100;  
testWaveform = [zeros(testShift, 1); ofdmWaveform; zeros(100, 1)];

snrLevel = 1000;
noiseLevel = 10^(-snrLevel/10);
testWaveform = testWaveform + sqrt(noiseLevel/2) * (randn(size(testWaveform)) + 1i*randn(size(testWaveform)));

refPSS = ltePSS(baseStation);

[foundShift, peakAmplitude] = syncPSS(testWaveform, refPSS, baseStation);

fprintf('True offset: %d samples\n', testShift);
fprintf('Detected offset: %d samples\n', foundShift);
fprintf('Error: %d samples\n', foundShift - testShift);
fprintf('Correlation peak amplitude: %.2f\n', peakAmplitude);

fprintf('Generating visualization plots\n');

figure('Name', 'PSS Frequency Domain', 'Position', [100 100 600 400]);
plot(real(pssData), 'm.-', 'LineWidth', 1.5);
hold on; plot(imag(pssData), 'b.-', 'LineWidth', 1.5);
title('PSS Signal in Frequency Domain');
xlabel('Subcarrier Index'); ylabel('Amplitude');
legend('Real Part', 'Imag Part'); grid on;

figure('Name', 'SSS Frequency Domain', 'Position', [150 150 600 400]);
plot(real(sssData), 'm.-', 'LineWidth', 1.5);
hold on; plot(imag(sssData), 'b.-', 'LineWidth', 1.5);
title('SSS Signal in Frequency Domain');
xlabel('Subcarrier Index'); ylabel('Amplitude');
legend('Real Part', 'Imag Part'); grid on;

figure('Name', 'PSS Constellation', 'Position', [200 200 600 400]);
plot(real(pssData), imag(pssData), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
title('PSS Constellation Diagram');
xlabel('Real Component'); ylabel('Imag Component'); grid on; axis equal;

figure('Name', 'PSS Time Domain', 'Position', [250 250 600 400]);
tPSS = (0:length(pssSignal)-1) / samplingFreq * 1e6;
plot(tPSS, real(pssSignal), 'm-', 'LineWidth', 1);
hold on; plot(tPSS, imag(pssSignal), 'b-', 'LineWidth', 1);
title('PSS Signal in Time Domain');
xlabel('Time (µs)'); ylabel('Amplitude');
legend('Real Part', 'Imag Part'); grid on;

figure('Name', 'SSS Time Domain', 'Position', [300 300 600 400]);
tSSS = (0:length(sssSignal)-1) / samplingFreq * 1e6;
plot(tSSS, real(sssSignal), 'm-', 'LineWidth', 1);
hold on; plot(tSSS, imag(sssSignal), 'b-', 'LineWidth', 1);
title('SSS Signal in Time Domain');
xlabel('Time (µs)'); ylabel('Amplitude');
legend('Real Part', 'Imag Part'); grid on;

figure('Name', 'OFDM Signal', 'Position', [350 350 600 400]);
tOFDM = (0:length(ofdmWaveform)-1) / samplingFreq * 1e6;
plot(tOFDM, abs(ofdmWaveform), 'm-', 'LineWidth', 1);
title('OFDM Signal Amplitude');
xlabel('Time (µs)'); ylabel('Amplitude Magnitude');
grid on;

figure('Name', 'Correlation Function', 'Position', [400 400 600 400]);
refGridVis = zeros(size(lteDLResourceGrid(baseStation)));
refLocationsVis = ltePSSIndices(baseStation, 0);
refGridVis(refLocationsVis) = refPSS;
refSignalVis = lteOFDMModulate(baseStation, refGridVis);

corrVis = xcorr(testWaveform, refSignalVis);

Ncorr = length(corrVis);
Mref = length(refSignalVis);
Nsig = length(testWaveform);
lags = (-Mref+1:Nsig-1);
corrTime = lags / samplingFreq * 1e6;

if length(corrTime) ~= length(corrVis)
    fprintf('Dimension mismatch: corrTime=%d, corrVis=%d\n', ...
            length(corrTime), length(corrVis));
    corrTime = (0:Ncorr-1) / samplingFreq * 1e6;
end

plot(corrTime, abs(corrVis), 'm-', 'LineWidth', 1);
title('PSS Correlation Function');
xlabel('Offset (µs)'); ylabel('Correlation Magnitude');
grid on;

figure('Name', 'Correlation Peak Zoom', 'Position', [450 450 600 400]);
[~, peakIndex] = max(abs(corrVis));
peakRange = max(1, peakIndex-50):min(length(corrVis), peakIndex+50);
plot(corrTime(peakRange), abs(corrVis(peakRange)), 'm-', 'LineWidth', 2);
title('Zoomed Correlation Peak');
xlabel('Offset (µs)'); ylabel('Correlation Magnitude');
grid on;

figure('Name', 'Power Spectral Density', 'Position', [500 500 600 400]);
[psdData, freq] = pwelch(ofdmWaveform, [], [], [], samplingFreq);
plot(freq/1e6, 10*log10(psdData), 'm-', 'LineWidth', 1);
title('Signal Power Spectral Density');
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
grid on;

fprintf('Converting to signed char format for HackRF\n');

maxAmplitude = max([max(abs(real(ofdmWaveform))), max(abs(imag(ofdmWaveform)))]);
normalizedOFDM = ofdmWaveform / maxAmplitude;

scaleValue = 127;
scaledOFDM = normalizedOFDM * scaleValue;

iqData = zeros(2*length(scaledOFDM), 1);
iqData(1:2:end) = real(scaledOFDM);   
iqData(2:2:end) = imag(scaledOFDM);   

iqInt8 = int8(round(iqData));

fprintf('Saving output files\n');

pssFile = sprintf('pss_id_%d.mat', cellIdentity);
save(pssFile, 'pssData', 'pssSignal', 'cellIdentity');
fprintf('PSS saved to: %s\n', pssFile);

sssFile = sprintf('sss_id_%d.mat', cellIdentity);
save(sssFile, 'sssData', 'sssSignal', 'cellIdentity');
fprintf('SSS saved to: %s\n', sssFile);

ofdmFile = sprintf('ofdm_wave_id_%d.mat', cellIdentity);
save(ofdmFile, 'ofdmWaveform', 'baseStation', 'samplingFreq');
fprintf('OFDM signal saved to: %s\n', ofdmFile);

hackrfFile = sprintf('ofdm_hackrf_id_%d.iq', cellIdentity);
fileHandle = fopen(hackrfFile, 'wb');
fwrite(fileHandle, iqInt8, 'int8');
fclose(fileHandle);
fprintf('HackRF data saved to: %s\n', hackrfFile);
fprintf('Format: signed 8-bit interleaved I/Q, sampling rate: %.0f Hz\n', samplingFreq);

function testRealSignal(fileName, samplingFreq, baseStation, refPSS)
    fprintf('\n=== REAL SIGNAL TESTING ===\n');
    fprintf('Loading file: %s\n', fileName);
    
    fileHandle = fopen(fileName, 'rb');
    if fileHandle == -1
        fprintf('Error: cannot open file %s\n', fileName);
        return;
    end
    
    iqRaw = fread(fileHandle, Inf, 'int8');
    fclose(fileHandle);
    
    if mod(length(iqRaw), 2) ~= 0
        iqRaw = iqRaw(1:end-1);
    end
    
    Icomp = double(iqRaw(1:2:end));
    Qcomp = double(iqRaw(2:2:end));
    inputSignal = complex(Icomp, Qcomp);
    
    inputSignal = inputSignal / 127;
    
    fprintf('Loaded %d complex samples\n', length(inputSignal));
    
    [timeShift, corrPeak] = syncPSS(inputSignal, refPSS, baseStation);
    
    fprintf('Synchronization results:\n');
    fprintf('Detected time shift: %d samples\n', timeShift);
    fprintf('Correlation peak value: %.4f\n', corrPeak);
    
    figure('Name', 'Real Signal Synchronization', 'Position', [200 200 1200 600]);
    
    subplot(2,3,1);
    plot(real(inputSignal(1:min(1000, end))), 'm-');
    hold on; plot(imag(inputSignal(1:min(1000, end))), 'b-');
    title('Received Signal (First 1000 Samples)');
    xlabel('Sample'); ylabel('Amplitude');
    legend('I Comp', 'Q Comp'); grid on;
    
    subplot(2,3,2);
    plot(abs(inputSignal(1:min(1000, end))), 'm-');
    title('Received Signal Amplitude');
    xlabel('Sample'); ylabel('Amplitude Magnitude');
    grid on;
    
    subplot(2,3,3);
    scatter(real(inputSignal(1:10:min(10000, end))), ...
            imag(inputSignal(1:10:min(10000, end))), 1, 'm.');
    title('Signal Constellation');
    xlabel('I Comp'); ylabel('Q Comp'); grid on; axis equal;
    
    refGrid = zeros(size(lteDLResourceGrid(baseStation)));
    refLocations = ltePSSIndices(baseStation, 0);
    refGrid(refLocations) = refPSS;
    refSignal = lteOFDMModulate(baseStation, refGrid);
    
    if length(inputSignal) > 50000
        corrResult = [];
        stepSize = 10000;
        for i = 1:stepSize:length(inputSignal)-length(refSignal)
            segment = inputSignal(i:min(i+stepSize+length(refSignal)-1, end));
            corrSegment = xcorr(segment, refSignal);
            corrResult = [corrResult; abs(corrSegment)];
        end
    else
        corrFull = xcorr(inputSignal, refSignal);
        corrResult = abs(corrFull);
    end
    
    subplot(2,3,[4,5,6]);
    plot(corrResult, 'm-');
    title('Synchronization Correlation Function');
    xlabel('Position'); ylabel('Correlation Magnitude');
    grid on;
    
    [peaks, locations] = findpeaks(corrResult, 'MinPeakHeight', max(corrResult)*0.3);
    hold on;
    plot(locations, peaks, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    legend('Correlation', 'Detected Peaks');
    
    fprintf('Found %d possible synchronization positions\n', length(peaks));
end