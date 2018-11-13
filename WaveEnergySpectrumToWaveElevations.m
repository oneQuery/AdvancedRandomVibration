%% Initialize
clear ; clc ; close all ;

%% Given for ISSC wave spectrum
% Conditions to represent the wave spectrum
angularFrequencyMax = 4 ;

% Variables to represent to the wave energy spectrum
periodAvegrage = 10 ; 
waveSignificantHeight = 6 ;
A = (173 * waveSignificantHeight^2) / (periodAvegrage^2) ; 
B = 691 / (waveSignificantHeight^4) ;
angularFrequency = 0:0.01:angularFrequencyMax ;

% Wave energy spectrum formula
waveEnergySpectrum = (A ./ angularFrequency.^5) .* exp(-B ./ angularFrequency.^5) ;

% Visualization of the wave energy spectrum
figure(1) ;
plot(angularFrequency, waveEnergySpectrum) ;
title('ISSC wave spectrum') ;
xlabel('\omega (rad/sec)') ; ylabel('S_{\eta}(\omega) (m^{2}/(rad/sec))') ;
grid on ;

%% Wave elevation time series according to discrete angular frequncy
% It Needs 
%  1. discrete angular wave frequency
%  2. discrete wave energy spectrum according to discrete angular wave frequecy
%  3. angular frequency interval
%  4. time step to represent the wave elevation in time series
%  5. random phase anlge according to discrete angular frequency

% the number of bins to devide the angular frequency domain for discrete anular frequnecy
numBins = 50 ;      

% 1. discrete angular wave frequency
angularFrequencyDiscrete = linspace(0, angularFrequencyMax, numBins) ;

% 2. discrete wave energy spectrum
waveEnergySpectrumDiscrete = interp1(angularFrequency, waveEnergySpectrum, angularFrequencyDiscrete) ;

% 3. Angular frequency interval
angularFreqeuncyInterval = angularFrequencyMax / numBins ;

% 4. Random phase angle according to discrete angular frequency
phaseAngleRandSeed = 1 ;
rng(phaseAngleRandSeed) ;
randomPhaseAngle = rand(1, numBins) * (2*pi) ;

% Wave elevation according to discrete angular frequncdy
% waveElevation across row: anglularFrequencyIndex, across column: timeStep
timeStep = 0:0.01:30 ;
for angularFrequencyIndex = 1:length(angularFrequencyDiscrete)
    waveElevation(:, angularFrequencyIndex) =...
        sqrt(2 * waveEnergySpectrumDiscrete(angularFrequencyIndex) * angularFreqeuncyInterval)...
        * cos(angularFrequencyDiscrete(angularFrequencyIndex) * timeStep...
        + randomPhaseAngle(angularFrequencyIndex)) ;
end

% Visualization of the wave elevation accordingto discrete angular wave frequncy
figure(2) ;
suptitle('Time singal of wave elevation for all discrete angular frequency')
NumSubplotColumn = 5 ;
for angularFrequencyIndex = 1:numBins

    subplot(numBins / NumSubplotColumn, NumSubplotColumn, angularFrequencyIndex) ;
    plot(timeStep, waveElevation(:, angularFrequencyIndex)) ;
    angularFrequencyIndexStr = num2str(angularFrequencyIndex) ;
    title(['\omega_', angularFrequencyIndexStr]) ;
    grid on ;
    
end