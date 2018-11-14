%% Initialize
clear ; clc ; close all ;

%% Load time signal of wave elevation
Bins = 500 ;
RandSeed = 2 ;

waveElevationMatName = sprintf('waveElevation_Bins%d_RandSeed%d.mat', Bins, RandSeed) ;
waveFreqeuncyDiscreteMatName = sprintf('waveFrequency_Bins%d_RandSeed%d.mat', Bins, RandSeed) ;
timeStepMatName = sprintf('timeStep_Bins%d_RandSeed%d.mat', Bins, RandSeed) ;

load(waveElevationMatName) ;
load(waveFreqeuncyDiscreteMatName) ;
load(timeStepMatName) ;

% Syncronize axis. across row: time, across column: wave frequency
waveElevation = waveElevation' ;
waveFrequencyDiscrete = waveFrequencyDiscrete' ;

periodWhole = max(timeStep) ;

%% Calculate coefficient A
for waveFrequencyDiscreteIndex = 1:length(waveFrequencyDiscrete)
    ADiscrete(waveFrequencyDiscreteIndex) =...
        (2 / periodWhole) *...
        trapz(timeStep, waveElevation .* cos(waveFrequencyDiscrete(waveFrequencyDiscreteIndex)*timeStep)) ;
end

%% Calculate coefficient B
for waveFrequencyDiscreteIndex = 1:length(waveFrequencyDiscrete)
    BDiscrete(waveFrequencyDiscreteIndex) =...
        (2 / periodWhole) *...
        trapz(timeStep, waveElevation .* sin(waveFrequencyDiscrete(waveFrequencyDiscreteIndex)*timeStep)) ;
end

%% Calculate wave height and phase angle
waveHeightDiscreteSquare = ADiscrete.^2 + BDiscrete.^2 ;
% phaseAngle = - BDiscrete ./ ADiscrete ;

%% Derive wave spectrum
waveFrequencyInterval = waveFrequencyDiscrete(2) - waveFrequencyDiscrete(1) ;
waveSpectrumDiscrete = waveHeightDiscreteSquare / (2*waveFrequencyInterval) ;

%% Visulization of the wave spectrum
waveSpectrumFig = figure ;
plot(waveFrequencyDiscrete, waveSpectrumDiscrete) ;
xlabel('\omega (rad/sec)') ; ylabel('S_{\eta}(\omega) (m^{2}/(rad/sec))') ;
title('Wave spectrum') ;
grid on ;
waveSpectrumFigName = sprintf('waveSpectrum_Bins%d_RandSeed%d.png', Bins, RandSeed) ;
saveas(waveSpectrumFig, waveSpectrumFigName) ;

%% Significant wave height
variance = trapz(waveFrequencyDiscrete, waveSpectrumDiscrete) ;
waveSignificantHeight = 4 * sqrt(variance)