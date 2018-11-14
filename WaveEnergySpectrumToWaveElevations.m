%% Initialize
clear ; clc ; close all ;

%% Given for ISSC wave spectrum
% Condition to represent the wave spectrum
waveFrequencyMax = 4 ;

% Variables to represent to the wave energy spectrum
periodAvegrage = 10 ; 
waveSignificantHeight = 6 ;
A = (173 * waveSignificantHeight^2) / (periodAvegrage^2) ; 
B = 691 / (waveSignificantHeight^4) ;
waveFrequency = 0:0.01:waveFrequencyMax ;

% Wave energy spectrum formula
waveEnergySpectrum = (A ./ waveFrequency.^5) .* exp(-B ./ waveFrequency.^5) ;

% Visualization of the wave energy spectrum
waveEnergySpectrumFig = figure ;
figure(waveEnergySpectrumFig) ;
plot(waveFrequency, waveEnergySpectrum) ;
title('ISSC wave spectrum') ;
xlabel('\omega (rad/sec)') ; ylabel('S_{\eta}(\omega) (m^{2}/(rad/sec))') ;
grid on ;
saveas(waveEnergySpectrumFig, 'WaveEnergySpectrumISSC.png') ;

%% Wave elevation time series according to discrete wave frequncy
% It Needs 
%  1. discrete wave frequency
%  2. discrete wave energy spectrum according to discrete wave frequecy
%  3. wave frequency interval
%  4. time step to represent the wave elevation in time series
%  5. random phase anlge according to discrete wave frequency

% the number of bins to devide the wave frequency domain for discrete anular frequnecy
numBins = 500 ;      

% 1. discrete wave frequency
waveFrequencyDiscrete = linspace(0, waveFrequencyMax, numBins) ;

% 2. discrete wave energy spectrum
waveEnergySpectrumDiscrete = interp1(waveFrequency, waveEnergySpectrum, waveFrequencyDiscrete) ;

% 3. Wave frequency interval
waveFreqeuncyInterval = waveFrequencyMax / numBins ;

% 4. Random phase angle according to discrete wave frequency
phaseAngleRandSeed = 2 ;
rng(phaseAngleRandSeed) ;
randomPhaseAngle = rand(1, numBins) * (2*pi) ;

% Wave elevation according to discrete wave frequncdy
% waveElevation:  across row: anglularFrequencyIndex, across column: timeStep
timeMax = 60 * 5 ;
timeStep = 0:0.01:timeMax ;
for waveFrequencyIndex = 1:length(waveFrequencyDiscrete)
    componentWaveElevation(:, waveFrequencyIndex) =...
        sqrt(2 * waveEnergySpectrumDiscrete(waveFrequencyIndex) * waveFreqeuncyInterval)...
        * cos(waveFrequencyDiscrete(waveFrequencyIndex) * timeStep...
        + randomPhaseAngle(waveFrequencyIndex)) ;
end


%% Visualization of the component wave(each regular wave) elevation according to discrete wave frequncy
showsComponentWaveElevation = false ;

if showsComponentWaveElevation
    numSubplotInOneFig = 24 ;
    numFig = ceil(numBins / numSubplotInOneFig) ;
    numSubplotColumn = 3 ;
    numSubplotRow = ceil(numSubplotInOneFig / numSubplotColumn) ;
    waveFrequencyIndex = 0 ;     

    for figIndex = 1:numFig
        componentWaveElevationFig(figIndex) = figure ;
        figure(componentWaveElevationFig(figIndex)) ;
        set(componentWaveElevationFig(figIndex), 'position', [0 0 (numSubplotColumn * 150 * 1.4) (numSubplotRow * 100 * 1.2)]) ;
        tempWaveFrequencyIndexHistory = [] ;     % reset the wave frequency index history

        for subplotCount = 1:numSubplotInOneFig
            waveFrequencyIndex = waveFrequencyIndex + 1 ;
            if waveFrequencyIndex > numBins
                break
            end
            tempWaveFrequencyIndexHistory(subplotCount) = waveFrequencyIndex ;    % Hold wave frequency index history for subgroup title

            subplot(numSubplotRow, numSubplotColumn, subplotCount) ;
            plot(timeStep, componentWaveElevation(:, waveFrequencyIndex)) ;
            axis([0 timeMax -max(max(componentWaveElevation)) max(max(componentWaveElevation))]) ;
            xlabel('t (sec)') ; ylabel('\eta(t) (m)') ;
            title(['\omega_{', num2str(waveFrequencyIndex), '}']) ;
            grid on ;
        end

        sgtitle(['Time singal of wave elevation for \omega_{',...
            num2str(min(tempWaveFrequencyIndexHistory)), '} _- _{',...
            num2str(max(tempWaveFrequencyIndexHistory)), '}']) ;

        componentWavefigSaveName = sprintf('WaveElevation%d_%d.png',...
            min(tempWaveFrequencyIndexHistory),...
            max(tempWaveFrequencyIndexHistory)) ;

        saveas(componentWaveElevationFig(figIndex), componentWavefigSaveName) ;
    end
end

%% Visualization of the irregular wave(sum of each regular wave) elevation
componentWaveElevation = fillmissing(componentWaveElevation, 'constant', 0) ;
waveElevation = sum(componentWaveElevation, 2) ;

waveElevationFig = figure ;
plot(timeStep, waveElevation) ;
axis([0 timeMax -15 15]) ;
xlabel('t (sec)') ; ylabel('\eta(t) (m)') ;
title('Wave elevation') ;
grid on ;
waveElevationSaveFigName = sprintf('waveElevation_Bins%d_RandSeed%d.png', numBins, phaseAngleRandSeed) ;
saveas(waveElevationFig, waveElevationSaveFigName) ;

waveElevationSaveMatName = sprintf('waveElevation_Bins%d_RandSeed%d.mat', numBins, phaseAngleRandSeed) ;
waveFreqeuncyDiscreteSaveMatName = sprintf('waveFrequency_Bins%d_RandSeed%d.mat', numBins, phaseAngleRandSeed) ;
timeStepSaveMatName = sprintf('timeStep_Bins%d_RandSeed%d.mat', numBins, phaseAngleRandSeed) ;
save(waveElevationSaveMatName, 'waveElevation') ;   
save(waveFreqeuncyDiscreteSaveMatName, 'waveFrequencyDiscrete') ;
save(timeStepSaveMatName, 'timeStep')