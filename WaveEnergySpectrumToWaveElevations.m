%% Initialize
clear ; clc ; close all ;

%% Given for ISSC wave spectrum
% Condition to represent the wave spectrum
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
waveEnergySpectrumFig = figure ;
figure(waveEnergySpectrumFig) ;
plot(angularFrequency, waveEnergySpectrum) ;
title('ISSC wave spectrum') ;
xlabel('\omega (rad/sec)') ; ylabel('S_{\eta}(\omega) (m^{2}/(rad/sec))') ;
grid on ;
saveas(waveEnergySpectrumFig, 'WaveEnergySpectrumISSC.png') ;

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
timeMax = 30 ;
timeStep = 0:0.01:timeMax ;
for angularFrequencyIndex = 1:length(angularFrequencyDiscrete)
    waveElevation(:, angularFrequencyIndex) =...
        sqrt(2 * waveEnergySpectrumDiscrete(angularFrequencyIndex) * angularFreqeuncyInterval)...
        * cos(angularFrequencyDiscrete(angularFrequencyIndex) * timeStep...
        + randomPhaseAngle(angularFrequencyIndex)) ;
end


%% Visualization of the wave elevation according to discrete angular wave frequncy
numSubplotInOneFig = 15 ;
numFig = ceil(numBins / numSubplotInOneFig) ;
numSubplotColumn = 3 ;
numSubplotRow = ceil(numSubplotInOneFig / numSubplotColumn) ;
angularFrequencyIndex = 0 ;     

for figIndex = 1:numFig
    waveElevationFig(figIndex) = figure ;
    figure(waveElevationFig(figIndex)) ;
    set(waveElevationFig(figIndex), 'position', [0 0 (numSubplotColumn * 300) (numSubplotRow * 150)]) ;
    tempAngularFrequencyIndexHistory = [] ;     % reset the angular frequency index history
    
    for subplotCount = 1:numSubplotInOneFig
        angularFrequencyIndex = angularFrequencyIndex + 1 ;
        if angularFrequencyIndex > numBins
            break
        end
        tempAngularFrequencyIndexHistory(subplotCount) = angularFrequencyIndex ;    % Save angular frequency index history for subgroup title
        
        subplot(numSubplotRow, numSubplotColumn, subplotCount) ;
        plot(timeStep, waveElevation(:, angularFrequencyIndex)) ;
        axis([0 timeMax -max(max(waveElevation)) max(max(waveElevation))]) ;
        xlabel('t (sec)') ; ylabel('\eta(t) (m)') ;
        title(['\omega_{', num2str(angularFrequencyIndex), '}']) ;
        grid on ;
    end

    sgtitle(['Time singal of wave elevation for \omega_{',...
        num2str(min(tempAngularFrequencyIndexHistory)), '} _- _{',...
        num2str(max(tempAngularFrequencyIndexHistory)), '}']) ;
    
    figSaveName = sprintf('waveElevation%d_%d.png',...
        min(tempAngularFrequencyIndexHistory),...
        max(tempAngularFrequencyIndexHistory)) ;
    
    saveas(waveElevationFig(figIndex), figSaveName) ;

end

