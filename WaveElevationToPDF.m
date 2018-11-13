%% Initialize
clear ; clc ; close all ;

%% Add path of function
addpath('rainflow') ;

%% Load the wave elevation
load('waveElevation_Bins500_RandSeed2.mat') ;
load('timeStep.mat') ;

%% Analyze the wave elevation
dataWaveElevation = rainflow(waveElevation, timeStep) ;

discreteApmlitude = dataWaveElevation(1, :) ;
discreteWaveHeightMean = dataWaveElevation(2, :) ;
numWaveHeightOccurance = dataWaveElevation(3, :) ;

%% Probability densify function
% Wave height average
waveHeightAverageSquare = sum(discreteWaveHeightMean.^2 .* numWaveHeightOccurance) / sum(numWaveHeightOccurance) ;

% Probability density function
probabilityDensityFunction = (2 * discreteWaveHeightMean) / (waveHeightAverageSquare) .* exp((-(discreteWaveHeightMean.^2) / waveHeightAverageSquare)) ;

%% Visualization of the probability density function
PDFFig = figure ;
plot(discreteWaveHeightMean, probabilityDensityFunction, '.') ;
axis([0 max(discreteWaveHeightMean)*1.2 0 0.25]) ;
title('Probability density function') ;
xlabel('Wave height (m)') ; ylabel('Probability density function (%/m)') ;
grid on ;
saveas(PDFFig, 'PDFCase6.png') ;
