%% Initialize
clear ; clc ; close all ;

%% Wave spectrum
% Wave scatter diagram
meanWavePeriodList = [10; 10; 10; 12; 12; 12; 14; 14; 14] ;
significantWaveHeightList = [4; 6; 8; 4; 6; 8; 4; 6; 8] ;
occurrenceOfSeaStateList = [500; 1000; 500; 250; 500; 250; 100; 200; 100] ;
T = table(meanWavePeriodList, significantWaveHeightList, occurrenceOfSeaStateList) ;

% Wave spectrum(row: spectrum index, column: frequency step)
A = 173 * (significantWaveHeightList.^2) ./ (meanWavePeriodList.^4) ;
B = 691 ./ (meanWavePeriodList.^4) ;
maxFrequency = 3 ;
frequencyInterval = 0.01 ;

frequencySteps = 0:frequencyInterval:maxFrequency ;
for waveSpectrumIndex = 1:height(T) 
    waveSpectrum(waveSpectrumIndex, :) =...
        A(waveSpectrumIndex) ./ (frequencySteps.^5)...
        .* exp(-B(waveSpectrumIndex) ./ (frequencySteps.^4)) ;
end
waveSpectrum(isnan(waveSpectrum)) = 0 ;

% Plot the wave spectrums
figure('Name', 'Wave spectrum') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(T)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(frequencySteps, waveSpectrum(waveSpectrumIndex, :)) ;
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([0 20]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\omega(rad/sec)', 'fontsize', 16) ;
    ylabel('S_\eta(\omega)(m^2/(rad/sec))', 'fontsize', 16) ;
end
saveas(gcf, 'waveSpectrum.png') ;
close ;

%% Wave amplitude
waveAmplitude = sqrt(2 .* waveSpectrum .* frequencyInterval) ;

% Plot the wave amplitudes
figure('Name', 'Wave amplitude') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(T)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(frequencySteps, waveAmplitude(waveSpectrumIndex, :)) ;
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([0 0.6]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\omega(rad/sec)', 'fontsize', 16) ;
    ylabel('A(\omega)(m)', 'fontsize', 16) ;
end
saveas(gcf, 'waveAmplitude.png') ;
close ;

%% Heave force amplitude
waterDensity = 1 ;
gravitationalAcceleration = 9.81 ;

% Cylinder characteristics
cylinderDiameter = 10 ;

% Heave force amplitude
heaveForceAmplitude =...
    waterDensity * gravitationalAcceleration * pi * (cylinderDiameter^2) ./ 4 .* waveAmplitude ;

% Plot the heave force amplitudes
figure('Name', 'Heave force amplitude') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(T)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(frequencySteps, heaveForceAmplitude(waveSpectrumIndex, :)) ;
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([0 480]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\omega(rad/sec)', 'fontsize', 16) ;
    ylabel('F_0(\omega)(N)', 'fontsize', 16) ;
end
saveas(gcf, 'heaveForceAmplitude.png') ;
close ;

%% Response amplitude
% Cylinder characteristics
cylinderMass = 5000 ;

% Pipe characteristics
pipeArea = 1 ;
pipeYoungsModulus = 20 ;        % Mpa
pipeLength = 100 ;
pipeDampingRatio = 0.1 ;
pipeElasticModulus = pipeArea * pipeYoungsModulus / pipeLength ;
pipeDampingCoefficient = 2 * pipeDampingRatio * sqrt(cylinderMass * pipeElasticModulus) ;

% Response amplitude
responseAmplitude = heaveForceAmplitude...
    ./ sqrt((pipeElasticModulus - cylinderMass .* (frequencySteps.^2)).^2 +...
    (pipeDampingCoefficient .* frequencySteps).^2) ;

% Plot the response amplitudes
figure('Name', 'Response amplitude') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(T)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(frequencySteps, responseAmplitude(waveSpectrumIndex, :)) ;
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([0 1]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\omega(rad/sec)', 'fontsize', 16) ;
    ylabel('X(\omega)(mm)', 'fontsize', 16) ;
end
saveas(gcf, 'responseAmplitude.png') ;
close ;

%% Stress amplitude
stressAmplitude = pipeYoungsModulus ./ pipeLength .* responseAmplitude ;

% Plot the stress amplitudes
figure('Name', 'Stress amplitude') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(T)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(frequencySteps, stressAmplitude(waveSpectrumIndex, :)) ;
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([0 0.2]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\omega(rad/sec)', 'fontsize', 16) ;
    ylabel('\Delta\sigma(\omega)(Mpa)', 'fontsize', 16) ;
end
saveas(gcf, 'stressAmplitude.png') ;
close ;

%% Probability density function
% Variance
variance = sum(stressAmplitude * frequencyInterval, 2) ;

% Probability density function
for waveSpectrumIndex = 1:height(T)
    PDF(waveSpectrumIndex, :) =...
        stressAmplitude(waveSpectrumIndex, :)...
        ./ (variance(waveSpectrumIndex).^2)...
        .* exp(-(stressAmplitude(waveSpectrumIndex, :).^2)...
        ./ (2 .* (variance(waveSpectrumIndex).^2))) ;
end
PDF = PDF ./ 100 ;

% Plot the probability density fucntions for stress amplitudes
figure('Name', 'Probability density function') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(T)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(stressAmplitude(waveSpectrumIndex, :), PDF(waveSpectrumIndex, :), '.');
    set(gca, 'FontSize', 16) ;
    grid on;
    xlim([0 0.2]) ;
    ylim([0 0.6]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\Delta\sigma(Mpa)', 'fontsize', 16) ;
    ylabel('f(/Mpa)', 'fontsize', 16) ;
end
saveas(gcf, 'ProbabilityDensityFunction.png') ;
close ;

%% Total damage
% predicted number of cycles to failure for stress range
loga = 12.164 ;
m = 3 ;
N = 10.^(loga - m .* log10(stressAmplitude)) ;

% Total damage
for waveSpectrumIndex = 1:height(T)
    tatalDamage(waveSpectrumIndex, :) =...
        occurrenceOfSeaStateList(waveSpectrumIndex)...
        * sum(PDF(waveSpectrumIndex, :) ./ N(waveSpectrumIndex, :) , 2) ;
end