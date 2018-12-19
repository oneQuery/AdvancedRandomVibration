%% Initialize
clc ; clear ; close all ; 

%% In frequency domain
%% Wave spectrum
% Wave scatter diagram
meanWavePeriodList = [10; 10; 10; 12; 12; 12; 14; 14; 14] ;
significantWaveHeightList = [4; 6; 8; 4; 6; 8; 4; 6; 8] ;
occurrenceOfSeaStateList = [500; 1000; 500; 250; 500; 250; 100; 200; 100] ;
waveScatterDiagram = table(meanWavePeriodList, significantWaveHeightList, occurrenceOfSeaStateList) ;

% Wave spectrum(row: spectrum index, column: frequency step)
A = 173 * (significantWaveHeightList.^2) ./ (meanWavePeriodList.^4) ;
B = 691 ./ (meanWavePeriodList.^4) ;
maxFrequency = 3 ;
frequencyInterval = 0.01 ;

frequencySteps = 0:frequencyInterval:maxFrequency ;
for waveSpectrumIndex = 1:height(waveScatterDiagram) 
    waveSpectrum(waveSpectrumIndex, :) =...
        A(waveSpectrumIndex) ./ (frequencySteps.^5)...
        .* exp(-B(waveSpectrumIndex) ./ (frequencySteps.^4)) ;
end
waveSpectrum(isnan(waveSpectrum)) = 0 ;

% Plot the wave spectrums
figure('Name', 'Wave spectrum') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(waveScatterDiagram)
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
for waveSpectrumIndex = 1:height(waveScatterDiagram)
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
waterDensity = 1025 ;
gravitationalAcceleration = 9.81 ;

% Cylinder characteristics
cylinderDiameter = 10 ;

% Heave force amplitude
heaveForceAmplitude =...
    waterDensity * gravitationalAcceleration * pi * (cylinderDiameter^2) ./ 4 .* waveAmplitude ;

% Plot the heave force amplitudes
figure('Name', 'Heave force amplitude') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(frequencySteps, heaveForceAmplitude(waveSpectrumIndex, :)) ;
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([0 5E+5]) ;
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
pipeYoungsModulus = 20E+6 ;        % pa
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
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(frequencySteps, responseAmplitude(waveSpectrumIndex, :)) ;
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([0 3]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\omega(rad/sec)', 'fontsize', 16) ;
    ylabel('X(\omega)(m)', 'fontsize', 16) ;
end
saveas(gcf, 'responseAmplitude.png') ;
close ;

%% Stress amplitude
stressAmplitude = pipeYoungsModulus ./ pipeLength .* responseAmplitude ;

% Plot the stress amplitudes
figure('Name', 'Stress amplitude') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(frequencySteps, stressAmplitude(waveSpectrumIndex, :)) ;
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([0 5E+5]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\omega(rad/sec)', 'fontsize', 16) ;
    ylabel('\Delta\sigma(\omega)(pa)', 'fontsize', 16) ;
end
saveas(gcf, 'stressAmplitude.png') ;
close ;

%% Probability density function
% Variance
variance = sum(stressAmplitude * frequencyInterval, 2) ;

% Probability density function
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    PDF(waveSpectrumIndex, :) =...
        stressAmplitude(waveSpectrumIndex, :)...
        ./ (variance(waveSpectrumIndex).^2)...
        .* exp(-(stressAmplitude(waveSpectrumIndex, :).^2)...
        ./ (2 .* (variance(waveSpectrumIndex).^2))) ;
end
% PDF = PDF ./ 100 ;

% Plot the probability density fucntions for stress amplitudes
figure('Name', 'Probability density function') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(stressAmplitude(waveSpectrumIndex, :), PDF(waveSpectrumIndex, :), '.');
    set(gca, 'FontSize', 16) ;
    grid on;
    xlim([0 1E+6]) ;
    ylim([0 0.7E-5]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\Delta\sigma(pa)', 'fontsize', 16) ;
    ylabel('f(/pa)', 'fontsize', 16) ;
end
saveas(gcf, 'probabilityDensityFunction.png') ;
close ;

%% Total damage
% predicted number of cycles to failure for stress range
loga = 12.164 ;
m = 3 ;
N = 10.^(loga - m .* log10(stressAmplitude)) ;

% Total damage
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    totalDamage(waveSpectrumIndex, :) =...
        occurrenceOfSeaStateList(waveSpectrumIndex)...
        * sum(PDF(waveSpectrumIndex, :) ./ N(waveSpectrumIndex, :) , 2) ;
end

%% In time domain
timeStepInterval = 0.1 ;
maxTime = 1000 ;
timeSteps = 0:timeStepInterval:maxTime ;
%% Inpulse response
cylinderNaturalFrequency = sqrt(pipeElasticModulus / cylinderMass) ;
dampedNaturalFrequency = cylinderNaturalFrequency * sqrt(1 - pipeDampingRatio^2) ;

impulseResponseFunction =...
    (exp(-pipeDampingRatio .* cylinderNaturalFrequency .* timeSteps)...
    ./ (cylinderMass .* dampedNaturalFrequency))...
    .* sin(dampedNaturalFrequency .* timeSteps) ;

%% Heave force(row: spectrum index, column: time step)
heaveForceForEachFrequency = cell(height(waveScatterDiagram), length(frequencySteps)) ;
rng('default') ; rng(1) ;
randomPhaseAngleForEachFrequency = 2 * pi * rand(height(waveScatterDiagram), length(frequencySteps)) ;

for waveSpectrumIndex = 1:height(waveScatterDiagram)
    for frequencyStepIndex = 1:length(frequencySteps)
        heaveForceForEachFrequency{waveSpectrumIndex, frequencyStepIndex} =...
            heaveForceAmplitude(waveSpectrumIndex, frequencyStepIndex)...
            .*cos(frequencySteps(frequencyStepIndex) .* timeSteps...
            + randomPhaseAngleForEachFrequency(waveSpectrumIndex, frequencyStepIndex)) ;
    end
end

for waveSpectrumIndex = 1:height(waveScatterDiagram)
    heaveForceInTimeDomain(waveSpectrumIndex, :) = sum(cell2mat(heaveForceForEachFrequency(waveSpectrumIndex, :)')) ;
end

% Plot the heave force in time domain
figure('Name', 'Heave force in time domain') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(timeSteps, heaveForceInTimeDomain(waveSpectrumIndex, :));
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([-5E+6 5E+6]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('t(sec)', 'fontsize', 16) ;
    ylabel('F_z(N)', 'fontsize', 16) ;
end
saveas(gcf, 'heaveForceInTimeDomain.png') ;
close ;

%% Stress in time domain
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    stressInTimeDomain(waveSpectrumIndex, :) = (pipeYoungsModulus / pipeLength)...
        * conv(heaveForceInTimeDomain(waveSpectrumIndex, :), impulseResponseFunction)...
        * timeStepInterval ;
end

% Plot the heave force in time domain
figure('Name', 'Stress in time domain') ;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    subplot(3, 3, waveSpectrumIndex) ;
    plot(timeSteps, stressInTimeDomain(waveSpectrumIndex, 1:length(timeSteps)));
    set(gca, 'FontSize', 16) ;
    grid on;
    ylim([-5E+6 5E+6]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('t(sec)', 'fontsize', 16) ;
    ylabel('\sigma(pa)', 'fontsize', 16) ;
end
saveas(gcf, 'stressInTimeDomain.png') ;
close ;

%% Occurrence
rf = cell(height(waveScatterDiagram), 1) ;
for waveSpectrumIndex = 1:height(waveScatterDiagram)
    rf{waveSpectrumIndex} = rainflow(stressInTimeDomain(waveSpectrumIndex, :), timeStepInterval) ;
end

% Plot the stress amplitude vs number of cycles
% figure('Name', 'Histogram') ;
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]) ;
% for waveSpectrumIndex = 1:height(waveScatterDiagram)
%     subplot(3, 3, waveSpectrumIndex) ;
%     plot(timeSteps, rf{waveSpectrumIndex});
%     set(gca, 'FontSize', 16) ;
%     grid on;
%     ylim([-5E+6 5E+6]) ;
%     title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
%         'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
%     xlabel('t(sec)', 'fontsize', 16) ;
%     ylabel('\sigma(pa)', 'fontsize', 16) ;
% end
% saveas(gcf, 'stressInTimeDomain.png') ;
% close ;