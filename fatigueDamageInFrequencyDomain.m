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
    ylim([0 1]) ;
    title(['H_1_/_3=' num2str(significantWaveHeightList(waveSpectrumIndex))...
        'm, T_m_e_a_n=' num2str(meanWavePeriodList(waveSpectrumIndex)) 'sec']) ;
    xlabel('\omega(rad/sec)', 'fontsize', 16) ;
    ylabel('A(\omega)(m)', 'fontsize', 16) ;
end
saveas(gcf, 'waveAmplitude.png') ;
close ;

%% Heave force
