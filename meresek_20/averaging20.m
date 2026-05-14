clear
clc
close all

% alapadatok, eleresi utak beallitasa
parentDir = 'D:\uni\s8\optomechatronika projekt\meresek_20';
filePattern = fullfile(parentDir, '**', '*.csv');
allFiles = dir(filePattern);

% kimeneti mappa az abraknak
outputDir = 'D:\uni\s8\optomechatronika projekt\kepek\v0.3\ido_frekvenciat_abrak';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% ures tablazat az eredmenyeknek
results = table();

fprintf('%d db fajl talalhato osszesen.\n', length(allFiles));

currentMouse = '';
subPlotIdx = 0;

% mintaveteli ido es frekvencia
dt = 0.4883 / 1000;   % s
Fs = 1 / dt;          % Hz

% frekvenciatartomanyos abra beallitas
maxFreqToPlot = 80;   % Hz
targetFreq = 4;       % Hz

% fix egerenkenti abra-elrendezes
numRows = 2;
numCols = 2;

for k = 1:length(allFiles)

    folderName = allFiles(k).folder;
    fileName = allFiles(k).name;
    fullFilePath = fullfile(folderName, fileName);

    % egerek azonositoja a mappa nevebol
    [~, mouseID] = fileparts(folderName);

    % uj eger eseten uj abrapar
    if ~strcmp(mouseID, currentMouse)

        % elozo eger kesz abrainak mentese
        if ~isempty(currentMouse)
            exportgraphics(timeFig, fullfile(outputDir, [currentMouse '_time.png']), 'Resolution', 300);
            exportgraphics(freqFig, fullfile(outputDir, [currentMouse '_frequency.png']), 'Resolution', 300);

            fprintf('Abrak mentve: %s\n', currentMouse);
        end

        currentMouse = mouseID;

        % idotartomanyos figure
        figure('Name', ['Idotartomany - Eger: ' mouseID], ...
               'NumberTitle', 'off', ...
               'Units', 'normalized', ...
               'Position', [0.05 0.05 0.85 0.8]);

        timeFig = gcf;

        % frekvenciatartomanyos figure
        figure('Name', ['Frekvenciatartomany - Eger: ' mouseID], ...
               'NumberTitle', 'off', ...
               'Units', 'normalized', ...
               'Position', [0.08 0.08 0.85 0.8]);

        freqFig = gcf;

        subPlotIdx = 1;
    end

    fprintf('Feldolgozas: [%s] -> %s\n', mouseID, fileName);

    try
        % beolvasas
        rawData = readmatrix(fullFilePath, 'Delimiter', ';', 'NumHeaderLines', 2);

        % utolso ures oszlop eltavolitasa, ha letezik
        if all(isnan(rawData(:, end)))
            rawData(:, end) = [];
        end

        % END OF FILE sor kiszurese
        rawData = rawData(~all(isnan(rawData), 2), :);

        % ellenorzes: 2048 mintanak kell lennie
        if size(rawData, 1) ~= 2048
            warning('%s fajlban %d sor maradt, nem 2048.', fileName, size(rawData, 1));
        end

        % idotengely
        t = (0:size(rawData,1)-1)' * dt;

        % elso oszlop kihagyasa, tobbi oszlop az egyedi gorbe
        curves = rawData(:, 2:end);
        numCurves = size(curves, 2);

        % averaging
        avgCurve = mean(curves, 2, 'omitnan');

        %% =========================
        %  Idotartomanyos vizualizacio
        % =========================

        figure(timeFig);
        subplot(numRows, numCols, subPlotIdx);
        hold on;

        % egyedi gorbek
        plot(t, curves, 'Color', [0.8 0.8 0.8]);

        % atlaggorbe
        plot(t, avgCurve, 'r', 'LineWidth', 1.5);

        title(sprintf('%s | n=%d', fileName, numCurves), ...
              'Interpreter', 'none', 'FontSize', 7);

        grid on;
        xlim([t(1) t(end)]);

        xlabel('Ido (s)');
        ylabel('V');

        %% =========================
        %  Frekvenciatartomanyos elokeszites
        % =========================

        % FFT elott NaN kezeles
        curvesClean = fillmissing(curves, 'linear', 1, 'EndValues', 'nearest');

        % DC komponens levonasa
        curvesForFFT = curvesClean - mean(curvesClean, 1, 'omitnan');

        avgCurveClean = fillmissing(avgCurve, 'linear', 'EndValues', 'nearest');
        avgCurveForFFT = avgCurveClean - mean(avgCurveClean, 'omitnan');

        N = size(curvesForFFT, 1);

        % frekvenciatengely
        f = Fs * (0:(N/2)) / N;

        % egyedi gorbek FFT-je
        Y = fft(curvesForFFT);

        P2 = abs(Y / N);
        P1 = P2(1:N/2+1, :);
        P1(2:end-1, :) = 2 * P1(2:end-1, :);

        % atlaggorbe FFT-je
        Yavg = fft(avgCurveForFFT);

        P2avg = abs(Yavg / N);
        P1avg = P2avg(1:N/2+1);
        P1avg(2:end-1) = 2 * P1avg(2:end-1);

        % 4 Hz-hez legkozelebbi frekvenciaindex
        [~, idx4Hz] = min(abs(f - targetFreq));
        amp4Hz = P1avg(idx4Hz);

        %% =========================
        %  Frekvenciatartomanyos vizualizacio
        % =========================

        figure(freqFig);
        subplot(numRows, numCols, subPlotIdx);
        hold on;

        % egyedi gorbek spektruma
        plot(f, P1, 'Color', [0.8 0.8 0.8]);

        % atlaggorbe spektruma
        plot(f, P1avg, 'r', 'LineWidth', 1.5);

        % 4 Hz jeloles
        xline(targetFreq, '--k', '4 Hz', ...
              'LabelVerticalAlignment', 'bottom', ...
              'LabelHorizontalAlignment', 'left', ...
              'FontSize', 7);

        title(sprintf('%s | A_{4Hz}=%.3g', fileName, amp4Hz), ...
              'Interpreter', 'none', 'FontSize', 7);

        grid on;
        xlim([0 maxFreqToPlot]);

        xlabel('Frekvencia (Hz)');
        ylabel('Amplitudo');

        %% =========================
        %  Mentes tablazatba
        % =========================

        newRow = table( ...
            {mouseID}, ...
            {fileName}, ...
            numCurves, ...
            {avgCurve'}, ...
            {f}, ...
            {P1avg'}, ...
            amp4Hz, ...
            'VariableNames', { ...
                'MouseId', ...
                'FileName', ...
                'NumCurves', ...
                'AverageSignal', ...
                'FrequencyVector', ...
                'AverageSpectrum', ...
                'Amp4Hz' ...
            });

        results = [results; newRow];

        subPlotIdx = subPlotIdx + 1;

        % ha veletlenul 4-nel tobb file lenne egy egerhez
        if subPlotIdx > 5
            warning('A(z) %s egerhez 4-nel tobb file tartozik. A 2x2 elrendezes megtelt.', mouseID);
        end

    catch exception
        fprintf('Hiba a fajlnal %s: %s\n', fileName, exception.message);
    end
end

%% =========================
%  Utolso eger abrainak mentese
% =========================

if ~isempty(currentMouse)
    exportgraphics(timeFig, fullfile(outputDir, [currentMouse '_time.png']), 'Resolution', 300);
    exportgraphics(freqFig, fullfile(outputDir, [currentMouse '_frequency.png']), 'Resolution', 300);

    fprintf('Abrak mentve: %s\n', currentMouse);
end

%% =========================
%  Feldolgozott adatok mentese
% =========================

save(fullfile(parentDir, 'processed_data.mat'), 'results');

disp('Kesz! Az egerenkenti 2x2-es idotartomanyos es frekvenciatartomanyos grafikonok elkeszultek es PNG formatumban mentve lettek.');