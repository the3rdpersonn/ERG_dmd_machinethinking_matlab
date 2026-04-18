% alapadatok, eleresi utak beallitasa
parentDir = 'D:\uni\s8\optomechatronika projekt\meresek';
filePattern = fullfile(parentDir, '**', '*.csv');
allFiles = dir(filePattern);

% ures tablazat az eredmenyeknek
results = table();

fprintf('%d db fajl talalhato osszesen.\n', length(allFiles));

currentMouse = '';
subPlotIdx = 0;

for k = 1:length(allFiles)
    folderName = allFiles(k).folder;
    fileName = allFiles(k).name;
    fullFilePath = fullfile(folderName, fileName);

    % egerek azonositoja a mappa nevebol
    [~, mouseID] = fileparts(folderName);

    % egerenkenti diagramok
    if ~strcmp(mouseID, currentMouse)
        currentMouse = mouseID;

        mouseFiles = dir(fullfile(folderName, '*.csv'));
        numFilesForThisMouse = length(mouseFiles);

        numCols = 3;
        numRows = ceil(numFilesForThisMouse / numCols);

        figure('Name', ['Eger: ' mouseID], ...
               'NumberTitle', 'off', ...
               'Units', 'normalized', ...
               'Position', [0.05 0.05 0.9 0.85]);

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
        % ez a sor teljesen NaN lesz readmatrix utan
        rawData = rawData(~all(isnan(rawData), 2), :);

        % ellenorzes: 2048 mintanak kell lennie
        if size(rawData, 1) ~= 2048
            warning('%s fajlban %d sor maradt, nem 2048.', fileName, size(rawData, 1));
        end

        % idotengely a fajl hossza alapjan
        dt = 0.4883 / 1000;   % s
        t = (0:size(rawData,1)-1)' * dt;

        % gorbek: 120 oszlop->valahol nem 120, hanem 80 vagy 40 oszlop
        % (curve) van
        curves = rawData(:, 2:end);

        numCurves = size(curves, 2);

        % averaging
        avgCurve = mean(curves, 2, 'omitnan');

        % vizualizacio
        subplot(numRows, numCols, subPlotIdx);
        hold on;

        % kulon gorbek
        plot(t, curves, 'Color', [0.8 0.8 0.8]);

        % atlaggörbe
        plot(t, avgCurve, 'r', 'LineWidth', 1.5);

        title(sprintf('%s | n=%d', fileName, numCurves), 'Interpreter', 'none', 'FontSize', 7);
        grid on;
        xlim([t(1) t(end)]);

        if subPlotIdx > (numRows - 1) * numCols
            xlabel('Ido (s)');
        end
        if mod(subPlotIdx, numCols) == 1
            ylabel('V');
        end

        subPlotIdx = subPlotIdx + 1;

        % mentes tablazatba
        newRow = table({mouseID}, {fileName}, {avgCurve'}, ...
            'VariableNames', {'MouseId', 'FileName', 'AverageSignal'});
        results = [results; newRow];

    catch exception
        fprintf('Hiba a fajlnal %s: %s\n', fileName, exception.message);
    end
end

save(fullfile(parentDir, 'processed_data.mat'), 'results');
disp('Kesz! A grafikonok elkeszultek.');

%kiegesziteni frekvenciatartomanyos grafikonrajzolassal