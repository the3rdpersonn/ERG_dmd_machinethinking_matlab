clear
clc

parentDir = 'D:\uni\s8\optomechatronika projekt\meresek_20';

filePattern = fullfile(parentDir, '**', '*.csv');
allFiles = dir(filePattern);

fprintf('%d fajl talalva\n', length(allFiles));

% üres feature table
featureTable = table();

for k = 1:length(allFiles)

    folderName = allFiles(k).folder;
    fileName   = allFiles(k).name;
    fullPath   = fullfile(folderName,fileName);

    % egér ID a mappa nevéből
    [~, mouseID] = fileparts(folderName);

    fprintf('Beolvasas: %s\n', fileName);

    try

        % --- CSV BEOLVASAS ---

        rawData = readmatrix(fullPath,'Delimiter',';','NumHeaderLines',2);

        % üres utolsó oszlop törlése
        if all(isnan(rawData(:,end)))
            rawData(:,end) = [];
        end

        % END OF FILE sor törlése
        rawData = rawData(~all(isnan(rawData),2),:);

        % idővektor
        dt = 0.4883/1000;
        t = (0:size(rawData,1)-1)' * dt;

        % görbék
        curves = rawData(:,2:end);

        numCurves = size(curves,2);

        % --- ÁTLAGGÖRBE ---

        avgCurve = mean(curves,2,'omitnan');

        % --- IDŐTARTOMÁNYI FEATURE-ÖK ---

        t_ms = t * 1000;

        mean_0_50 = mean(avgCurve(t_ms>=0 & t_ms<50));
        mean_50_100 = mean(avgCurve(t_ms>=50 & t_ms<100));
        mean_100_200 = mean(avgCurve(t_ms>=100 & t_ms<200));

        [peakMax, idxMax] = max(avgCurve);
        peakMaxTime = t_ms(idxMax);

        [peakMin, idxMin] = min(avgCurve);
        peakMinTime = t_ms(idxMin);

        rmsVal = rms(avgCurve);
        implicitTime = peakMaxTime;
        totalEnergy = sum(avgCurve.^2);

        % --- FREKVENCIATARTOMÁNY ---

        Fs = 1/dt;

        x = detrend(avgCurve);

        N = length(x);

        Y = fft(x);

        P2 = abs(Y/N).^2;
        P1 = P2(1:floor(N/2)+1);

        f = Fs*(0:floor(N/2))/N;

        band_0_10 = sum(P1(f>=0 & f<10));
        band_10_30 = sum(P1(f>=10 & f<30));
        band_30_80 = sum(P1(f>=30 & f<80));

        [~,idxDom] = max(P1);
        domFreq = f(idxDom);

        [~, idx4] = min(abs(f - 4));
        amp4Hz = P1(idx4);
        phase4Hz = angle(Y(idx4));

        % --- ÚJ SOR A FEATURE TABLE-BE ---

        newRow = table( ...
            string(mouseID), ...
            string(fileName), ...
            numCurves, ...
            mean_0_50, ...
            mean_50_100, ...
            mean_100_200, ...
            peakMax, ...
            peakMaxTime, ...
            peakMin, ...
            peakMinTime, ...
            rmsVal, ...
            implicitTime, ...
            totalEnergy, ...
            band_0_10, ...
            band_10_30, ...
            band_30_80, ...
            domFreq, ...
            amp4Hz, ...
            phase4Hz, ...
            'VariableNames',{ ...
            'MouseID','FileName','NumCurves', ...
            'Mean_0_50','Mean_50_100','Mean_100_200', ...
            'PeakMax','PeakMaxTime', ...
            'PeakMin','PeakMinTime', ...
            'RMS','ImplicitTime','TotalEnergy', ...
            'Band_0_10','Band_10_30','Band_30_80', ...
            'DominantFreq','Amp4Hz','Phase4Hz'});

        featureTable = [featureTable; newRow];

    catch exception

        fprintf('Hiba: %s\n',exception.message);

    end

end

save(fullfile(parentDir,'feature_table.mat'),'featureTable')

disp('Feature table kesz.')

% --- FEATURE MATRIX LETREHOZASA ---

% numerikus feature-ök kiválasztása
X = featureTable{:,4:end};

% normalizálás
X = zscore(X);

% mentés későbbi AI lépésekhez
save(fullfile(parentDir,'feature_matrix.mat'),'X','featureTable')

disp('Feature matrix kesz.')

size(X)

disp('Uj feature-ok ellenorzese:')
disp(featureTable(:, {'MouseID','FileName','ImplicitTime','TotalEnergy','Amp4Hz','Phase4Hz'}))

% --- PCA ---

[coeff, score, latent, ~, explained] = pca(X);

figure
scatter(score(:,1), score(:,2), 60, 'filled')
xlabel(['PC1 (' num2str(explained(1), '%.1f') '%)'])
ylabel(['PC2 (' num2str(explained(2), '%.1f') '%)'])
title('PCA of ERG feature matrix')
grid on

figure
tiledlayout(2,2)

nexttile
histogram(featureTable.ImplicitTime)
xlabel('ms')
title('ImplicitTime')

nexttile
histogram(featureTable.TotalEnergy)
xlabel('Energy')
title('TotalEnergy')

nexttile
histogram(featureTable.Amp4Hz)
xlabel('4 Hz amplitude')
title('Amp4Hz')

nexttile
histogram(featureTable.RMS)
xlabel('RMS')
title('RMS')