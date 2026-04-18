clear
clc

parentDir = 'D:\uni\s8\optomechatronika projekt\meresek_20';

load(fullfile(parentDir,'feature_matrix.mat'),'featureTable')

featureNamesForClustering = { ...
    'NumCurves', ...
    'Mean_0_50','Mean_50_100','Mean_100_200', ...
    'PeakMax','PeakMin','PeakMinTime', ...
    'RMS', ...
    'Band_0_10','Band_10_30','Band_30_80', ...
    'DominantFreq','Amp4Hz','Phase4Hz'};

X = featureTable{:, featureNamesForClustering};
X = zscore(X);

colStd = std(X, 0, 1, 'omitnan');
validCols = colStd > 0 & all(isfinite(X), 1);

X = X(:, validCols);
featureNamesForClustering = featureNamesForClustering(validCols);

validRows = all(isfinite(X), 2);
X = X(validRows, :);
featureTable = featureTable(validRows, :);

% --- PCA a vizualizaciohoz ---
[coeff, score, latent, ~, explained] = pca(X);

% --- KLASZTEREZES: K-MEANS ---
k = 4; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bestIdx = [];
bestSumd = inf;
maxAttempts = 50;

for attempt = 1:maxAttempts
    rng(attempt)
    [tmpIdx, ~, sumd] = kmeans(X, k, 'Replicates', 20);

    clusterSizes = groupcounts(tmpIdx);

    if all(clusterSizes > 1) %garancia hogy nem lesz 1 elemu cluster, cserebe limitalja a machine thinkinget
        bestIdx = tmpIdx;
        bestSumd = sum(sumd);
        break
    end

    if sum(sumd) < bestSumd
        bestIdx = tmpIdx;
        bestSumd = sum(sumd);
    end
end

idx = bestIdx;
featureTable.ClusterID = idx;

clusterSizes = groupcounts(idx);

if any(clusterSizes == 1)
    fprintf('\nFigyelem: %d db 1 elemu cluster maradt. Valoszinuleg a k tul nagy.\n', sum(clusterSizes == 1));
end
% klaszterazonosito hozzadasa a tablazathoz
featureTable.ClusterID = idx;

% --- PCA abra klaszterekkel ---
figure
gscatter(score(:,1), score(:,2), idx)
xlabel(['PC1 (' num2str(explained(1), '%.1f') '%)'])
ylabel(['PC2 (' num2str(explained(2), '%.1f') '%)'])
title(['K-means clustering, k = ' num2str(k)])
grid on

% --- SILHOUETTE SCORE ---
clusterSizes = groupcounts(idx);
validClusterIDs = find(clusterSizes > 1);
validRows = ismember(idx, validClusterIDs);

if numel(validClusterIDs) >= 2
    figure
    [silhValues, h] = silhouette(X(validRows,:), idx(validRows));
    title(['Silhouette plot | k = ' num2str(k) ])

    meanSilh = mean(silhValues, 'omitnan');
    fprintf('\nAtlagos silhouette score: %.3f\n', meanSilh);
else
    meanSilh = NaN;
    fprintf('\nA silhouette score nem ertekelheto: nincs eleg nem-singleton cluster.\n');
end

% --- klaszteronkenti elemszam ---
disp('Clusterek elemszamai:')
disp(groupcounts(featureTable.ClusterID))

% --- klaszteronkenti fajlok listazasa ---
for c = 1:k
    fprintf('\n===== Cluster %d =====\n', c);
    disp(featureTable(featureTable.ClusterID == c, {'MouseID','FileName'}))
end

% --- CLUSTER FEATURE SUMMARY ---
summaryTable = groupsummary(featureTable, 'ClusterID', 'mean', ...
    {'ImplicitTime','TotalEnergy','Amp4Hz','RMS','PeakMax','PeakMin'});

disp(' ')
disp('===== Cluster feature summary =====')
disp(summaryTable)

% --- FEATURE DISTRIBUTION PLOTOK ---
figure
tiledlayout(2,2)

nexttile
boxchart(categorical(featureTable.ClusterID), featureTable.ImplicitTime)
xlabel('Cluster')
ylabel('ms')
title('Implicit time')

nexttile
boxchart(categorical(featureTable.ClusterID), featureTable.TotalEnergy)
xlabel('Cluster')
ylabel('Energy')
title('Signal energy')

nexttile
boxchart(categorical(featureTable.ClusterID), featureTable.Amp4Hz)
xlabel('Cluster')
ylabel('4 Hz amplitude')
title('4 Hz response')

nexttile
boxchart(categorical(featureTable.ClusterID), featureTable.RMS)
xlabel('Cluster')
ylabel('RMS')
title('Signal RMS')

% --- mentes ---
save(fullfile(parentDir,'clustering_results.mat'), ...
    'featureTable', 'idx', 'k', 'score', 'explained', 'summaryTable', ...
    'meanSilh', 'featureNamesForClustering', 'X')

disp('Klaszterezes kesz.')

% --- KLASZTERENKENTI ATLAGGORBEK EGY FIGURE-BAN ---
dt = 0.4883 / 1000;

figure('Name', 'Cluster atlaggorbek', ...
       'NumberTitle', 'off', ...
       'Units', 'normalized', ...
       'Position', [0.05 0.05 0.9 0.85]);

numCols = 2;
numRows = ceil(k / numCols);

for c = 1:k
    clusterRows = featureTable(featureTable.ClusterID == c, :);

    if isempty(clusterRows)
        continue
    end

    subplot(numRows, numCols, c)
    hold on

    clusterCurves = [];

    for i = 1:height(clusterRows)
        mouseID = char(clusterRows.MouseID(i));
        fileName = char(clusterRows.FileName(i));

        fullPath = fullfile(parentDir, mouseID, fileName);

        rawData = readmatrix(fullPath, 'Delimiter', ';', 'NumHeaderLines', 2);

        if all(isnan(rawData(:,end)))
            rawData(:,end) = [];
        end

        rawData = rawData(~all(isnan(rawData),2), :);

        curves = rawData(:,2:end);
        avgCurve = mean(curves, 2, 'omitnan');

        t = (0:length(avgCurve)-1)' * dt;

        plot(t, avgCurve, 'Color', [0.75 0.75 0.75], 'LineWidth', 1);

        clusterCurves(:, end+1) = avgCurve;
    end

    clusterMean = mean(clusterCurves, 2, 'omitnan');
    plot(t, clusterMean, 'r', 'LineWidth', 2.5);

    xlabel('Ido (s)');
    ylabel('V');
    title(['Cluster ' num2str(c) ' | n = ' num2str(height(clusterRows))]);
    grid on
    hold off
end

% --- KLASZTER FREKVENCIA SPEKTRUM ---
Fs = 1/dt;

for c = 1:k
    clusterRows = featureTable(featureTable.ClusterID == c, :);

    if isempty(clusterRows)
        continue
    end

    figure('Name',['Cluster ' num2str(c) ' spectrum'], ...
           'NumberTitle', 'off')

    hold on
    spectra = [];

    for i = 1:height(clusterRows)
        mouseID = char(clusterRows.MouseID(i));
        fileName = char(clusterRows.FileName(i));

        fullPath = fullfile(parentDir, mouseID, fileName);

        rawData = readmatrix(fullPath,'Delimiter',';','NumHeaderLines',2);

        if all(isnan(rawData(:,end)))
            rawData(:,end) = [];
        end

        rawData = rawData(~all(isnan(rawData),2),:);

        curves = rawData(:,2:end);
        avgCurve = mean(curves,2,'omitnan');

        x = detrend(avgCurve);
        N = length(x);

        Y = fft(x);
        P2 = abs(Y/N).^2;
        P1 = P2(1:floor(N/2)+1);

        f = Fs*(0:floor(N/2))/N;

        plot(f, P1, 'Color', [0.8 0.8 0.8])

        spectra(:, end+1) = P1;
    end

    meanSpectrum = mean(spectra, 2, 'omitnan');
    plot(f, meanSpectrum, 'r', 'LineWidth', 2)

    xline(4, '--')
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title(['Cluster ' num2str(c) ' frequency spectrum'])
    xlim([0 40])
    grid on
    hold off
end