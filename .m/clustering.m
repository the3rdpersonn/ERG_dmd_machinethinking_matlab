clear
clc

parentDir = 'D:\uni\s8\optomechatronika projekt\meresek';

load(fullfile(parentDir,'feature_matrix.mat'),'X','featureTable')

% --- PCA a vizualizaciohoz ---
[coeff, score, latent, ~, explained] = pca(X);

% --- KLASZTEREZES: K-MEANS ---
k = 6;   % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rng(1)   % reprodukalhatosag
idx = kmeans(X, k, 'Replicates', 20);

% klaszterazonosito hozzadasa a tablazathoz
featureTable.ClusterID = idx;

% --- PCA abra klaszterekkel ---
figure
gscatter(score(:,1), score(:,2), idx)
xlabel(['PC1 (' num2str(explained(1), '%.1f') '%)'])
ylabel(['PC2 (' num2str(explained(2), '%.1f') '%)'])
title(['K-means clustering, k = ' num2str(k)])
grid on

% --- klaszteronkenti elemszam ---
disp('Clusterek elemszamai:')
disp(groupcounts(featureTable.ClusterID))

% --- klaszteronkenti fajlok listazasa ---
for c = 1:k
    fprintf('\n===== Cluster %d =====\n', c);
    disp(featureTable(featureTable.ClusterID == c, {'MouseID','FileName'}))
end

% --- mentes ---
save(fullfile(parentDir,'clustering_results.mat'), ...
    'featureTable', 'idx', 'k', 'score', 'explained')

disp('Klaszterezes kesz.')

% --- KLASZTERENKENTI ATLAGGORBEK KIRAJZOLASA ---

dt = 0.4883 / 1000;   % s

for c = 1:k
    clusterRows = featureTable(featureTable.ClusterID == c, :);

    if isempty(clusterRows)
        continue
    end

    figure('Name', ['Cluster ' num2str(c)], ...
           'NumberTitle', 'off', ...
           'Units', 'normalized', ...
           'Position', [0.1 0.1 0.8 0.7]);

    hold on

    clusterCurves = [];

    for i = 1:height(clusterRows)
        mouseID = char(clusterRows.MouseID(i));
        fileName = char(clusterRows.FileName(i));

        fullPath = fullfile(parentDir, mouseID, fileName);

        rawData = readmatrix(fullPath, 'Delimiter', ';', 'NumHeaderLines', 2);

        % üres utolsó oszlop törlése
        if all(isnan(rawData(:,end)))
            rawData(:,end) = [];
        end

        % teljesen NaN sorok törlése
        rawData = rawData(~all(isnan(rawData),2), :);

        curves = rawData(:,2:end);
        avgCurve = mean(curves, 2, 'omitnan');

        t = (0:length(avgCurve)-1)' * dt;

        % egyedi átlaggörbék halványan
        plot(t, avgCurve, 'Color', [0.75 0.75 0.75], 'LineWidth', 1);

        clusterCurves(:, end+1) = avgCurve;
    end

    % klaszter összátlag
    clusterMean = mean(clusterCurves, 2, 'omitnan');
    plot(t, clusterMean, 'r', 'LineWidth', 2.5);

    xlabel('Ido (s)');
    ylabel('V');
    title(['Cluster ' num2str(c) ' atlaggorbei | n = ' num2str(height(clusterRows))]);
    grid on
    hold off
end