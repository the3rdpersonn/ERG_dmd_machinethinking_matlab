clear
clc

parentDir = 'D:\uni\s8\optomechatronika projekt\meresek_20';

load(fullfile(parentDir, 'feature_matrix.mat'), 'featureTable')

% -------------------------------------------------------------
% FEATURE VALASZTAS
% sawtooth core + altalanos tamogato feature-ok
% flicker/H1 kiveve a clusteringbol
% -------------------------------------------------------------
featureNamesForClustering = { ...
    'NumCurves', ...
    'BaselineStd', ...
    'NONAmp','NONTime', ...
    'PONAmp','PONTime', ...
    'LNONAmp','LNONTime', ...
    'NOFFAmp','NOFFTime', ...
    'POFFAmp','POFFTime', ...
    'AWaveAmp','AWaveTime', ...
    'BWaveAmp','BWaveTime','BARatio', ...
    'OPRMS','OPPeakToPeak', ...
    'Mean_0_50','Mean_50_100','Mean_100_200', ...
    'RMS','TotalEnergy', ...
    'Band_0_10','Band_10_30','Band_30_80'};

Xraw = featureTable{:, featureNamesForClustering};

% -------------------------------------------------------------
% Z-SCORE NaN-TURESSel
% -------------------------------------------------------------
X = nan(size(Xraw));
mu = nan(1, size(Xraw,2));
sigma = nan(1, size(Xraw,2));

for j = 1:size(Xraw,2)
    col = Xraw(:,j);
    mu(j) = mean(col, 'omitnan');
    sigma(j) = std(col, 0, 'omitnan');

    if ~isnan(sigma(j)) && sigma(j) > 0
        X(:,j) = (col - mu(j)) / sigma(j);
    else
        X(:,j) = NaN(size(col));
    end
end

validCols = ~all(isnan(X),1) & (sigma > 0);
X = X(:, validCols);
featureNamesForClustering = featureNamesForClustering(validCols); %????????????

% -------------------------------------------------------------
% FEATURE WEIGHTS
% sawtooth core eros
% a/b-wave + baseline + OP kozepes
% shape/frequency gyengebb
% -------------------------------------------------------------
featureWeights = ones(1, numel(featureNamesForClustering));

for j = 1:numel(featureNamesForClustering)
    name = featureNamesForClustering{j};

    if ismember(name, {'NONAmp','NONTime','PONAmp','PONTime','LNONAmp','LNONTime','NOFFAmp','NOFFTime','POFFAmp','POFFTime'})
        featureWeights(j) = 1.5;
    elseif ismember(name, {'AWaveAmp','AWaveTime','BWaveAmp','BWaveTime','BARatio','BaselineStd','OPRMS','OPPeakToPeak'})
        featureWeights(j) = 1;
    elseif ismember(name, {'Mean_0_50','Mean_50_100','Mean_100_200','RMS','TotalEnergy','Band_0_10','Band_10_30','Band_30_80'})
        featureWeights(j) = 0.8;
    else
        featureWeights(j) = 1;
    end
end

% -------------------------------------------------------------
% PCA CSAK VIZUALIZACIORA
% -------------------------------------------------------------
Xviz = X;
for j = 1:size(Xviz,2)
    col = Xviz(:,j);
    med = median(col, 'omitnan');
    if isnan(med)
        med = 0;
    end
    col(isnan(col)) = med;
    Xviz(:,j) = col;
end

[coeff, score, latent, ~, explained] = pca(Xviz);

% -------------------------------------------------------------
% MASZKOLT K-MEANS
% -------------------------------------------------------------
k = 4;
maxAttempts = 40;
maxIter = 100;

bestIdx = [];
bestCentroids = [];
bestObjective = inf;

for attempt = 1:maxAttempts
    rng(attempt)

    [tmpIdx, tmpC, tmpObj] = masked_kmeans(X, k, featureWeights, maxIter);

    clusterSizes = groupcounts(tmpIdx);

    penalty = 0;
    if any(clusterSizes <= 1)
        penalty = 1000 * sum(clusterSizes <= 1);
    end

    scoreAttempt = tmpObj + penalty;

    if scoreAttempt < bestObjective
        bestObjective = scoreAttempt;
        bestIdx = tmpIdx;
        bestCentroids = tmpC;
    end
end

idx = bestIdx;
centroids = bestCentroids;
featureTable.ClusterID = idx;

clusterSizes = groupcounts(idx);

if any(clusterSizes == 1)
    fprintf('\nFigyelem: %d db 1 elemu cluster maradt. Valoszinuleg a k tul nagy.\n', sum(clusterSizes == 1));
end

% -------------------------------------------------------------
% PCA abra klaszterekkel
% -------------------------------------------------------------
figure
gscatter(score(:,1), score(:,2), idx)
xlabel(['PC1 (' num2str(explained(1), '%.1f') '%)'])
ylabel(['PC2 (' num2str(explained(2), '%.1f') '%)'])
title(['Masked k-means clustering, k = ' num2str(k)])
grid on

% -------------------------------------------------------------
% pseudo-silhouette maszkolt tavolsaggal
% -------------------------------------------------------------
silhValues = masked_silhouette(X, idx, featureWeights);
meanSilh = mean(silhValues, 'omitnan');

figure
hold on

colors = lines(k);
yStart = 0;

for c = 1:k
    vals = silhValues(idx == c);
    vals = vals(~isnan(vals));
    vals = sort(vals, 'descend');

    if isempty(vals)
        continue
    end

    y = yStart + (1:numel(vals));
    barh(y, vals, 1.0, 'FaceColor', colors(c,:), 'EdgeColor', 'none');

    % klaszterszam kiirasa bal oldalra
    text(-0.08, mean(y), ['C' num2str(c)], 'FontWeight', 'bold');

    yStart = y(end) + 2;
end

xline(meanSilh, '--r', ['mean = ' num2str(meanSilh, '%.3f')], 'LineWidth', 1.5);
xline(0, '-k');

xlabel('Masked silhouette value')
ylabel('Mintak (csoportositva klaszterenkent)')
title(['Masked silhouette diagram | k = ' num2str(k)])
xlim([-1 1])
grid on
hold off

fprintf('\nAtlagos masked silhouette score: %.3f\n', meanSilh);

% -------------------------------------------------------------
% klaszteronkenti elemszam
% -------------------------------------------------------------
disp('Clusterek elemszamai:')
disp(groupcounts(featureTable.ClusterID))

% -------------------------------------------------------------
% klaszteronkenti fajlok listazasa
% -------------------------------------------------------------
for c = 1:k
    fprintf('\n===== Cluster %d =====\n', c);
    disp(featureTable(featureTable.ClusterID == c, {'MouseID','FileName'}))
end

% -------------------------------------------------------------
% CLUSTER FEATURE SUMMARY
% -------------------------------------------------------------
summaryVars = intersect( ...
    {'NONAmp','NONTime','PONAmp','PONTime','LNONAmp','LNONTime','NOFFAmp','NOFFTime','POFFAmp','POFFTime', ...
     'AWaveAmp','AWaveTime','BWaveAmp','BWaveTime','BARatio','BaselineStd','OPRMS','RMS','TotalEnergy'}, ...
    featureTable.Properties.VariableNames);

summaryTable = groupsummary(featureTable, 'ClusterID', 'mean', summaryVars);

disp(' ')
disp('===== Cluster feature summary =====')
disp(summaryTable)

% -------------------------------------------------------------
% FEATURE DISTRIBUTION PLOTOK
% -------------------------------------------------------------
figure
tiledlayout(2,3)

plotBoxIfExists(featureTable, 'NONAmp', 'Cluster', '\muV', 'NON amplitude')
plotBoxIfExists(featureTable, 'PONAmp', 'Cluster', '\muV', 'PON amplitude')
plotBoxIfExists(featureTable, 'LNONAmp', 'Cluster', '\muV', 'LNON amplitude')
plotBoxIfExists(featureTable, 'NOFFAmp', 'Cluster', '\muV', 'NOFF amplitude')
plotBoxIfExists(featureTable, 'POFFAmp', 'Cluster', '\muV', 'POFF amplitude')
plotBoxIfExists(featureTable, 'BARatio', 'Cluster', 'ratio', 'B/A ratio')

% -------------------------------------------------------------
% mentes
% -------------------------------------------------------------
save(fullfile(parentDir, 'clustering_results.mat'), ...
    'featureTable', 'idx', 'k', 'score', 'explained', 'summaryTable', ...
    'meanSilh', 'featureNamesForClustering', 'X', 'centroids', 'featureWeights')

disp('Klaszterezes kesz.')

% -------------------------------------------------------------
% klaszterenkenti atlaggorbek
% -------------------------------------------------------------
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

% =============================================================
% HELPER FUNCTIONS
% =============================================================

function plotBoxIfExists(featureTable, varName, xlab, ylab, ttl)
    nexttile
    if ismember(varName, featureTable.Properties.VariableNames)
        boxchart(categorical(featureTable.ClusterID), featureTable.(varName))
        xlabel(xlab)
        ylabel(ylab)
        title(ttl)
    else
        axis off
        title([varName ' hianyzik'])
    end
end

function [idx, centroids, totalObj] = masked_kmeans(X, k, featureWeights, maxIter)
    n = size(X,1);

    validRows = find(any(~isnan(X), 2));
    if numel(validRows) < k
        error('Tul sok az ures sor / NaN.');
    end

    seedRows = validRows(randperm(numel(validRows), k));
    centroids = X(seedRows, :);

    idx = ones(n,1);

    for iter = 1:maxIter
        oldIdx = idx;

        for i = 1:n
            d = nan(k,1);
            for c = 1:k
                d(c) = masked_distance(X(i,:), centroids(c,:), featureWeights);
            end

            if all(isnan(d))
                idx(i) = randi(k);
            else
                [~, idx(i)] = min(d);
            end
        end

        for c = 1:k
            members = X(idx == c, :);

            if isempty(members)
                ridx = validRows(randi(numel(validRows)));
                centroids(c,:) = X(ridx,:);
            else
                centroids(c,:) = mean(members, 1, 'omitnan');
            end
        end

        if isequal(idx, oldIdx)
            break
        end
    end

    totalObj = 0;
    for i = 1:n
        di = masked_distance(X(i,:), centroids(idx(i),:), featureWeights);
        if ~isnan(di)
            totalObj = totalObj + di^2;
        end
    end
end

function d = masked_distance(x, y, w)
    valid = ~isnan(x) & ~isnan(y) & ~isnan(w);

    if ~any(valid)
        d = NaN;
        return
    end

    diff2 = w(valid) .* (x(valid) - y(valid)).^2;
    d = sqrt(sum(diff2) / sum(w(valid)));
end

function silh = masked_silhouette(X, idx, featureWeights)
    n = size(X,1);
    silh = nan(n,1);
    clusters = unique(idx(:))';

    for i = 1:n
        same = find(idx == idx(i));
        same(same == i) = [];

        if isempty(same)
            silh(i) = NaN;
            continue
        end

        aVals = nan(numel(same),1);
        for j = 1:numel(same)
            aVals(j) = masked_distance(X(i,:), X(same(j),:), featureWeights);
        end
        a = mean(aVals, 'omitnan');

        b = inf;
        for c = clusters
            if c == idx(i)
                continue
            end

            other = find(idx == c);
            if isempty(other)
                continue
            end

            bVals = nan(numel(other),1);
            for j = 1:numel(other)
                bVals(j) = masked_distance(X(i,:), X(other(j),:), featureWeights);
            end

            bc = mean(bVals, 'omitnan');
            if ~isnan(bc)
                b = min(b, bc);
            end
        end

        if isinf(b) || isnan(a) || isnan(b) || max(a,b) == 0
            silh(i) = NaN;
        else
            silh(i) = (b - a) / max(a, b);
        end
    end
end