clear
clc

parentDir = 'D:\uni\s8\optomechatronika projekt\meresek_20';

filePattern = fullfile(parentDir, '**', '*.csv');
allFiles = dir(filePattern);

fprintf('%d fajl talalva\n', length(allFiles));

dt = 0.4883 / 1000;
Fs = 1 / dt;

featureTable = table();

for k = 1:length(allFiles)

    folderName = allFiles(k).folder;
    fileName   = allFiles(k).name;
    fullPath   = fullfile(folderName, fileName);

    [~, mouseID] = fileparts(folderName);

    fprintf('Beolvasas: %s\n', fileName);

    try
        % -------------------------------------------------------------
        % CSV BEOLVASAS
        % -------------------------------------------------------------
        rawData = readmatrix(fullPath, 'Delimiter', ';', 'NumHeaderLines', 2);

        if isempty(rawData)
            error('Ures fajl vagy olvashatatlan adat.');
        end

        if all(isnan(rawData(:,end)))
            rawData(:,end) = [];
        end

        rawData = rawData(~all(isnan(rawData),2), :);

        if size(rawData,2) < 2
            error('Nincs eleg oszlop a gorbekhez.');
        end

        t = (0:size(rawData,1)-1)' * dt;
        t_ms = t * 1000;

        curves = rawData(:,2:end);
        numCurves = size(curves, 2);

        avgCurve = mean(curves, 2, 'omitnan');

        % -------------------------------------------------------------
        % 0) BASELINE / QUALITY
        % -------------------------------------------------------------
        baseMask = t_ms >= 0 & t_ms < 10;
        if any(baseMask)
            baselineMean = mean(avgCurve(baseMask), 'omitnan');
            baselineStd  = std(avgCurve(baseMask), 0, 'omitnan');
        else
            baselineMean = NaN;
            baselineStd  = NaN;
        end

        x0 = avgCurve - baselineMean;

        if any(isnan(x0))
            x0 = fillmissing(x0, 'linear', 'EndValues', 'nearest');
        end

        minAmp = 1e-7;

        % -------------------------------------------------------------
        % 1) ALTALANOS SHAPE FEATURE-OK (meghagyva)
        % -------------------------------------------------------------
        mean_0_50    = mean(x0(t_ms >= 0   & t_ms < 50),  'omitnan');
        mean_50_100  = mean(x0(t_ms >= 50  & t_ms < 100), 'omitnan');
        mean_100_200 = mean(x0(t_ms >= 100 & t_ms < 200), 'omitnan');

        [peakMax, idxMax] = max(x0);
        peakMaxTime = t_ms(idxMax);

        [peakMin, idxMin] = min(x0);
        peakMinTime = t_ms(idxMin);

        rmsVal = rms(x0);
        totalEnergy = sum(x0.^2, 'omitnan');
        peakToPeak = peakMax - peakMin;

        % -------------------------------------------------------------
        % 2) SAWTOOTH-SPECIFIKUS FEATURE-OK (fo blokk)
        % Tsai: sawtooth low-pass 40 Hz; NON/NOFF baseline-tol,
        % PON/LNON/POFF peak-to-peak az elozo troughhoz kepest.
        % -------------------------------------------------------------
        sawCurve = lowpass_local(x0, Fs, 40);

        NONAmp  = NaN; NONTime  = NaN;
        PONAmp  = NaN; PONTime  = NaN;
        LNONAmp = NaN; LNONTime = NaN;
        NOFFAmp = NaN; NOFFTime = NaN;
        POFFAmp = NaN; POFFTime = NaN;

        period_ms = 250;          % 4 Hz sawtooth
        offChange_ms = 125;       % kb. periodus fele

        % --- ON baseline: rapid valtozas utani elso 5 ms
        onBaseMask = t_ms >= 0 & t_ms < 5;
        if any(onBaseMask)
            onBaseline = mean(sawCurve(onBaseMask), 'omitnan');
        else
            onBaseline = mean(sawCurve(1:min(10,end)), 'omitnan');
        end

        % NON
        nonMask = t_ms >= 5 & t_ms <= 80;
        if any(nonMask)
            nonSegment = sawCurve(nonMask);
            nonIdxAll = find(nonMask);

            [nonVal, nonLocal] = min(nonSegment);
            nonIdx = nonIdxAll(nonLocal);

            NONAmp  = onBaseline - nonVal;
            NONTime = t_ms(nonIdx);

            if NONAmp < minAmp || NONAmp < 0
                NONAmp = NaN;
                NONTime = NaN;
            end

            % PON
            if ~isnan(NONTime)
                ponMask = t_ms >= NONTime & t_ms <= 140;
                if any(ponMask)
                    ponSegment = sawCurve(ponMask);
                    ponIdxAll = find(ponMask);

                    [pks, locs] = findpeaks(ponSegment);

                    if ~isempty(pks)
                        [ponVal, bestPk] = max(pks);
                        ponIdx = ponIdxAll(locs(bestPk));

                        PONAmp  = ponVal - nonVal;
                        PONTime = t_ms(ponIdx);

                        if PONAmp < minAmp || PONAmp < 0
                            PONAmp = NaN;
                            PONTime = NaN;
                        end

                        % LNON
                        if ~isnan(PONTime)
                            lnonMask = t_ms >= PONTime & t_ms <= 240;
                            if any(lnonMask)
                                lnonSegment = sawCurve(lnonMask);
                                lnonIdxAll = find(lnonMask);

                                [lnonVal, lnonLocal] = min(lnonSegment);
                                lnonIdx = lnonIdxAll(lnonLocal);

                                LNONAmp  = ponVal - lnonVal;
                                LNONTime = t_ms(lnonIdx);

                                if LNONAmp < minAmp || LNONAmp < 0
                                    LNONAmp = NaN;
                                    LNONTime = NaN;
                                end
                            end
                        end
                    end
                end
            end
        end

        % --- OFF baseline: masodik rapid valtozas utani elso 5 ms
        offBaseMask = t_ms >= offChange_ms & t_ms < offChange_ms + 5;
        if any(offBaseMask)
            offBaseline = mean(sawCurve(offBaseMask), 'omitnan');
        else
            offBaseline = NaN;
        end

        % NOFF
        noffMask = t_ms >= offChange_ms + 5 & t_ms <= offChange_ms + 90;
        if any(noffMask)
            noffSegment = sawCurve(noffMask);
            noffIdxAll = find(noffMask);

            [noffVal, noffLocal] = min(noffSegment);
            noffIdx = noffIdxAll(noffLocal);

            if ~isnan(offBaseline)
                NOFFAmp = offBaseline - noffVal;
            end
            NOFFTime = t_ms(noffIdx) - offChange_ms;   % offsettol mert ido

            if NOFFAmp < minAmp || NOFFAmp < 0
                NOFFAmp = NaN;
                NOFFTime = NaN;
            end

            % POFF
            if ~isnan(NOFFTime)
                poffMask = t_ms >= t_ms(noffIdx) & t_ms <= period_ms;
                if any(poffMask)
                    poffSegment = sawCurve(poffMask);
                    poffIdxAll = find(poffMask);

                    [pks, locs] = findpeaks(poffSegment);

                    if ~isempty(pks)
                        [poffVal, bestPk] = max(pks);
                        poffIdx = poffIdxAll(locs(bestPk));

                        POFFAmp  = poffVal - noffVal;
                        POFFTime = t_ms(poffIdx) - offChange_ms;

                        if POFFAmp < minAmp || POFFAmp < 0
                            POFFAmp = NaN;
                            POFFTime = NaN;
                        end
                    end
                end
            end
        end

        % -------------------------------------------------------------
        % 3) NEM SAW-SPECIFIKUS, DE INFORMATIV FEATURE-OK
        % (bent maradnak, de nem ez lesz a fo resz)
        % -------------------------------------------------------------
        aWaveAmp = NaN;
        aWaveTime = NaN;
        aWaveSlope = NaN;
        aIdx = NaN;

        aMask = t_ms >= 0 & t_ms <= 50;
        if any(aMask)
            [aVal, aLocal] = min(x0(aMask));
            aIdxCandidates = find(aMask);
            aIdx = aIdxCandidates(aLocal);

            aWaveAmp = -aVal;
            aWaveTime = t_ms(aIdx);

            if aWaveAmp < minAmp || aWaveAmp < 0
                aWaveAmp = NaN;
                aWaveTime = NaN;
            end

            if ~isnan(aWaveTime)
                leadMask = t_ms >= 0 & t_ms <= aWaveTime;
                if nnz(leadMask) >= 3
                    pLead = polyfit(t_ms(leadMask), x0(leadMask), 1);
                    aWaveSlope = pLead(1);
                end
            end
        end

        bWaveAmp = NaN;
        bWaveTime = NaN;
        b_a_ratio = NaN;

        if ~isnan(aWaveTime)
            bMask = t_ms >= aWaveTime & t_ms <= 120;
            if any(bMask)
                [bVal, bLocal] = max(x0(bMask));
                bIdxCandidates = find(bMask);
                bIdx = bIdxCandidates(bLocal);

                bWaveAmp = bVal - x0(aIdx);
                bWaveTime = t_ms(bIdx);

                if bWaveAmp < minAmp || bWaveAmp < 0
                    bWaveAmp = NaN;
                    bWaveTime = NaN;
                end

                if ~isnan(aWaveAmp) && ~isnan(bWaveAmp) && abs(aWaveAmp) > eps
                    b_a_ratio = bWaveAmp / aWaveAmp;
                end
            end
        end

        % OP feature-ok bent maradnak
        opBand = bandpass_local(x0, Fs, 50, 300);
        opMask = t_ms >= 10 & t_ms <= 120;
        if any(opMask) && ~all(isnan(opBand(opMask)))
            opRMS = rms(opBand(opMask));
            opPeakToPeak = max(opBand(opMask)) - min(opBand(opMask));
        else
            opRMS = NaN;
            opPeakToPeak = NaN;
        end

        % -------------------------------------------------------------
        % 4) FREKVENCIATARTOMANY (shape / masodlagos)
        % Flicker-specifikus H1 marad a tablaban, de clusteringbol kimegy
        % -------------------------------------------------------------
        xDetr = detrend(x0);
        N = length(xDetr);

        Y = fft(xDetr);
        A = abs(Y / N);
        A1 = A(1:floor(N/2)+1);
        if numel(A1) > 2
            A1(2:end-1) = 2 * A1(2:end-1);
        end

        P2 = abs(Y / N).^2;
        P1 = P2(1:floor(N/2)+1);

        f = Fs * (0:floor(N/2)) / N;

        band_0_10  = sum(P1(f >= 0  & f < 10), 'omitnan');
        band_10_30 = sum(P1(f >= 10 & f < 30), 'omitnan');
        band_30_80 = sum(P1(f >= 30 & f < 80), 'omitnan');

        [~, idxDom] = max(P1);
        domFreq = f(idxDom);

        [~, idx4] = min(abs(f - 4));
        h1AmpRaw_4Hz = A1(idx4);
        h1Phase_4Hz = angle(Y(idx4));

        noiseMask4 = ((f >= 3 & f < 4) | (f > 4 & f <= 5));
        if any(noiseMask4)
            noiseAmp4 = mean(A1(noiseMask4), 'omitnan');
        else
            noiseAmp4 = NaN;
        end

        h1AmpCorr_4Hz = h1AmpRaw_4Hz - noiseAmp4;

        if ~isnan(noiseAmp4) && noiseAmp4 > 0
            snr4Hz = h1AmpRaw_4Hz / noiseAmp4;
        else
            snr4Hz = NaN;
        end

        validH1_4Hz = ~isnan(snr4Hz) && snr4Hz > 2;

        if ~validH1_4Hz
            h1AmpCorr_4Hz = NaN;
            h1Phase_4Hz = NaN;
        end

        % -------------------------------------------------------------
        % 5) FEATURE ROW
        % -------------------------------------------------------------
        newRow = table( ...
            string(mouseID), ...
            string(fileName), ...
            numCurves, ...
            baselineMean, ...
            baselineStd, ...
            mean_0_50, ...
            mean_50_100, ...
            mean_100_200, ...
            peakMax, ...
            peakMaxTime, ...
            peakMin, ...
            peakMinTime, ...
            peakToPeak, ...
            rmsVal, ...
            totalEnergy, ...
            NONAmp, NONTime, ...
            PONAmp, PONTime, ...
            LNONAmp, LNONTime, ...
            NOFFAmp, NOFFTime, ...
            POFFAmp, POFFTime, ...
            aWaveAmp, aWaveTime, aWaveSlope, ...
            bWaveAmp, bWaveTime, b_a_ratio, ...
            opRMS, opPeakToPeak, ...
            band_0_10, band_10_30, band_30_80, ...
            domFreq, ...
            h1AmpRaw_4Hz, noiseAmp4, h1AmpCorr_4Hz, h1Phase_4Hz, snr4Hz, double(validH1_4Hz), ...
            'VariableNames', { ...
                'MouseID','FileName','NumCurves', ...
                'BaselineMean','BaselineStd', ...
                'Mean_0_50','Mean_50_100','Mean_100_200', ...
                'PeakMax','PeakMaxTime', ...
                'PeakMin','PeakMinTime', ...
                'PeakToPeak','RMS','TotalEnergy', ...
                'NONAmp','NONTime', ...
                'PONAmp','PONTime', ...
                'LNONAmp','LNONTime', ...
                'NOFFAmp','NOFFTime', ...
                'POFFAmp','POFFTime', ...
                'AWaveAmp','AWaveTime','AWaveSlope', ...
                'BWaveAmp','BWaveTime','BARatio', ...
                'OPRMS','OPPeakToPeak', ...
                'Band_0_10','Band_10_30','Band_30_80', ...
                'DominantFreq', ...
                'H1AmpRaw_4Hz','NoiseAmp_4Hz','H1AmpCorr_4Hz','H1Phase_4Hz','SNR_4Hz','ValidH1_4Hz' ...
            });

        featureTable = [featureTable; newRow];

    catch exception
        fprintf('Hiba: %s\n', exception.message);
    end
end

save(fullfile(parentDir, 'feature_table.mat'), 'featureTable')
disp('Feature table kesz.')

X = featureTable{:, 4:end};
save(fullfile(parentDir, 'feature_matrix.mat'), 'X', 'featureTable')
disp('Feature matrix kesz.')

disp(size(X))

disp('Uj feature-ok ellenorzese:')
disp(featureTable(:, {'MouseID','FileName','NONAmp','PONAmp','LNONAmp','NOFFAmp','POFFAmp','AWaveAmp','BWaveAmp','BARatio'}))

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

Xviz = zscore(Xviz);

[coeff, score, latent, ~, explained] = pca(Xviz);

figure
scatter(score(:,1), score(:,2), 60, 'filled')
xlabel(['PC1 (' num2str(explained(1), '%.1f') '%)'])
ylabel(['PC2 (' num2str(explained(2), '%.1f') '%)'])
title('PCA of ERG feature matrix (visualization only)')
grid on

figure
tiledlayout(2,3)

nexttile
histogram(featureTable.NONAmp)
xlabel('\muV')
title('NON amplitude')

nexttile
histogram(featureTable.PONAmp)
xlabel('\muV')
title('PON amplitude')

nexttile
histogram(featureTable.LNONAmp)
xlabel('\muV')
title('LNON amplitude')

nexttile
histogram(featureTable.NOFFAmp)
xlabel('\muV')
title('NOFF amplitude')

nexttile
histogram(featureTable.POFFAmp)
xlabel('\muV')
title('POFF amplitude')

nexttile
histogram(featureTable.OPRMS)
xlabel('RMS')
title('OP RMS')

% =============================================================
% HELPER
% =============================================================
function y = bandpass_local(x, Fs, f1, f2)
    y = NaN(size(x));

    if numel(x) < 10 || Fs <= 2*f2
        return
    end

    try
        x = fillmissing(x, 'linear', 'EndValues', 'nearest');
        [b, a] = butter(3, [f1 f2] / (Fs/2), 'bandpass');
        y = filtfilt(b, a, x);
    catch
        y = NaN(size(x));
    end
end

function y = lowpass_local(x, Fs, fc)
    y = NaN(size(x));

    if numel(x) < 10 || Fs <= 2*fc
        return
    end

    try
        x = fillmissing(x, 'linear', 'EndValues', 'nearest');
        [b, a] = butter(3, fc/(Fs/2), 'low');
        y = filtfilt(b, a, x);
    catch
        y = NaN(size(x));
    end
end