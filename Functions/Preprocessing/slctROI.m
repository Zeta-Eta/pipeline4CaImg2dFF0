function [decayCurve, SNR, outlierInd, theta, STD] = slctROI(dFF0, decayWindow, decayTheta)
% select ROI based on `findpeaks`
roisNum = size(dFF0, 1);
frmsNum = size(dFF0, 2);

ngtvInd = dFF0 <= 0;
ngtvNum = sum(ngtvInd, 2);
ngtvSTD = NaN(roisNum, 1);
ngtvTheta = NaN(roisNum, 1);
pstvInd = dFF0 >= 0;
pstvNum = sum(pstvInd, 2);
pstvSTD = NaN(roisNum, 1);
pstvTheta = NaN(roisNum, 1);
tic;
parfor i = 1:roisNum
    dFF0temp = dFF0(i, :);
    ngtvSTD(i) = sqrt(sum(dFF0temp(ngtvInd(i, :)).^2, 'omitnan')./(ngtvNum(i) - 1));
    %     ngtvPeaks = -findpeaks(-dFF0temp, 'SortStr', 'descend', 'MinPeakHeight', ngtvSTD(i), 'MinPeakDistance', minPeakDist);
    %     ngtvTheta(i) = (mean(ngtvPeaks) - std(ngtvPeaks) + 0.5*min(dFF0temp))./1.5;
    ngtvTheta(i) = -3*ngtvSTD(i);

    pstvSTD(i) = sqrt(sum(dFF0temp(pstvInd(i, :)).^2, 'omitnan')./(pstvNum(i) - 1));
%     pstvPeaks = findpeaks(dFF0temp, 'SortStr', 'descend', 'MinPeakHeight', pstvSTD(i), 'MinPeakDistance', minPeakDist);
%     pstvTheta(i) = (sqrt(sum(pstvPeaks.^2, 'omitnan')./(length(pstvPeaks) - 1)) + max(dFF0temp))./2;
    pstvTheta(i) = (pstvSTD(i) + max(dFF0temp))./2;
%     pstvTheta(i) = (mean(pstvPeaks) + std(pstvPeaks) + max(dFF0temp))./2;
end
toc;

outlierInd = dFF0 <= ngtvTheta;
noiseSTD = NaN(roisNum, 1);
dM = diff([zeros(roisNum, 1), dFF0 >= pstvTheta , zeros(roisNum, 1)], 1, 2);
% cID = cell(roisNum, 1);
peakID = cell(roisNum, 1);
fullID2use = cell(roisNum, 1);
fullData = cell(roisNum, 1);
meanData = NaN(roisNum, decayWindow + 1);

ops = fitoptions('Method', 'NonlinearLeastSquares');
decayCurveModel2Fit = fittype( @(a, b, c, x) a.*(c.^x) + b, ...
    'independent', 'x', ...
    'dependent', 'y', ...
    'options', ops);
decayCurveFit = @(x, y) fit(x, y, decayCurveModel2Fit, ...
    'StartPoint', [y(1) - y(end), y(end), 0.75], ...
    'Lower', [0, -y(1), 0], ...
    'Upper', [2*y(1), y(1), 1], ...
    'MaxFunEvals', 1000, 'MaxIter', 500, 'TolFun', 1e-7, 'TolX', 1e-7);

tempXdata = 0:decayWindow;
modelParams = NaN(roisNum, 3);
Rsqrd = NaN(roisNum, 1);
tic;
parfor i = 1:roisNum
    if ~isnan(pstvTheta(i))
        dMtemp = dM(i, :);
        dFF0temp = dFF0(i, :);

        cIDtemp = [find(dMtemp == 1); find(dMtemp == -1) - 1];
        cIDnum = size(cIDtemp, 2);
        peakIDtemp = [NaN(1, cIDnum), frmsNum];
        for j = cIDnum:-1:1
            [~, peakIDtemp(j)] = max(dFF0temp(cIDtemp(1, j):cIDtemp(2, j)));
            peakIDtemp(j) = cIDtemp(1, j) - 1 + peakIDtemp(j);
            if peakIDtemp(j + 1) - peakIDtemp(j) <= decayWindow
                peakIDtemp(j) = [];
                cIDtemp(:, j) = [];
            end
        end
        peakID{i} = peakIDtemp(1:end-1);
        fullData{i} = dFF0temp(peakID{i}.' + tempXdata);
        if length(peakID{i}) > 10
            [clstrID, clustrCntr] = kmeans(fullData{i}(:, 1), 2);
            [~, clstrOdr] = max(clustrCntr);
            tempYdata = mean(double(fullData{i}(clstrID == clstrOdr, :)), 1, 'omitnan');
        else
            tempYdata = mean(double(fullData{i}), 1, 'omitnan');
        end
%         meanData(i, :) = tempYdata;
        modelParamsTemp = coeffvalues(decayCurveFit(tempXdata.', tempYdata.'));
        tempXtheta = max(min(ceil(log(decayTheta)./log(modelParamsTemp(3))), decayWindow - round(0.1*decayWindow)), 15);
        tempData = fullData{i}(:, tempXtheta + 2:end);
        tempData = tempData(:);
        tempDataPstv = tempData >= 0;
        pstvTheta(i) = 3*sqrt(sum(tempData(tempDataPstv).^2, 'omitnan')./(sum(tempDataPstv) - 1));
        
        tempXdata2 = 0:tempXtheta;
        dMtemp = diff([0, dFF0temp >= pstvTheta(i) , 0], 1, 2);
        cIDtemp = [find(dMtemp == 1); find(dMtemp == -1) - 1];
        cIDnum = size(cIDtemp, 2);
        peakIDtemp = [NaN(1, cIDnum), frmsNum];
        for j = cIDnum:-1:1
            [~, peakIDtemp(j)] = max(dFF0temp(cIDtemp(1, j):cIDtemp(2, j)));
            peakIDtemp(j) = cIDtemp(1, j) - 1 + peakIDtemp(j);
            if peakIDtemp(j + 1) - peakIDtemp(j) <= tempXtheta
                peakIDtemp(j) = [];
                cIDtemp(:, j) = [];
            end
        end
%         cID{i} = cIDtemp;
        peakID{i} = peakIDtemp(1:end-1);
        fullData{i} = dFF0temp(peakID{i}.' + tempXdata2);
        if length(peakID{i}) > 10
            [clstrID, clustrCntr] = kmeans(fullData{i}(:, 1), 2);
            [~, clstrOdr] = max(clustrCntr);
            peakID2use = clstrID == clstrOdr;
            fullID2use{i} = reshape(peakID{i}(peakID2use) + tempXdata2.', 1, []);
            tempYdata = mean(double(fullData{i}(peakID2use, :)), 1, 'omitnan');
        else
            fullID2use{i} = reshape(peakID{i} + tempXdata2.', 1, []);
            tempYdata = mean(double(fullData{i}), 1, 'omitnan');
        end
        tempYdata2 = NaN(1, decayWindow + 1);
        tempYdata2(tempXdata2 + 1) = tempYdata;
        meanData(i, :) = tempYdata2;
        [decayCurveModelFitd, GoF] = decayCurveFit(tempXdata2.', tempYdata.');
        modelParamsTemp = coeffvalues(decayCurveModelFitd);
        modelParams(i, :) = modelParamsTemp;
        Rsqrd(i) = GoF.rsquare;
        cIDtemp = setdiff(setdiff(1:frmsNum, find(ngtvInd(i, :) | dFF0temp >= pstvTheta(i))), fullID2use{i});
        noiseSTD(i) = sqrt(sum(dFF0temp(cIDtemp).^2, 'omitnan')./(length(cIDtemp) - 1));
    end
end
toc;
SNR = meanData(:, 1)./noiseSTD./3;
decayCurve = table(fullID2use, peakID, fullData, meanData, modelParams, Rsqrd);

STD.ngtv = ngtvSTD;
STD.pstv = pstvSTD;
theta.ngtv = ngtvTheta;
theta.pstv = pstvTheta;


end

