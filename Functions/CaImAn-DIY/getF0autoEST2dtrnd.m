function [F0, windowInfo] = getF0autoEST2dtrnd(F, frmsNum, windowInfo, n, prctileMethod)
%% get the automatic estimation of F0 to detrend
roisNum = size(F, 1);
windowNum = windowInfo.Num;
windowSTT2END = windowInfo.stt2end;

windowPrctile = NaN(roisNum, windowNum);
windowPrctileMode = NaN(roisNum, 1);
prsistntParams1 = setPrsistntParams(n(1), 7);
prsistntParams2 = setPrsistntParams(n(2), 7);
tic
parfor i = 1:roisNum
    [windowPrctile(i, :), windowPrctileMode(i)] = getPautoEST(F(i, :), windowNum, windowSTT2END, prsistntParams1, prsistntParams2);
end
toc
windowInfo.prctile = windowPrctile;
windowInfo.prctileMode = min(100, max(0, windowPrctileMode));
windowInfo.prctileMean = mean(windowPrctile, 2);
windowInfo.prctileMedian = median(windowPrctile, 2);
clear windowPrctileMode;
if strcmp(prctileMethod, 'mode')
    windowPrctile = windowInfo.prctileMode;
elseif strcmp(prctileMethod, 'mean')
    windowPrctile = windowInfo.prctileMean;
else
    windowPrctile = windowInfo.prctileMedian;
end

tic
F0 = applyPautoEST2getF0(F, frmsNum, roisNum, windowNum, windowSTT2END, windowPrctile, windowInfo.cntr);
toc






end

