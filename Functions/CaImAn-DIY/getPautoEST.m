function [windowPrctile, windowPrctileMode] = getPautoEST(F, windowNum, windowSTT2END, prsistntParams1, prsistntParams2)
%% get the automatic estimation of percentile
windowPrctile = NaN(1, windowNum);
for j = 1:windowNum
    [~, pdf, ~, cdf] = kde_DIY(F(windowSTT2END(1, j):windowSTT2END(2, j)), prsistntParams1);
    [~, ind] = max(pdf);
    windowPrctile(j) = cdf(ind)*100;
end
windowPrctile = min(100, max(0, windowPrctile));
[~, pdf, xgrid, ~] = kde_DIY(windowPrctile, prsistntParams2);
[~, ind] = max(pdf);
windowPrctileMode = xgrid(ind);


end

