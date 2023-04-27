function F0 = applyPautoEST2getF0(F, frmsNum, roisNum, windowNum, windowSTT2END, windowPrctile, cntr)
% apply the automatic estimation of percentile to get F0
F0 = NaN(roisNum, frmsNum);
windowSTT2END = cell2mat(cellfun(@(x) x(1):x(2), mat2cell(windowSTT2END.', ones(1, windowNum)), 'UniformOutput', false));
frmsID = 1:frmsNum;
parfor i = 1:roisNum
    tempF = F(i, :);
    interpFunc4F0 = griddedInterpolant(cntr, prctile(tempF(windowSTT2END), windowPrctile(i), 2), 'makima', 'nearest');
    F0(i, :) = interpFunc4F0(frmsID);
end

end

