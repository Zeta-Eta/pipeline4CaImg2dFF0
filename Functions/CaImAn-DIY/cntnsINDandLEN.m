function [cInd, cLen] = cntnsINDandLEN(M)
%% continuous length and index
rowNum = size(M, 1);
dM = diff([zeros(rowNum, 1), M, zeros(rowNum, 1)], 1, 2);
cInd = cell(rowNum, 1);
for rowID = 1:rowNum
    dMtemp = dM(rowID, :);
    cInd{rowID} = [find(dMtemp == 1); find(dMtemp == -1)];
end
cLen = cellfun(@(x) diff(x, 1, 1), cInd, 'UniformOutput', false);

end

