function [F, ROIsInfo] = reExtractF(stat, orgID, frmsNum, imgSize, batchSize, overlapON, binInfo)
% re-extract F from data.bin
stat = reshape(stat(orgID), [], 1);
roisNum = length(stat);

%% indices to binary images
bwI0 = false(imgSize);
bwI = cellfun(@(roi) XY2I(roi.xpix + 1, roi.ypix + 1, imgSize, bwI0, 1), stat, 'UniformOutput', false);
bwImat = cell2mat(reshape(bwI, 1, 1, []));

lambdaI0 = zeros(imgSize, binInfo.dataTypeOUT);
lambdaI = cellfun(@(roi) XY2I(roi.xpix + 1, roi.ypix + 1, imgSize, lambdaI0, roi.lam), stat, 'UniformOutput', false);

%% get overlap info
overlapI = sum(bwImat, 3) > 1;
overlapROIind = cell(imgSize);
tempADDind = (0:roisNum - 1).*prod(imgSize);
for i = find(overlapI).'
    overlapROIind{i} = find(bwImat(i + tempADDind));
end
overlapInfo = cellfun(@(roi) bwI2overlapInfo(roi, overlapI, overlapROIind), bwI, 'UniformOutput', false);

%% get pixID
if overlapON
    pixInd4extract = bwImat;
else
    pixInd4extract = ~overlapI & bwImat;
end

pixIndCell = cell(roisNum, 1);
wCell = cell(roisNum, 1);
for roisInd = 1:roisNum
    pixInd = pixInd4extract(:, :, roisInd).';
    pixIndCell{roisInd} = find(pixInd);
    lambdaItemp = lambdaI{roisInd}.';
    w = max(lambdaItemp(pixInd), 0);
    wCell{roisInd} = w.'./sum(w);
end

%% re-extract F from data.bin
[N, tempID, tempID2useNum] = autoGroup(frmsNum, batchSize);
FOVpixDim = flip(imgSize);
F = NaN(roisNum, frmsNum, binInfo.dataTypeOUT);
for n = 1:N
    frmsTemp = read_bin_file_DIY(binInfo.binPath, tempID(1, n), tempID2useNum(n), ...
        FOVpixDim, binInfo.bitSize, binInfo.dataTypeIN, binInfo.dataTypeOUT, 2);
    frmsInd = tempID(1, n):tempID(2, n);
    for roisInd = 1:roisNum
        F(roisInd, frmsInd) = wCell{roisInd} * frmsTemp(pixIndCell{roisInd}, :);
    end
end

%% summary output
ROIsInfo = table(orgID, bwI, lambdaI, overlapInfo);


end


%% sub-functions

function I = XY2I(X, Y, dimYX, I, input)
I(sub2ind(dimYX, Y, X)) = input;
end

function overlapInfo = bwI2overlapInfo(bwI, overlapI, overlapROIind)
tempID = find(overlapI & bwI);
overlapInfo = [num2cell(tempID), overlapROIind(tempID)];
end

