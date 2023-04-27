function [groupNum, tempID, overlapSize] = autoGroupOverlap(totalSize, groupSize, overlapSize)
%% Auto Grouping
groupNum = 1 + round((totalSize - groupSize)/(groupSize - overlapSize));

tempID = round(linspace(0, totalSize - groupSize, groupNum));
tempID = [tempID + 1; tempID + groupSize];

overlapSize = tempID(2, 1:end - 1) - tempID(1, 2:end) + 1;

end

