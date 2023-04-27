function [N, tempID, tempID2useNum] = autoGroup(totalSize, groupSize)
%% Auto Grouping

N = ceil(totalSize/groupSize);
if N == totalSize/groupSize
    tempID = [1:groupSize:totalSize;
        groupSize:groupSize:totalSize];
else
    tempID = [1:groupSize:totalSize;
        groupSize:groupSize:totalSize, totalSize];
end

tempID2useNum = tempID(2, :) - tempID(1, :) + 1;

end

