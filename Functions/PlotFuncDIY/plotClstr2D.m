function [clstrID, clstrCntr] = plotClstr2D(data, clstrNum, method)
%% plot results of clustering (2D)

%% Clustering
if strcmp(method, 'k-means')
    [clstrID, clstrCntr] = kmeans(data, clstrNum);
elseif strcmp(method, 'k-medoids')
    [clstrID, clstrCntr] = kmedoids(data, clstrNum);
elseif strcmp(method, 'GMM')
    optns = statset('Display', 'off', 'MaxIter', 1e5, 'TolFun', 1e-8);
    fitdGMM = fitgmdist(data, clstrNum, 'Options', optns, 'Start', 'plus');
    clstrID = cluster(fitdGMM, data);
    clstrCntr = fitdGMM.mu;
end
[~, clstrCntrID] = sort(vecnorm(clstrCntr - min(data), 2, 2));
clstrCntr = clstrCntr(clstrCntrID, :);
[~, clstrCntrID2] = sort(clstrCntrID);
clstrID = clstrCntrID2(clstrID);

%% Color Setting
color = [...
    180  65  55; ...
    155  90  60; ...
    220 139  55; ...
    239 195  90; ...
    100 175  75; ...
     65 135 125; ...
     50 110 175; ...
    150 100 150]./255;
% red brown orange yellow green cyan blue purple
%  o    x     o      o      o    o    o     x
clr = color([1, 3:7], :);
clrID = round(linspace(1, size(clr, 1), clstrNum));

%% Plotting
hold on;
gscatter(data(:, 1), data(:, 2), clstrID, clr(clrID, :), '.', 16);
hold off;

end

