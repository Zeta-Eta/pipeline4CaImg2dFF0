function [clstrID, clstrCntr, clstrBndr] = plotClstr1D(data, clstrNum, binWidth, faceAlpha, method)
%% plot results of clustering (1D)

%% Clustering
if strcmp(method, 'k-means')
    [clstrID, clstrCntr] = kmeans(data, clstrNum);
    [clstrCntr, ~] = sort(clstrCntr);
    clstrBndr = mean([clstrCntr(1:end-1, 1), clstrCntr(2:end, 1)], 2);
elseif strcmp(method, 'k-medoids')
    [clstrID, clstrCntr] = kmedoids(data, clstrNum);
    [clstrCntr, ~] = sort(clstrCntr);
    clstrBndr = mean([clstrCntr(1:end-1, 1), clstrCntr(2:end, 1)], 2);
elseif strcmp(method, 'GMM')
    modelParams = fitGMMviaCDF(data, clstrNum, 'k-medoids');
    clstrCntr = modelParams.w.';
    [~, tempID] = sort(modelParams.mu);
    w = modelParams.w(tempID);
    mu = modelParams.mu(tempID);
    sigma = modelParams.sigma(tempID);
    clstrBndr = NaN(clstrNum - 1, 1);
    for i = 1:clstrNum - 1
        clstrBndr(i) = findBndr4normPDF(w(i:i+1), mu(i:i+1), sigma(i:i+1));
    end
end
binLimits = [min(data), clstrBndr.', max(data)];
for i = 1:clstrNum
    clstrID(binLimits(i) <= data & data <= binLimits(i+1)) =  i;
end

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
if strcmp(method, 'GMM')
    x = linspace(binLimits(1), binLimits(end), 2^8);
    histogram(data, 'BinLimits', binLimits([1, end]), 'BinWidth', binWidth, 'Normalization', 'pdf', ...
            'EdgeColor', 'none', 'FaceColor', 0.75*[1 1 1], 'FaceAlpha', faceAlpha);
    for i = 1:clstrNum
%         x = linspace(binLimits(i), binLimits(i+1), 2^7);
%         plot(x, w(i).*normpdf(x, mu(i), sigma(i)), 'Color', clr(clrID(i), :));
        patch([x(1), x, x(end)], [0, w(i).*normpdf(x, mu(i), sigma(i)), 0], ...
            clr(clrID(i), :), 'EdgeColor', 'none', 'FaceAlpha', faceAlpha);
    end
    plot(x.', sum(w.*normpdf(x.', mu, sigma), 2), 'r');
else
    for i = 1:clstrNum
        histogram(data(clstrID == i), 'BinLimits', binLimits(i:i+1), 'BinWidth', binWidth, ...
            'EdgeColor', 'none', 'FaceColor', clr(clrID(i), :), 'FaceAlpha', faceAlpha);
    end
end
plot(repmat(clstrBndr.', 2, 1), repmat(ylim.', 1, clstrNum - 1), 'k');
hold off;

end


%% functions

function bndr = findBndr4normPDF(w, mu, sigma)
problem.objective = @(x) w(2).*normpdf(x, mu(2), sigma(2)) - w(1).*normpdf(x, mu(1), sigma(1));
problem.x0 = mean(mu);
problem.solver = 'fzero';
problem.options = optimset('fzero');
bndr = fzero(problem);
end
