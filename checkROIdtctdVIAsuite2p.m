clc;
clear;
close all hidden;

addpath(genpath('Functions'));
pythonPath = 'C:\Users\Zeta_Eta\.conda\envs\suite2p_v0.12.1\python.exe';

MaxCPUcoreNum = 10;
saveON = 0;

minPixNum = 10;
minPixLen = 100;
maxPixLen = 300;
maxPixPrp = 0.75;

if ne(pyenv().Status, 'Loaded')
    pyenv('Version', pythonPath);
end

recInf = readtable('preprocessing\Recording_Information_Summary.xlsx');
recInf = recInf(recInf{:, 1} > 0, [1, 2, 3, 8]);

%% Global Workflow Control
recID = 46;

Date = num2str(recInf.Date(recID));
FOV = recInf.FOV{recID};
TaskInfo = recInf.TaskInfo{recID};
layoutNum = TaskInfo(1:2);

FPS = 15;
FOVpixDim = [512, 512]; % [dimX, dimY] | [Width, Height] | [row, col]
FOVpixNum = prod(FOVpixDim);

loadPath = ['MOTcrctd&ROIdtctdData\', Date, '\', FOV, '\'];
% loadPath = ['D:\Monkey_Bubble\Bubble_2-photon\', loadPath];

loadPathROI = fullfile(loadPath, 'frms4ROIdtct', 'plane0');


%% Create Pool on Local Machine

if isempty(gcp('nocreate'))
    parpool("local", MaxCPUcoreNum);
elseif gcp('nocreate').NumWorkers ~= MaxCPUcoreNum
    delete(gcp('nocreate'));
    parpool("local", MaxCPUcoreNum);
end


%% Remove outlier ROIs
load(fullfile(loadPathROI, 'Fall.mat'), 'F', 'stat');
% statNum = cellfun(@(x) length(fieldnames(x)), [tempS.stat]);
stat = struct2table(cell2mat(stat));

roisIDerr = cellfun(@range, stat.xpix) >= minPixLen | ...
    cellfun(@range, stat.ypix) >= minPixLen;

roisIDerrID = find(roisIDerr);
a = NaN(length(roisIDerrID), 1);
bwI0 = zeros(flip(FOVpixDim));
for i = 1:length(roisIDerrID)
    bwI = XY2bwI(double(stat.xpix{roisIDerrID(i)} + 1), double(stat.ypix{roisIDerrID(i)} + 1), flip(FOVpixDim), bwI0);
    bwI = regionprops(bwI, 'Image').Image;
    bwI = bwmorph(bwI, 'thicken', 1);
    if mean(bwI, 'all') < maxPixPrp
        roisIDerr(roisIDerrID(i)) = false;
    end
end

roisIDerr = cellfun(@(x) sum(~x), stat.overlap) < minPixNum | ...
    cellfun(@range, stat.xpix) >= maxPixLen | ...
    cellfun(@range, stat.ypix) >= maxPixLen | roisIDerr;
%% Get dFF0 and decayRate
[roisNum, frmsNum] = size(F);
windowInfo.Size = FPS*60;
windowInfo.stepSize = FPS*10;

frmsID = 1:frmsNum;
frmsIDout = [];
% frmsIDorigin = setdiff(frmsID, frmsIDout);
frmsID(frmsIDout) = [];
frmsNum = length(frmsID);
[windowInfo.Num, windowInfo.stt2end] = autoGroupOverlap(frmsNum, windowInfo.Size, windowInfo.Size - windowInfo.stepSize);
windowInfo.cntr = mean(frmsID(windowInfo.stt2end(1, :) - 1 + (1:windowInfo.Size)'), 1);

[~, zeroLen] = cntnsINDandLEN(F <= 0);
zeroLen(cellfun(@isempty, zeroLen)) = {0};
roisIDerr = cellfun(@(x) max(x) >= windowInfo.Size, zeroLen) | roisIDerr;

F0 = zeros(roisNum, frmsNum);
[F0(~roisIDerr, :), windowInfo] = getF0autoEST2dtrnd(F(~roisIDerr, frmsID), frmsNum, windowInfo, 2.^[9, 11], 'mode');

F0(roisIDerr, :) = repmat(mean(F(roisIDerr, :), 2), 1, frmsNum);
[~, zeroLen] = cntnsINDandLEN(F0 <= 0);
zeroLen(cellfun(@isempty, zeroLen)) = {0};
roisIDerr = cellfun(@(x) sum(x) > 0, zeroLen) | roisIDerr;
% roisIDerr = table(find(roisIDerr), zeroLen(roisIDerr), 'VariableNames', {'ID', 'zeroLen'});

F0(F0 <= 0) = eps;

dFF0 = F./F0 - 1;

validID = find(~roisIDerr);
[decayCurve, SNR, outlierInd, theta] = slctROI(dFF0(validID, :), FPS*10, 0.05);
find(isnan(SNR))

decayRate = max(decayCurve.modelParams(:, 3), ...
    (decayCurve.meanData(:, 2) - decayCurve.modelParams(:, 2))./(decayCurve.meanData(:, 1) - decayCurve.modelParams(:, 2)));
figure;
[clstrID, clstrCntr] = plotClstr1D(decayRate, 2, 0.025, 0.5, 'GMM'); % 'GMM' | 'k-means' | 'k-medoids'
iscell = zeros(roisNum, 2);
iscell(validID(clstrID == 2), 1) = 1;
iscell(validID, 2) = decayRate;

dFF0scld = dFF0.*mean(F0, 2);
if saveON
    py.numpy.save(fullfile(loadPathROI, 'iscell.npy'), py.numpy.array(iscell));
    py.numpy.save(fullfile(loadPathROI, 'dFF0scld.npy'), py.numpy.array(dFF0scld));
end

%% Plot F, F0, F-F0, dFF0scld, dFF0 & decayCurve
% close all;
n = 409 + 1;
n2 = find(validID == n);

figure;
zoom xon;
pan xon;
pan off;
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(F(n, :), 'c');
hold on;
plot(F(n, :) - F0(n, :), 'r');
plot(dFF0scld(n, :), 'g');
plot(F0(n, :), 'b');
hold off;

nexttile;
plot(dFF0(n, :), 'g');
hold on;
dFF0peak = zeros(1, frmsNum);
dFF0peak(decayCurve.peakID{n2}) = dFF0(n, decayCurve.peakID{n2});
plot(dFF0peak, 'Color', 0.7*[1 1 1]);
dFF0peak = NaN(1, frmsNum);
dFF0peak(decayCurve.fullID2use{n2}) = dFF0(n, decayCurve.fullID2use{n2});
plot(dFF0peak, 'k');

plot([1, frmsNum], repmat(theta.ngtv(n2), 1, 2), 'r');
plot([1, frmsNum], repmat(theta.pstv(n2), 1, 2), 'b');
hold off;


figure;
y = decayCurve.meanData(n2, :);
y(isnan(y)) = [];
x = 0:length(y) - 1;
plot(x, y, '.:b');
hold on;
xModel = linspace(x(1), x(end), 2^8);
tempParams = decayCurve.modelParams(n2, :);
yModel = tempParams(1).*(tempParams(3).^xModel) + tempParams(2);
plot(xModel, yModel, 'r');
hold off;

%% Functions
function bwI = XY2bwI(X, Y, dimYX, bwI)
bwI(sub2ind(dimYX, Y, X)) = 1;
end
