clc;
clear;
close all hidden;

addpath(genpath('Functions'));
pythonPath = 'C:\Users\Zeta_Eta\.conda\envs\suite2p_v0.12.1\python.exe';

MaxCPUcoreNum = 10;

reExtractON = 1; % re-extract F from data.bin
saveON = 0;
saveON4suite2p = 0;
ROIsIDunstbl = [] + 1; % python Index + 1 | default: first two merged ROIs

if ne(pyenv().Status, 'Loaded')
    pyenv('Version', pythonPath);
end

recInf = readtable('preprocessing\Recording_Information_Summary.xlsx');
recInf = recInf(recInf{:, 1} > 0, [1, 2, 3, 8]);

%% Global Workflow Control
recID = 44;

Date = num2str(recInf.Date(recID));
FOV = recInf.FOV{recID};
TaskInfo = recInf.TaskInfo{recID};
layoutNum = TaskInfo(1:2);

FPS = 15;
FOVpixDim = [512, 512]; % [dimX, dimY] | [Width, Height] | [row, col]
FOVpixNum = prod(FOVpixDim);

loadPath = ['MOTcrctd&ROIdtctdData\', Date, '\', FOV, '\'];
% loadPath = ['D:\Monkey_Bubble\Bubble_2-photon\', loadPath];

% loadPathROI = fullfile(loadPath, 'frms4ROIdtct', 'suite2p\plane0');
loadPathROI = fullfile(loadPath, 'frms4ROIdtct', 'plane0');

savePathROI = ['roiData\', layoutNum, '\'];
if ~exist(savePathROI, 'dir') && saveON
    mkdir(savePathROI);
end

%% Create Pool on Local Machine

if isempty(gcp('nocreate'))
    parpool("local", MaxCPUcoreNum);
elseif gcp('nocreate').NumWorkers ~= MaxCPUcoreNum
    delete(gcp('nocreate'));
    parpool("local", MaxCPUcoreNum);
end


%%
load(fullfile(loadPathROI, 'Fall.mat'), 'stat', 'iscell', 'ops');
stat = stat.';

ROIsIDchsn = iscell(:, 1);
mrgdROIsFull = cellfun(@(x) double(x.imerge.'), stat, 'UniformOutput', false);

ROIsID4mrg = unique(cell2mat(mrgdROIsFull)) + 1;
ROIsIDchsn(ROIsID4mrg) = 0;

ROIsIDmrgd = find(~cellfun(@isempty, mrgdROIsFull));
if isempty(ROIsIDunstbl)
    ROIsIDunstbl = ROIsIDmrgd(1:2);
end
ROIsIDchsn(ROIsIDunstbl) = 0;
ROIsIDmrgd2 = setdiff(ROIsIDmrgd, ROIsIDunstbl);

mrgdROIsPart = mrgdROIsFull(ROIsIDmrgd2);
tempID = tabulate(cell2mat(mrgdROIsPart));
tempID = tempID(tempID(:, 2) > 1, 1);
tempID = cellfun(@(x) any(ismember(x, tempID)), mrgdROIsPart);
tempIDall = 1:sum(tempID);
tempROIs = mrgdROIsPart(tempID);
tempID = ROIsIDmrgd2(tempID);
for i = tempIDall
    tempID(i) = tempID(i).*any(cellfun(@(x) all(ismember(tempROIs{i}, x)), tempROIs(setdiff(tempIDall, i))));
end
ROIsIDovrlpd = tempID(tempID > 0);
ROIsIDchsn(ROIsIDovrlpd) = 0;
ROIsIDmrgd2 = setdiff(ROIsIDmrgd2, ROIsIDovrlpd);
ROIsIDchsn(ROIsIDmrgd2) = 1;

if saveON4suite2p
    iscell(:, 1) = ROIsIDchsn;
    py.numpy.save(fullfile(loadPathROI, 'iscell.npy'), py.numpy.array(iscell));
end

ROIsIDchsn = find(ROIsIDchsn);

ROIsID.unstbl = ROIsIDunstbl;
ROIsID.ovrlpd = ROIsIDovrlpd;
ROIsID.mrgd = ROIsIDmrgd;
ROIsID.chsn = ROIsIDchsn;

%%
load(fullfile(loadPathROI, 'Fall.mat'), 'F');
[~, frmsNum] = size(F);

if reExtractON
    orgID = [ROIsIDchsn; ROIsIDunstbl];
    batchSize = 120;
    binInfo.binPath = fullfile(loadPathROI, 'data.bin');
    binInfo.bitSize = 2;
    binInfo.dataTypeIN = 'int16';
    binInfo.dataTypeOUT = 'single';
    tic;
    [F, ROIsInfo] = reExtractF(stat, orgID, frmsNum, flip(FOVpixDim), batchSize, 0, binInfo);
    toc;
    Fchsn = F(1:length(ROIsIDchsn), :);
    Funstbl = F(length(ROIsIDchsn) + 1:end, :);
else
    Fchsn = F(ROIsIDchsn, :);
    Funstbl = F(ROIsIDunstbl, :);
end

windowInfo.Size = FPS*60;
windowInfo.stepSize = FPS*10;

[windowInfo.Num, windowInfo.stt2end] = autoGroupOverlap(frmsNum, windowInfo.Size, windowInfo.Size - windowInfo.stepSize);
frmsID = 1:frmsNum;
windowInfo.cntr = mean(frmsID(windowInfo.stt2end(1, :) - 1 + (1:windowInfo.Size)'), 1);

[F0, windowInfo] = getF0autoEST2dtrnd(Fchsn, frmsNum, windowInfo, 2.^[9, 11], 'mode');
F0(F0 <= 0) = eps;
dFF0 = Fchsn./F0 - 1;

[F0unstbl, windowInfo] = getF0autoEST2dtrnd(Funstbl, frmsNum, windowInfo, 2.^[9, 11], 'mode');
F0unstbl(F0unstbl <= 0) = eps;
dFF0unstbl = Funstbl./F0unstbl - 1;

if saveON
    load(fullfile(loadPath, 'frmsAvgInfo_Final.mat'), 'outlierID');
    if reExtractON
        save(fullfile(savePathROI, ['roiData-', Date, '-', FOV, '_', TaskInfo, '.mat']), ...
            'dFF0', 'dFF0unstbl', 'Fchsn', 'Funstbl', 'ROIsInfo', 'stat', 'ROIsID', 'ops', 'outlierID');
    else
        save(fullfile(savePathROI, ['roiData-', Date, '-', FOV, '_', TaskInfo, '.mat']), ...
            'dFF0', 'dFF0unstbl', 'Funstbl', 'stat', 'ROIsID', 'ops', 'outlierID');
    end
end

