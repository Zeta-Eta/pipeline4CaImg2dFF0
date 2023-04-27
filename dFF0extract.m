clc;
clear;
close all hidden;

addpath(genpath('Functions'));
pythonPath = 'C:\Users\Zeta_Eta\.conda\envs\suite2p_v0.12.1\python.exe';

MaxCPUcoreNum = 6;
saveON = 0;
saveON4suite2p = 0;
ROIsIDunstbl = [] + 1; % python Index + 1 | 44: 1734, 1746, 1755

if ne(pyenv().Status, 'Loaded')
    pyenv('Version', pythonPath);
end

recInf = readtable('preprocessing\Recording_Information_Summary.xlsx');
recInf = recInf(recInf{:, 1} > 0, [1, 2, 3, 8]);

%% Global Workflow Control
recID = 45;

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
load(fullfile(loadPathROI, 'Fall.mat'), 'F', 'stat', 'iscell', 'ops');

ROIsIDchsn = iscell(:, 1);
mrgdROIsFull = cellfun(@(x) double(x.imerge), stat, 'UniformOutput', false);

ROIsID4mrg = unique(cell2mat(mrgdROIsFull)) + 1;
ROIsIDchsn(ROIsID4mrg) = 0;

ROIsIDmrgd = find(~cellfun(@isempty, mrgdROIsFull));
if isempty(ROIsIDunstbl)
    ROIsIDunstbl = ROIsIDmrgd(1:2);
end
ROIsIDchsn(ROIsIDunstbl) = 0;
ROIsIDmrgd = setdiff(ROIsIDmrgd, ROIsIDunstbl);

mrgdROIsPart = mrgdROIsFull(ROIsIDmrgd);
tempID = tabulate(cell2mat(mrgdROIsPart));
tempID = tempID(tempID(:, 2) > 1, 1);
tempID = cellfun(@(x) any(ismember(x, tempID)), mrgdROIsPart);
tempIDall = 1:sum(tempID);
tempROIs = mrgdROIsPart(tempID);
tempID = ROIsIDmrgd(tempID);
for i = tempIDall
    tempID(i) = tempID(i).*any(cellfun(@(x) all(ismember(tempROIs{i}, x)), tempROIs(setdiff(tempIDall, i))));
end
ROIsIDovrlpd = tempID(tempID > 0);
ROIsIDchsn(ROIsIDovrlpd) = 0;
ROIsIDmrgd = setdiff(ROIsIDmrgd, ROIsIDovrlpd);
ROIsIDchsn(ROIsIDmrgd) = 1;

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
[~, frmsNum] = size(F);
windowInfo.Size = FPS*60;
windowInfo.stepSize = FPS*10;

[windowInfo.Num, windowInfo.stt2end] = autoGroupOverlap(frmsNum, windowInfo.Size, windowInfo.Size - windowInfo.stepSize);
frmsID = 1:frmsNum;
windowInfo.cntr = mean(frmsID(windowInfo.stt2end(1, :) - 1 + (1:windowInfo.Size)'), 1);

[F0, windowInfo] = getF0autoEST2dtrnd(F(ROIsIDchsn, :), frmsNum, windowInfo, 2.^[9, 11], 'mode');
F0(F0 <= 0) = eps;
dFF0 = F(ROIsIDchsn, :)./F0 - 1;

Funstbl = F(ROIsIDunstbl, :);
[F0unstbl, windowInfo] = getF0autoEST2dtrnd(Funstbl, frmsNum, windowInfo, 2.^[9, 11], 'mode');
F0unstbl(F0unstbl <= 0) = eps;
dFF0unstbl = Funstbl./F0unstbl - 1;

if saveON
    load(fullfile(loadPath, 'frmsAvgInfo_Final.mat'), 'outlierID');
    save(fullfile(savePathROI, ['roiData-', Date, '-', FOV, '_', TaskInfo, '.mat']), 'dFF0', 'dFF0unstbl', 'Funstbl', 'stat', 'ROIsID', 'ops', 'outlierID');
end

