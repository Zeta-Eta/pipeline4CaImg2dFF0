clc;
clear;
close all hidden;

addpath(genpath('Functions'));
pythonPath = 'C:\Users\Zeta_Eta\.conda\envs\suite2p_v0.12.1\python.exe';
suite2pVersion = 'v0.12.1';

MaxCPUcoreNum = 10;
MaxCPUcoreNum4SPMD1 = 5; % Step 2|4
MaxCPUcoreNum4SPMD2 = 3; % Step 6
MaxCPUcoreNum4SPMD3 = 8; % Step 8
deletePoolON = 0;

recInf = readtable('preprocessing\Recording_Information_Summary.xlsx');
recInf = recInf(recInf{:, 1} > 0, [1, 2, 3, 8]);

peakBRT = 4000; % peak brightness you set to see during recording

%% Global Workflow Control
% Step 0: TIFF to HDF5
% Step 1: Uniform Downsampling for Generating an Average Frame (the 1st Template)
% Step 2: The 1st Motion Correction for Generating an Average Frame (the 2nd Template)
% Step 3: Check the Results of the 1st Motion Correction
% Step 4: The 2nd Motion Correction to Top 50 Frames in each Trial (rigid)
% Step 5: Check the Results of the 2nd Motion Correction
% Step 6: The 3rd Motion Correction and Interpolation for Generating Non-rigid Displacement Fields
% Step 7: Check the Results of the 3rd Motion Correction
% Step 8: The 4th Motion Correction to All Frames (non-rigid)
% Step 9: Check the Results of the 4th Motion Correction
% Step 10: Prepare `ops.npy` for suite2p
recID = 55;
stepNum = 1:10; % from 0 to 10
saveImageON = 1;
closeFigureON = 1;
playMovieON = [0, 1]; % [Step 3, Step 5|9]
moveOriH5ON = [0, 1];
trlNum4avg = 200; % Step 1
thresh4avg = 0.95; % Step 1
frmsResampNumByTrl = 50; % Step 1|4|5|6
frmsResampIDsprs = [3, 9, 18, 31, 44]; % Step 1
fitgmdistTolFun = [1e-8, 1e-8]; % Step 4|8
fitgmdistRegVal = [1e-3, 1e-3]; % Step 4|8
trlsBinSize = 50; % Step 6
trlsBinDist =  2; % Step 6
CLAHEutmostON = 1; % Step 6
smoothMethod = 'rloess'; % Step 7 | 'rloess' or 'rlowess'
smoothWindow = repmat(round(trlsBinSize/trlsBinDist), 1, 2); % Step 7 | `repmat(round(trlsBinSize/trlsBinDist), 1, 2)`
extrapMethod = 'linearDIY'; % Step 7 | 'linearDIY' or 'nearest'

iterNum4step1 =  5; % iteration number for Step 1, set it to be at least 2. 
iterNum4step2 =  2; % iteration number for Step 2, set it to be at least 2. 
iterNum4step6 =  2; % iteration number for Step 6, set it to be at least 2 (but not too large). 

Date = num2str(recInf.Date(recID));
FOV = recInf.FOV{recID};
TaskInfo = recInf.TaskInfo{recID};
layoutNum = TaskInfo(1:2);

deleteH5ON = 1; % At the end of Step 3, delete H5 generated in Step 2. 
fpsCrctON = 0; % Step 0: transform X fps to X/2 fps

% Step 5: 
MSCmethod = 'BIC'; % Model Select Criteria -> 'RSS'(recommend) or 'BIC'
CompNum2use0 = []; % Number of Components to Use -> an int from `1:8` or just `[]`(calculated automatically)

tiff2hdf5ON       = ismember(0, stepNum);
UNIFdnsampON      = ismember(1, stepNum);
motCrct1stON      = ismember(2, stepNum);
motCrct1stCheckON = ismember(3, stepNum);
motCrct2ndON      = ismember(4, stepNum);
motCrct2ndCheckON = ismember(5, stepNum);
motCrct3rdON      = ismember(6, stepNum);
motCrct3rdCheckON = ismember(7, stepNum);
motCrct4thON      = ismember(8, stepNum);
motCrct4thCheckON = ismember(9, stepNum);
opsPrpr4suite2pON = ismember(10, stepNum);


FPS = 15;
FOVpixDim = [512, 512]; % [dimX, dimY] | [Width, Height] | [colNum, rowNum]
FOVpixNum = prod(FOVpixDim);
inputDataType = 'uint16';
outputFormat = 'h5';
% ↑ Time x Height x Width in HDF5   (C order)
%  Width x Height x Time  in MATLAB (Fortran order)
h5GroupName = 'data'; % set to 'data' for the convenience of suite2p users.
outputFormatFinal = 'bin'; % could skip the process of converting data from h5 file to bin file
outputDataTypeFinal = 'int16'; % make sure the max of dynamic range will never be out of the datatype range. (in our data, the max is 2^13)
bitSize = 2;

tempH5Name = 'tempH5\temp.h5';
loadPathRAW = ['rawData\', Date, '\'];
loadPathBHV = ['bhvData\', layoutNum, '\'];
savePathH5 = 'rawDataH5\';
savePath = ['MOTcrctd&ROIdtctdData\', Date, '\', FOV, '\'];

% loadPath = ['D:\Monkey_Bubble\Bubble_2-photon\', loadPath];
savePathH5 = ['R:\Monkey_Bubble\Bubble_2-photon\', savePathH5];
% loadPathH5 = savePathH5;
loadPathH5 = 'tempH5\';

% savePath = ['D:\Monkey_Bubble\Bubble_2-photon\', savePath];

frmsFolderDir = dir([loadPathRAW, 'TSeries*-', FOV, '-*']);
if ~exist(savePath, 'dir') && motCrct1stON
    mkdir(savePath);
end

MaxRAM2use = 40; % unit: GB
% ↑ recommended range: 30%-70% of RAM
%   and better to be the integral multiples of 5
%   e.g.: 5/8GB | 5~10/16GB | 10~20/32GB | 20~45/64GB | 40~90/128GB
scaleInt = 1000; % or 1024
MaxFrmsNum = @(byteNum, pixNum) scaleInt*round(MaxRAM2use*1024.^2/(byteNum*pixNum));
MaxFrmsNum4double = MaxFrmsNum(8, FOVpixNum); % 'double' memSize: 8 byte
MaxFrmsNum4single = MaxFrmsNum(4, FOVpixNum); % 'single' memSize: 4 byte
MaxFrmsNum4uint16 = MaxFrmsNum(2, FOVpixNum); % 'uint16' memSize: 2 byte
MaxFrmsNum4dAND16 = MaxFrmsNum(8 + 2, FOVpixNum);

batchSize4load1st = [1,  20]; % batch size for load/run in Step 2
batchSize4load2nd = [1,  20]; % batch size for load/run in Step 4
batchSize4load3rd = [1,   1]; % batch size for load/run in Step 6 (the 2nd is actually bin_num)
batchSize4load4th = [1, 100]; % batch size for load/run in Step 8
% ↑ higher batchSize4load will make the Motion Correction load/run slightly faster, 
%                                   [↑ before attaining CPU performance peak]
%   but also needs slightly larger RAM. 
batchSize4load5th = MaxFrmsNum4dAND16; % chunk size for load/run in Step 9

chunkSize4save1st = 100; % chunk size for save in Step 0, 1, 2
chunkSize4save2nd = 100; % chunk size for save in Step 8
% ↑ higher chunkSize4save will make the Motion Correction save slightly faster, 
%   but also needs larger RAM and usually output larger files.

method2getTemplate = {'mean'; 'mean'}; % 'median' or 'mean'

suffix0 = '_Rigid';
options4motCrct = NoRMCorreSetParms('d1', FOVpixDim(1), 'd2', FOVpixDim(2), ...
    'grid_size', FOVpixDim, 'overlap_pre', 0, ...
    'init_batch', batchSize4load1st(1), 'bin_width', batchSize4load1st(2), ...
    'mem_batch_size', chunkSize4save1st, 'method', method2getTemplate, ...
    'us_fac', 100, 'max_dev', [Inf, Inf], 'max_shift', [Inf, Inf], ...
    'output_type', outputFormat, 'h5_groupname', h5GroupName, ...
    'col_shift', 0, 'upd_template', false, 'print_msg', 10);

suffix1 = '_nonRigid';


screenPixSize = get(0, 'ScreenSize');
basePixSize = 1080; % never change it
scaleFactor = screenPixSize(4)/basePixSize;

CLAHEslight = @(I) adapthisteq(rescale(I), 'NumTiles', repmat(2^5, 1, 2), ...
    'Distribution', 'rayleigh', 'Alpha', 0.39, 'NBins',  2^9, 'Range', 'original');
CLAHEstrong = @(I) adapthisteq(rescale(I), 'NumTiles', repmat(2^5, 1, 2), ...
    'Distribution', 'rayleigh', 'Alpha', 0.32, 'NBins', 2^16, 'Range', 'original');
ORGNLgreen = @(I) cat(3, zeros(size(I)), rescale(min(I, peakBRT)), zeros(size(I)));

CLAHEutmost = @(I) adapthisteq(rescale(I), 'NumTiles', repmat(2^5, 1, 2), ...
    'Distribution', 'rayleigh', 'Alpha', 0.11777, 'NBins', 2^16, 'Range', 'original');
% ↑ maxRange = Alpha*sqrt(-2*log(eps))
% If Alpha < 1/sqrt(-2*log(eps)) ≈ 0.11777, 
% `CLAHEutmost` won't change any more after being rescaled again. 


%% Create Pool on Local Machine
if any(stepNum ~= 10)
    if isempty(gcp('nocreate'))
        parpool("local", MaxCPUcoreNum);
    elseif gcp('nocreate').NumWorkers ~= MaxCPUcoreNum
        delete(gcp('nocreate'));
        parpool("local", MaxCPUcoreNum);
    end
end

tStart = tic;

%% Color Defination
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
%  o    o     o      o      o    o    o     o
clr = color(:, :);
clrN = size(clr, 1);



%% Set Load Info
if tiff2hdf5ON
    frmsFolderNum = length(frmsFolderDir);
    frmsDiffNum = NaN(frmsFolderNum, 1);
    frmsFullDir = cell(frmsFolderNum, 1);
    parfor n = 1:frmsFolderNum
        frmsFullDir{n} = dir(fullfile(frmsFolderDir(n).folder, frmsFolderDir(n).name, '*.tif'));
        frmsDiffNum(n) = length(frmsFullDir{n});
    end
    frmsFullDir = cat(1, frmsFullDir{:});
    frmsFullNum = sum(frmsDiffNum);

    saveName4rawDataInfo = fullfile(savePathH5, ['rawDataInfo-', Date, '-', FOV, '.mat']);

    if fpsCrctON
        frmsFullDirORGNL = frmsFullDir;
        frmsDiffNumORGNL = frmsDiffNum;
        frmsFullNumORGNL = frmsFullNum;
        frmsDiffNum = floor(frmsDiffNumORGNL./2);
        frmsFullNum = sum(frmsDiffNum);
        
        frmsFullDir = frmsFullDirORGNL;
        tempID = cumsum(frmsDiffNumORGNL);
        frmsFullDir(tempID(mod(frmsDiffNumORGNL, 2) == 1)) = [];
    end

    cutON = 0;
    if mod(frmsFullNum, chunkSize4save1st) == 1
        frmsFullNum = frmsFullNum - 1;
        frmsDiffNum(end) = frmsDiffNum(end) - 1;
        frmsFullDir(end) = [];
        cutON = 1;
    end

    save(saveName4rawDataInfo, 'frmsFullDir', 'frmsDiffNum', 'frmsFullNum');

    if fpsCrctON
        save(saveName4rawDataInfo, 'frmsFullDirORGNL', 'frmsDiffNumORGNL', ...
            'frmsFullNumORGNL', '-append');
    end

    if cutON
        save(saveName4rawDataInfo, 'cutON', '-append');
    end
else
    loadName4rawDataInfo = fullfile(savePathH5, ['rawDataInfo-', Date, '-', FOV, '.mat']);
    load(loadName4rawDataInfo, 'frmsFullDir', 'frmsDiffNum', 'frmsFullNum');
end

if moveOriH5ON(1)
    copyfile(fullfile(savePathH5, ['rawData-', Date, '-', FOV, '.', outputFormat]), loadPathH5);
end

%% Step 0: TIFF to HDF5
if tiff2hdf5ON
    saveName4rawData = fullfile(savePathH5, ['rawData-', Date, '-', FOV, '.', outputFormat]);
    h5create(saveName4rawData, ['/', h5GroupName], [FOVpixDim, Inf], ...
        'ChunkSize', [FOVpixDim, chunkSize4save1st], 'Datatype', inputDataType);

    [N, tempID, tempID2useNum] = autoGroup(frmsFullNum, chunkSize4save1st);
    if fpsCrctON
        frmsFullDirODD = frmsFullDir(1:2:end);
        frmsFullDirEVEN = frmsFullDir(2:2:end);
        for n = 1:N
            tempIDusedNum = tempID(1, n) - 1;

            frmsTemp = zeros([FOVpixDim, tempID2useNum(n)], inputDataType);
            parfor m = tempID(1, n):tempID(2, n)
                frmsTemp(:, :, m - tempIDusedNum) = ...
                    (imread(fullfile(frmsFullDirODD(m).folder, frmsFullDirODD(m).name)) ...
                    + imread(fullfile(frmsFullDirEVEN(m).folder, frmsFullDirEVEN(m).name))).'./2;
            end
            h5write(saveName4rawData, ['/', h5GroupName], frmsTemp, ...
                [1, 1, tempID(1, n)], [FOVpixDim, tempID2useNum(n)]);
        end
    else
        for n = 1:N
            tempIDusedNum = tempID(1, n) - 1;

            frmsTemp = zeros([FOVpixDim, tempID2useNum(n)], inputDataType);
            parfor m = tempID(1, n):tempID(2, n)
                frmsTemp(:, :, m - tempIDusedNum) = imread(fullfile(frmsFullDir(m).folder, frmsFullDir(m).name)).';
            end
            h5write(saveName4rawData, ['/', h5GroupName], frmsTemp, ...
                [1, 1, tempID(1, n)], [FOVpixDim, tempID2useNum(n)]);
        end
    end
    clear frmsTemp;

end

%% Step 1: Uniform Downsampling for Generating an Average Frame (the 1st Template)
if UNIFdnsampON

    load(fullfile(loadPathBHV, ['bhvData-', Date, '-', FOV, '_', TaskInfo, '.mat']), 'dataset');
    crctTrlID = find(~dataset.errorType & ~any(isnan(dataset.eventFrms{:, [1 end]}), 2));
    trlID4avg = round(linspace(1, length(crctTrlID), trlNum4avg + 2));


    trlID4avg([1, end]) = [];
    trlID4avg = crctTrlID(trlID4avg);
    tempID = dataset.eventFrms.fixON(trlID4avg).' - 1;
    frmsID4avgFull = reshape(tempID + (1:frmsResampNumByTrl).', 1, []);
    frmsID4avgSprs = reshape(tempID + frmsResampIDsprs.', 1, []);

    frmsID4avg = frmsID4avgSprs;
    frmsNum4avg = length(frmsID4avg);
    tempNum = round(thresh4avg*frmsNum4avg);

    loadName4rawData = fullfile(loadPathH5, ['rawData-', Date, '-', FOV, '.', outputFormat]);

    frms4avgORGNL = zeros([flip(FOVpixDim), frmsNum4avg], inputDataType);
    parfor n = 1:frmsNum4avg
        frms4avgORGNL(:, :, n) = h5read(loadName4rawData, ...
            ['/', h5GroupName], [1, 1, frmsID4avg(n)], [FOVpixDim, 1]);
    end

    frmsAvgdORGNL = mean(frms4avgORGNL, 3);
    [~, tempID] = sort(corr(double(reshape(frms4avgORGNL, FOVpixNum, [])), frmsAvgdORGNL(:)));
    frmsAvgdORGNL = mean(frms4avgORGNL(:, :, tempID(tempNum + 1:end)), 3);

    h5create(tempH5Name, ['/', h5GroupName], [FOVpixDim, Inf], ...
        'Chunksize', [FOVpixDim, frmsNum4avg], 'Datatype', inputDataType);
    h5write(tempH5Name, ['/', h5GroupName], ...
        frms4avgORGNL, [1, 1, 1], [FOVpixDim, frmsNum4avg]);
    [~, ~, frmsAvgdORGNL, ~] = ...
        normcorre_batch_Rigid_4GPU(tempH5Name, ...
        NoRMCorreSetParms(options4motCrct), ...
        frmsAvgdORGNL, MaxCPUcoreNum4SPMD1, 1, 0, iterNum4step1, 0, 0);
    frmsAvgdORGNL = frmsAvgdORGNL.';
    delete(tempH5Name);


    frmsID4avg = frmsID4avgFull;
    frmsNum4avg = length(frmsID4avg);
    
    frms4avgORGNL = zeros([flip(FOVpixDim), frmsNum4avg], inputDataType);
    tic;
    parfor n = 1:frmsNum4avg
        frms4avgORGNL(:, :, n) = h5read(loadName4rawData, ...
            ['/', h5GroupName], [1, 1, frmsID4avg(n)], [FOVpixDim, 1]);
    end
    toc;
    h5create(tempH5Name, ['/', h5GroupName], [FOVpixDim, Inf], ...
        'Chunksize', [FOVpixDim, chunkSize4save1st], 'Datatype', inputDataType);
    h5write(tempH5Name, ['/', h5GroupName], ...
        frms4avgORGNL, [1, 1, 1], [FOVpixDim, frmsNum4avg]);

    corrORGNL = corr(double(reshape(frms4avgORGNL, FOVpixNum, [])), frmsAvgdORGNL(:));
    clear frms4avgORGNL;

    if ~motCrct1stCheckON
        figure;
        set(gcf, 'Position', ...
            [10, 50, ...
            560, 315].*scaleFactor);
        xLim = [0 frmsNum4avg + 1];

        t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;
        plot(corrORGNL, 'Color', '#0072BD');
        title('Chronological Order');
        xlim(xLim);

        nexttile;
        plot(sort(corrORGNL), 'Color', '#0072BD');
        title('ZNCC Ascending Order');
        xlim(xLim);

        xlabel(t, 'frame number')
        ylabel(t, 'ZNCC ( a.k.a. Pearson''s {\itr} )')

        figure;
        set(gca, 'Position', [0, 0, 1, 1]);
        set(gcf, 'Position', ...
            [10, basePixSize - FOVpixDim(2) - 110, ...
            FOVpixDim].*scaleFactor);
        imshow(CLAHEslight(frmsAvgdORGNL));
    end


    saveNameAvgInfo_Rigid = fullfile(savePath, 'frmsAvgInfo_Rigid.mat');
    save(saveNameAvgInfo_Rigid, 'thresh4avg', 'frmsNum4avg', 'frmsID4avg', ...
        'corrORGNL', 'frmsAvgdORGNL');

    if closeFigureON && max(stepNum) ~= 1
        close all;
    end
end

%% Step 2: The 1st Motion Correction for Generating an Average Frame (the 2nd Template)
if motCrct1stON
    loadNameAvgInfo_Rigid = fullfile(savePath, 'frmsAvgInfo_Rigid.mat');
    load(loadNameAvgInfo_Rigid, 'frmsNum4avg');

    saveName4avgULTMT = fullfile(savePath, ['frms4avgULTMT', suffix0, '.', outputFormat]);
    [DIYparams, ~, ~, ~] = ...
        normcorre_batch_Rigid_4GPU(tempH5Name, ...
        NoRMCorreSetParms(options4motCrct, 'output_filename', saveName4avgULTMT), ...
        frmsAvgdORGNL.', MaxCPUcoreNum4SPMD1, 1, 0, iterNum4step2, ...
        max(round(frmsNum4avg/(MaxFrmsNum4dAND16*MaxCPUcoreNum4SPMD1)), 1), 0);
    delete(tempH5Name);

    corrULTMT = DIYparams.ZNCC;
    frmsAvgdULTMT = DIYparams.templateOutput.';

    save(loadNameAvgInfo_Rigid, 'corrULTMT', 'frmsAvgdULTMT', '-append');

    if any(stepNum >= 4)
        delete(gcp('nocreate'));
        parpool("local", MaxCPUcoreNum);
    end

    if closeFigureON && max(stepNum) ~= 2
        close all;
    end
end

%% Step 3: Check the Results of the 1st Motion Correction
if motCrct1stCheckON
    loadNameAvgInfo_Rigid = fullfile(savePath, 'frmsAvgInfo_Rigid.mat');
    load(loadNameAvgInfo_Rigid, 'frmsNum4avg', 'corrORGNL', 'frmsAvgdORGNL', 'corrULTMT', 'frmsAvgdULTMT');

    figTEMP = figure;
    set(gcf, 'Position', ...
        [10, 50, ...
        560, 315].*scaleFactor);
    xLim = [0 frmsNum4avg + 1];
    t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(corrORGNL, 'Color', '#0072BD');
    hold on;
    plot(corrULTMT, 'Color', '#77AC30');
    hold off;
    title('Chronological Order');
    xlim(xLim);

    nexttile;
    lgd(1) = plot(sort(corrORGNL), 'Color', '#0072BD');
    hold on;
    lgd(2) = plot(sort(corrULTMT), 'Color', '#77AC30');
    hold off;
    legend(lgd, {'b4motCrct'; 'motCrcted'}, 'box', 'off', ...
        'Location', 'east');
    title('ZNCC Ascending Order');
    xlim(xLim);

    xlabel(t, 'frame number')
    ylabel(t, 'ZNCC ( a.k.a. Pearson''s {\itr} )')
    
    figTEMP.Renderer = 'painters';

    if saveImageON
        imwrite(CLAHEslight(frmsAvgdORGNL), ...
            fullfile(savePath, 'frmsAvgdORGNL_CLAHEslight_uint8.tif'), ...
            'Compression', 'none');
        imwrite(CLAHEstrong(frmsAvgdORGNL), ...
            fullfile(savePath, 'frmsAvgdORGNL_CLAHEstrong_uint8.tif'), ...
            'Compression', 'none');
        imwrite(ORGNLgreen(frmsAvgdORGNL), ...
            fullfile(savePath, 'frmsAvgdORGNL_ORGNLgreen_uint8.tif'), ...
            'Compression', 'none');

        saveas(figTEMP, fullfile(savePath, 'ZNCC4avg_Rigid.svg'));

        imwrite(CLAHEslight(frmsAvgdULTMT), ...
            fullfile(savePath, 'frmsAvgdULTMT_CLAHEslight_uint8.tif'), ...
            'Compression', 'none');
        imwrite(CLAHEstrong(frmsAvgdULTMT), ...
            fullfile(savePath, 'frmsAvgdULTMT_CLAHEstrong_uint8.tif'), ...
            'Compression', 'none');
        imwrite(ORGNLgreen(frmsAvgdULTMT), ...
            fullfile(savePath, 'frmsAvgdULTMT_ORGNLgreen_uint8.tif'), ...
            'Compression', 'none');
    end

    figure;
    set(gca, 'Position', [0, 0, 1, 1]);
    set(gcf, 'Position', ...
        [10, basePixSize - FOVpixDim(2) - 110, ...
        FOVpixDim].*scaleFactor);
    imshow(CLAHEslight(frmsAvgdORGNL));

    figure;
    set(gca, 'Position', [0, 0, 1, 1]);
    set(gcf, 'Position', ...
        [20 + FOVpixDim(1), basePixSize - FOVpixDim(2) - 110, ...
        FOVpixDim].*scaleFactor);
    imshow(CLAHEslight(frmsAvgdULTMT));

    if playMovieON(1)
        loadName4avgULTMT = fullfile(savePath, ['frms4avgULTMT', suffix0, '.', outputFormat]);
        frms4avgULTMT = permute(h5read(loadName4avgULTMT, ['/', h5GroupName]), ...
            [2, 1, 3]);
        movieTEMP = implay(single(rescale(frms4avgULTMT)), 100);
        movieTEMP.Parent.Position = ...
            [30 + 2*FOVpixDim(1), basePixSize - FOVpixDim(2) - 134, ...
            FOVpixDim + [0, 24]].*scaleFactor;
        clear frms4avgULTMT;
    end

    if deleteH5ON
        loadName4avgULTMT = fullfile(savePath, ['frms4avgULTMT', suffix0, '.', outputFormat]);
        delete(loadName4avgULTMT);
    end

    if closeFigureON && max(stepNum) ~= 3
        close all;
    end
end

%% Step 4: The 2nd Motion Correction to Top 50 Frames in each Trial (rigid)
if motCrct2ndON
    loadNameAvgInfo_Rigid = fullfile(savePath, 'frmsAvgInfo_Rigid.mat');
    load(loadNameAvgInfo_Rigid, 'frmsAvgdULTMT');
    loadName4rawData = fullfile(loadPathH5, ['rawData-', Date, '-', FOV, '.', outputFormat]);
    
    load(fullfile(loadPathBHV, ['bhvData-', Date, '-', FOV, '_', TaskInfo, '.mat']), 'dataset');
    trlsResampID = find(~any(isnan(dataset.eventFrms{:, [1 end]}), 2));
    frmsResampIDstt = dataset.eventFrms.fixON(trlsResampID).' - 1;
    frmsResampID = reshape(frmsResampIDstt + (1:frmsResampNumByTrl).', 1, []);
%     frmsResampNum = length(frmsResampID);
    trlsResampNum = length(trlsResampID);

    saveNameMOTcrctd_Rigid = fullfile(savePath, ['frmsMOTcrctd-', ...
        num2str(frmsResampNumByTrl), 'avgd', suffix0, '.', outputFormat]);
    [DIYparams_Rigid, ~, ~, ~] = ...
        normcorre_batch_Rigid_4GPU_outputDIY(loadName4rawData, ...
        NoRMCorreSetParms(options4motCrct, 'output_filename',  saveNameMOTcrctd_Rigid, ...
        'init_batch', batchSize4load2nd(1), 'bin_width', batchSize4load2nd(2), ...
        'mem_batch_size', frmsResampNumByTrl, 'print_msg', 100), ...
        frmsAvgdULTMT.', MaxCPUcoreNum4SPMD1, 1, 0, 1, 1, ...
        max(round(trlsResampNum*frmsResampNumByTrl/(MaxFrmsNum4dAND16*MaxCPUcoreNum4SPMD1)), 1), ...
        1, frmsResampID);

    trlsResamp = double(h5read(saveNameMOTcrctd_Rigid, ['/', h5GroupName]));

    corrTrl2Trl = corr(reshape(trlsResamp, FOVpixNum, []));
    corrTrl2Avg = DIYparams_Rigid.ZNCC;
    frmsTrlAvgd = DIYparams_Rigid.templateOutput;

    clear trlsResamp;
    
    MaxCompNum = 8; % change to 8
    MaxRepeatNum = 100;
    optns = statset('Display', 'off', 'MaxIter', 1000, 'TolFun', fitgmdistTolFun(1));
    fitdGMM = cell(MaxRepeatNum, MaxCompNum);
    for n = 1:MaxCompNum
        fprintf(['------------------------------\n', ...
            'GMMs compNum = %d \n'], n);
        tic;
        parfor m = 1:MaxRepeatNum
            warning('off');
            fitdGMM{m, n} = fitgmdist(atanh(corrTrl2Avg), n, 'Options', optns, 'Start', 'plus', ...
                'RegularizationValue', fitgmdistRegVal(1));
            warning('on');
        end
        toc;
    end

    frmsTrlAvgd = frmsTrlAvgd.';
    save(loadNameAvgInfo_Rigid, 'corrTrl2Trl', 'corrTrl2Avg', ...
        'frmsTrlAvgd', 'fitdGMM', 'trlsResampID', 'trlsResampNum', 'frmsResampIDstt', 'DIYparams_Rigid', '-append');

    if any(stepNum >= 6)
        delete(gcp('nocreate'));
        parpool("local", MaxCPUcoreNum);
    end

    if closeFigureON && max(stepNum) ~= 4
        close all;
    end
end

%% Step 5: Check the Results of the 2nd Motion Correction
if motCrct2ndCheckON
    loadNameAvgInfo_Rigid = fullfile(savePath, 'frmsAvgInfo_Rigid.mat');
    load(loadNameAvgInfo_Rigid, 'corrTrl2Trl', 'corrTrl2Avg', 'frmsTrlAvgd', 'frmsAvgdULTMT', 'trlsResampNum', 'fitdGMM');
    [MaxRepeatNum, MaxCompNum] = size(fitdGMM);

    corrTrl2Avg_atanh = atanh(corrTrl2Avg);
    x = linspace(min(corrTrl2Avg_atanh), max(corrTrl2Avg_atanh), 2^8).';
    ksd = ksdensity(corrTrl2Avg_atanh, x, 'Bandwidth', 0.03);

    alpha = 0.05;
    thresBIC = 1/5;

    k1 = round(0.5*MaxRepeatNum);
    k2 = ceil(0.1*k1);

    if strcmp(MSCmethod, 'RSS')
        MSC = cellfun(@(GMM) sum((pdf(GMM, x) - ksd).^2), fitdGMM);
    elseif strcmp(MSCmethod, 'BIC')
        MSC = cellfun(@(x) x.BIC, fitdGMM);
    end
    [MSC, minkMSCid] = mink(MSC, k1, 1);
    fitdGMM = fitdGMM(minkMSCid + (0:MaxCompNum - 1).*MaxRepeatNum);

    if isempty(CompNum2use0)
        if strcmp(MSCmethod, 'RSS')
%             CompNum2use = find((mean(MSC, 1) - min(mean(MSC, 1))) < ...
%                 thresBIC*(max(mean(MSC, 1)) - min(mean(MSC, 1))), 1);
            CompNum2use = find((median(MSC, 1) - min(median(MSC, 1))) < ...
                thresBIC*(max(median(MSC, 1)) - min(median(MSC, 1))), 1);
        elseif strcmp(MSCmethod, 'BIC')
            [~, CompNum2use] = min(mean(MSC, 1));
%             [~, CompNum2use] = min(median(MSC, 1));
        end
    else
        CompNum2use = CompNum2use0;
    end

    [~, ~, s] = anova1(MSC, [], 'off');
    figure;
    set(gcf, 'Position', ...
        [20 + 560, 50, ...
        315, 315].*scaleFactor);
    c = multcompare(s, 'Alpha', alpha);
    xLim = xlim;
    hold on;
    tempY = linspace(1, MaxCompNum, 2^8);
%     tempPlot(1) = plot(expfit(tempY).*temp, flip(tempY), 'Color', 'r');
    tempPlot = plot(xLim, MaxCompNum - CompNum2use + ones(1, 2), 'Color', '#77AC30');
    hold off;
    uistack(tempPlot, 'bottom');
    xlim(xLim);


%     [~, minBICid] = min(BIC(:, CompNum2use), [], 1);
%     fitdGMM2use = fitdGMM{minBICid, CompNum2use};
%     fitdParams2use = [fitdGMM2use.mu, fitdGMM2use.Sigma(:), fitdGMM2use.ComponentProportion.'];

    fitdParams = reshape(cell2mat(cellfun(@(x) sortrows( ...
        [x.mu, x.Sigma(:), x.ComponentProportion.'], 1), ...
        fitdGMM(1:k2, CompNum2use).', ...
        'UniformOutput', false)), CompNum2use, 3, []); % Top k2
%     fitdParams2use = median(fitdParams, 3);
    fitdParams2use = mean(fitdParams, 3);

    optns = statset('Display', 'off', 'MaxIter', 1000, 'TolFun', fitgmdistTolFun(1));
    initParams = struct('mu', fitdParams2use(:, 1), ...
        'Sigma', reshape(fitdParams2use(:, 2), 1, 1, []), ...
        'ComponentProportion', fitdParams2use(:, 3));
    fitdGMM2use = fitgmdist(corrTrl2Avg_atanh, CompNum2use, 'Options', optns, 'Start', initParams, ...
        'RegularizationValue', fitgmdistRegVal(1));
    fitdParams2use = sortrows([fitdGMM2use.mu, fitdGMM2use.Sigma(:), ...
        fitdGMM2use.ComponentProportion.'], 1);

%     fitdGMM2use = gmdistribution(fitdParams2use(:, 1), ...
%         reshape(fitdParams2use(:, 2), 1, 1, []), fitdParams2use(:, 3));

    fitdParams2use(:, 2) = sqrt(fitdParams2use(:, 2));
    tempID = fitdParams2use(:, 3) > 0.5*alpha;
    threshAll = tanh(min(fitdParams2use(tempID, 1) - 3*fitdParams2use(tempID, 2)));

    fig4save{1} = figure;
    set(gcf, 'Position', ...
        [30 + 560 + 315, 50, ...
        560, 315].*scaleFactor);
    histogram(corrTrl2Avg_atanh, 50, 'Normalization', 'pdf', 'DisplayStyle', 'bar', ...
        'EdgeColor', 'none', 'FaceColor', 0.75*[1 1 1]);
    hold on;
    clrID = round(linspace(1, clrN, CompNum2use));
    for n = 1:CompNum2use
        patch([x(1); x; x(end)], [0; fitdParams2use(n, 3).*pdf('normal', x, fitdParams2use(n, 1), fitdParams2use(n, 2)); 0], ...
            clr(clrID(n), :), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end

    tempPlot = plot(x, ksd, 'b', 'LineWidth', 6);
    tempPlot.Color(4) = 0.2;

    plot(x, pdf(fitdGMM2use, x), 'k');
    outlierID = find(corrTrl2Avg < threshAll);
    tempN = length(outlierID);
    save(loadNameAvgInfo_Rigid, 'outlierID', 'threshAll', 'fitdGMM2use', '-append');

    plot(repmat(atanh(threshAll), 1, 2), ylim, 'r');
    hold off;
    text(atanh(threshAll) + diff(xlim)*0.01, mean(ylim), {['\theta = ', num2str(threshAll, '%.4f')]; ...
        ['{\itN}(ZNCC < \theta) = ', num2str(tempN, '%d')]}, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
    ylabel('Probability Density');
    xlabel('ZNCC ( a.k.a. Pearson''s {\itr} )');
    xTick = linspace(min(corrTrl2Avg), max(corrTrl2Avg), 10);
    set(gca, 'XTick', atanh(xTick), 'XTickLabel', num2str(xTick.', '%.3f'), 'XTickLabelRotation', 45);

    fig4save{2} = figure;
    set(gcf, 'Position', ...
        [10, 50, ...
        560, 315].*scaleFactor);
    xLim = [0, trlsResampNum + 1];
    t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    plot(corrTrl2Avg, 'Color', '#0072BD');

    hold on;
    plot(xlim, repmat(threshAll, 1, 2), 'r');
    hold off;
    title('Chronological Order');
    xlim(xLim);

    nexttile;
    semilogx(sort(corrTrl2Avg), 'Color', '#0072BD');
    yLim = ylim;
    hold on;
    semilogx([xLim(1) + 1, repmat(tempN + 0.5, 1, 2)], [repmat(threshAll, 1, 2), yLim(1)], 'r');
    hold off;
    title('ZNCC Ascending Order');
    xlim(xLim + [1 0]);
    xLim = xlim;
    yLim = ylim;
    text(10^(0.025*diff(log10(xLim))), yLim(1) + 0.9*diff(yLim), {['\theta = ', num2str(threshAll, '%.4f')]; ...
        ['{\itN}(ZNCC < \theta) = ', num2str(tempN, '%d')]}, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    xlabel(t, 'trial number');
    ylabel(t, 'ZNCC ( a.k.a. Pearson''s {\itr} )');


    fig4save{3} = figure;
    set(gcf, 'Position', ...
        [40 + 2*560 + 315, 50, ...
        400, 315].*scaleFactor);
    imagesc(corrTrl2Trl);
    cb = colorbar;
    title(cb, 'ZNCC', 'FontAngle', 'normal', 'VerticalAlignment', 'baseline');
    colormap(flip(gray));
    axis tight;
    axis equal;

    for n = 1:length(fig4save)
        fig4save{n}.Renderer = 'painters';
    end
    if saveImageON
        saveas(fig4save{1}, fullfile(savePath, 'ZNCC_trlPDF_Rigid.svg'));
        saveas(fig4save{2}, fullfile(savePath, 'ZNCC_trl_Rigid.svg'));
        saveas(fig4save{3}, fullfile(savePath, 'ZNCC_trl2trl_Rigid.svg'));
        imwrite(ORGNLgreen(frmsTrlAvgd), ...
            fullfile(savePath, 'frmsTrlAvgd_ORGNLgreen_uint8_Rigid.tif'), ...
            'Compression', 'none');
        imwrite(CLAHEslight(frmsTrlAvgd), ...
            fullfile(savePath, 'frmsTrlAvgd_CLAHEslight_uint8_Rigid.tif'), ...
            'Compression', 'none');
        imwrite(CLAHEstrong(frmsTrlAvgd), ...
            fullfile(savePath, 'frmsTrlAvgd_CLAHEstrong_uint8_Rigid.tif'), ...
            'Compression', 'none');
    end

    figure;
    set(gca, 'Position', [0, 0, 1, 1]);
    set(gcf, 'Position', ...
        [10, basePixSize - FOVpixDim(2) - 110, ...
        FOVpixDim].*scaleFactor);
    imshow(CLAHEslight(frmsAvgdULTMT));

    figure;
    set(gca, 'Position', [0, 0, 1, 1]);
    set(gcf, 'Position', ...
        [20 + FOVpixDim(1), basePixSize - FOVpixDim(2) - 110, ...
        FOVpixDim].*scaleFactor);
    imshow(CLAHEslight(frmsTrlAvgd));

    if playMovieON(2)
        loadNameMOTcrctd_Rigid = fullfile(savePath, ['frmsMOTcrctd-', ...
            num2str(frmsResampNumByTrl), 'avgd', suffix0, '.', outputFormat]);
        trlsResamp = permute(h5read(loadNameMOTcrctd_Rigid, ['/', h5GroupName]), ...
            [2, 1, 3]);
        movieTEMP = implay(single(rescale(trlsResamp)), 10);
        movieTEMP.Parent.Position = ...
            [30 + 2*FOVpixDim(1), basePixSize - FOVpixDim(2) - 134, ...
            FOVpixDim + [0, 24]].*scaleFactor;
        clear frms4resamp;
    end

    if closeFigureON && max(stepNum) ~= 5
        close all;
    end
end


%% Step 6: The 3rd Motion Correction and Interpolation for Generating Non-rigid Displacement Fields
if motCrct3rdON
    loadNameMOTcrctd_Rigid = fullfile(savePath, ['frmsMOTcrctd-', ...
        num2str(frmsResampNumByTrl), 'avgd', suffix0, '.', outputFormat]);
    trlsResamp = h5read(loadNameMOTcrctd_Rigid, ['/', h5GroupName]);
    loadNameAvgInfo_Rigid = fullfile(savePath, 'frmsAvgInfo_Rigid.mat');
    load(loadNameAvgInfo_Rigid, 'outlierID', 'trlsResampNum', 'frmsResampIDstt');
    if ~isempty(outlierID)
        frmsResampIDstt(outlierID) = [];
        trlsResamp(:, :, outlierID) = [];
    end

    [trlsBinNum, trlsBinID, ~] = autoGroupOverlap(trlsResampNum - length(outlierID), trlsBinSize, trlsBinSize - trlsBinDist);
    trlsBin = NaN([FOVpixDim, trlsBinNum]);
    for i = 1:trlsBinNum
        trlsBin(:, :, i) = sum(trlsResamp(:, :, trlsBinID(1, i):trlsBinID(2, i)), 3);
    end
    trlsBin = trlsBin./trlsBinSize;
    trlsBinfrmsID = NaN(1, trlsBinNum);
    for i = 1:trlsBinNum
        trlsBinfrmsID(i) = mean(frmsResampIDstt(trlsBinID(1, i):trlsBinID(2, i)));
    end
    trlsBinfrmsID = round(trlsBinfrmsID + (frmsResampNumByTrl + 1)./2);

    if CLAHEutmostON
        tic;
        parfor i = 1:trlsBinNum
            trlsBin(:, :, i) = CLAHEutmost(trlsBin(:, :, i));
        end
        toc;
    end

    patchSize = [64, 64];
    overlapSize = [48, 48];
    patchDist = patchSize - overlapSize;
    [patchDim(1), patchID_X, ~] = autoGroupOverlap(FOVpixDim(1), patchSize(1), overlapSize(1));
    [patchDim(2), patchID_Y, ~] = autoGroupOverlap(FOVpixDim(2), patchSize(2), overlapSize(2));
    patchNum = prod(patchDim);
    patchID = NaN(patchDim(1), patchDim(2), 2, 2);
    [patchID(:, :, 1, 1), patchID(:, :, 2, 1)] = ndgrid(patchID_X(1, :), patchID_Y(1, :));
    [patchID(:, :, 1, 2), patchID(:, :, 2, 2)] = ndgrid(patchID_X(2, :), patchID_Y(2, :));
    patchIDcntr = mean(patchID, 4);
    patchID = reshape(permute(patchID, [4 3 1 2]), 4, []);

    trlsBinAvgd = mean(trlsBin, 3);
    figure;
    set(gca, 'Position', [0, 0, 1, 1]);
    set(gcf, 'Position', ...
        [10, basePixSize - FOVpixDim(2) - 110, ...
        FOVpixDim].*scaleFactor);
    imshow(trlsBinAvgd.');

    tmplts4DF = NaN([patchSize, patchNum]);
    trlsBinPatch = NaN([patchSize, patchNum*trlsBinNum], 'single');
    for i = 1:patchNum
        trlsBinPatch(:, :, (i - 1)*trlsBinNum + (1:trlsBinNum)) = single(trlsBin(patchID(1, i):patchID(2, i), patchID(3, i):patchID(4, i), :));
        tmplts4DF(:, :, i) = trlsBinAvgd(patchID(1, i):patchID(2, i), patchID(3, i):patchID(4, i));
    end
    saveName4DF = fullfile(savePath, ['frms4DF', suffix1, '.', outputFormat]);
    h5create(saveName4DF, ['/', h5GroupName], [patchSize, Inf], ...
        'Chunksize', [patchSize, trlsBinNum], 'Datatype', 'single');
    h5write(saveName4DF, ['/', h5GroupName], trlsBinPatch, ...
        [1, 1, 1], [patchSize, patchNum*trlsBinNum]);
    clear trlsBinPatch;

    [DIYparams_nonRigid, ~, ~, ~] = ...
        normcorre_batch_nonRigid_4GPU(saveName4DF, ...
        NoRMCorreSetParms(options4motCrct, 'd1', patchSize(1), 'd2', patchSize(2), ...
        'grid_size', patchSize, 'init_batch', batchSize4load3rd(1), ...
        'bin_width', batchSize4load3rd(2), 'mem_batch_size', trlsBinNum, 'print_msg', 50), ...
        tmplts4DF, MaxCPUcoreNum4SPMD2, 1, 0, iterNum4step6, 0, 0);

    delete(saveName4DF);
    saveNameAvgInfo_nonRigid = fullfile(savePath, 'frmsAvgInfo_nonRigid.mat');
    save(saveNameAvgInfo_nonRigid, 'trlsBin', 'trlsBinID', 'trlsBinNum', 'trlsBinfrmsID', ...
        'patchSize', 'overlapSize', 'patchDist', 'patchDim', 'patchID_X', 'patchID_Y', 'patchIDcntr', ...
        'DIYparams_nonRigid');
    clear trlsBin;

    if any(stepNum >= 8)
        delete(gcp('nocreate'));
        parpool("local", MaxCPUcoreNum);
    end

    if closeFigureON && max(stepNum) ~= 6
        close all;
    end
end


%% Step 7: Check the Results of the 3rd Motion Correction
if motCrct3rdCheckON
    loadNameAvgInfo_nonRigid = fullfile(savePath, 'frmsAvgInfo_nonRigid.mat');
    load(loadNameAvgInfo_nonRigid, 'trlsBinNum', 'trlsBinfrmsID', ...
        'patchSize', 'patchDist', 'patchDim', 'patchID_X', 'patchID_Y', 'patchIDcntr', ...
        'DIYparams_nonRigid');

    XgridVec = mean(patchID_X);
    YgridVec = mean(patchID_Y);
    ZgridVec = trlsBinfrmsID;

    %     nonRigidDF = permute(reshape(DIYparams.shifts(:, 1:2), trlsBinNum, patchDim(1), patchDim(2), 2) + ...
    %         reshape(DIYparams.shifts2(:, 1:2), 1, patchDim(1), patchDim(2), 2), [3 2 4 1]); % [Y X dim trl]
    nonRigidDF = smoothdata( ...
        permute( ...
        reshape(DIYparams_nonRigid.shifts(:, 1:2), trlsBinNum, patchDim(1), patchDim(2), 2), ...
        [2 3 1 4]), ...
        3, smoothMethod, smoothWindow); % [X Y trl dim]

    maxPix = max(FOVpixDim./(FOVpixDim - patchSize))*max(abs(nonRigidDF), [], 'all');
    sclFctr = min(patchDist)./maxPix;
    FOVpixDimADD = FOVpixDim + 2*max(patchDist);

    figDF = figure;
    set(gcf, 'Position', [10 + 200, 50, FOVpixDimADD]);
    patch('XData', [1, FOVpixDim(1), FOVpixDim(1), 1] + 0.5*[-1 1 1 -1], ...
        'YData', [1, 1, FOVpixDim(2), FOVpixDim(2)] + 0.5*[-1 -1 1 1], ...
        'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 1);
    trlsBinDF = zeros(FOVpixDimADD(2), FOVpixDimADD(1), 3, trlsBinNum, 'uint8');
    box on;
    hold on;
    for trl = 1:trlsBinNum
        tempDF = quiver(XgridVec, YgridVec, ...
            sclFctr*nonRigidDF(:, :, trl, 1).', sclFctr*nonRigidDF(:, :, trl, 2).', 0, ...
            'Color', [0.3 0.6 0.9], 'LineWidth', 1);
        
        axis([1, FOVpixDim(1), 1, FOVpixDim(2)] + (0.5 + max(patchDist))*[-1 1 -1 1]);
        set(gca, 'XTick', XgridVec, 'YTick', YgridVec, 'XGrid', 'on', 'YGrid', 'on');
        set(gca, 'Position', [0, 0, 1, 1], 'YDir', 'reverse');

        trlsBinDF(:, :, :, trl) = frame2im(getframe(figDF));
        delete(tempDF);
    end
    hold off;

    movieTEMP = implay(trlsBinDF, 60);
    movieTEMP.Parent.Position = ...
        [30 + 1*FOVpixDimADD(1), basePixSize - FOVpixDimADD(2) - 134, ...
        FOVpixDimADD + [0, 24]].*scaleFactor;

    if ~strcmp(extrapMethod, 'linearDIY')
        XYZgridVecs = {XgridVec, YgridVec, ZgridVec};
        interpX = griddedInterpolant(XYZgridVecs, nonRigidDF(:, :, :, 1), 'nearest', extrapMethod);
        interpY = griddedInterpolant(XYZgridVecs, nonRigidDF(:, :, :, 2), 'nearest', extrapMethod);
    end

    thk = round([(YgridVec(1) - 0.5)./patchDist(2), ...
        (FOVpixDim(1) + 0.5 - XgridVec(end))./patchDist(1), ...
        (FOVpixDim(2) + 0.5 - YgridVec(end))./patchDist(2), ...
        (XgridVec(1) - 0.5)./patchDist(1)]);
    XgridVec = [linspace(0.5, XgridVec(1), 1 + thk(4)), XgridVec(2:end - 1), ...
        linspace(XgridVec(end), FOVpixDim(1) + 0.5, 1 + thk(2))];
    YgridVec = [linspace(0.5, YgridVec(1), 1 + thk(1)), YgridVec(2:end - 1), ...
        linspace(YgridVec(end), FOVpixDim(2) + 0.5, 1 + thk(3))];

    if strcmp(extrapMethod, 'linearDIY')
        XgridNum = length(XgridVec);
        YgridNum = length(YgridVec);
        XYgridNum = XgridNum*YgridNum;

        tempXYin = reshape(patchIDcntr, [], 2);
        tempV = reshape(nonRigidDF, [], trlsBinNum, 2);

        edgeIDin = MatrixEdgeIndex([patchDim(1), patchDim(2)], ones(1, 4));
        edgeIDout = MatrixEdgeIndex([XgridNum, YgridNum], thk);
        cntrIDout = find(~ismember(1:XYgridNum, edgeIDout));
        nonRigidDFpadded = NaN(XYgridNum, trlsBinNum, 2, 'single');
        nonRigidDFpadded(cntrIDout, :, :) = tempV;

        tempXYin = tempXYin(edgeIDin, :);
        tempV = double(tempV(edgeIDin, :, :));

        tempXYout = NaN(XgridNum, YgridNum, 2);
        [tempXYout(:, :, 1), tempXYout(:, :, 2)] = ndgrid(XgridVec, YgridVec);
        tempXYout = reshape(tempXYout, [], 2);
        tempXYout = tempXYout(edgeIDout, :);

        [extrapID, W] = LinearExtrap(tempXYin, tempXYout, (1 + FOVpixDim)./2);

        for trl = 1:trlsBinNum
            tempV2 = tempV(:, trl, 1);
            nonRigidDFpadded(edgeIDout, trl, 1) = sum(tempV2(extrapID).*W, 2);
            tempV2 = tempV(:, trl, 2);
            nonRigidDFpadded(edgeIDout, trl, 2) = sum(tempV2(extrapID).*W, 2);
        end
        nonRigidDFpadded = reshape(nonRigidDFpadded, XgridNum, YgridNum, trlsBinNum, 2);
    else
        XYZgridVecs = {XgridVec, YgridVec, ZgridVec};
        nonRigidDFpadded = cat(4, interpX(XYZgridVecs), interpY(XYZgridVecs));
    end

    trlsBinDFpadded = zeros(FOVpixDimADD(2), FOVpixDimADD(1), 3, trlsBinNum, 'uint8');
    hold on;
    for trl = 1:trlsBinNum
        tempDF = quiver(XgridVec, YgridVec, ...
            sclFctr*nonRigidDFpadded(:, :, trl, 1).', sclFctr*nonRigidDFpadded(:, :, trl, 2).', 0, ...
            'Color', [0.3 0.6 0.9], 'LineWidth', 1);

        axis([1, FOVpixDim(1), 1, FOVpixDim(2)] + (0.5 + max(patchDist))*[-1 1 -1 1]);
        set(gca, 'XTick', XgridVec, 'YTick', YgridVec, 'XGrid', 'on', 'YGrid', 'on');
        set(gca, 'Position', [0, 0, 1, 1], 'YDir', 'reverse');

        trlsBinDFpadded(:, :, :, trl) = frame2im(getframe(figDF));
        delete(tempDF);
    end
    hold off;
    close(figDF);

    movieTEMP = implay(trlsBinDFpadded, 60);
    movieTEMP.Parent.Position = ...
        [30 + 2*FOVpixDimADD(1), basePixSize - FOVpixDimADD(2) - 134, ...
        FOVpixDimADD + [0, 24]].*scaleFactor;

    save(loadNameAvgInfo_nonRigid, 'XgridVec', 'YgridVec', 'nonRigidDFpadded', '-append');

    if closeFigureON && max(stepNum) ~= 7
        close all;
    end
end


%% Step 8: The 4th Motion Correction to All Frames (non-rigid)
if motCrct4thON
    loadNameAvgInfo_Rigid = fullfile(savePath, 'frmsAvgInfo_Rigid.mat');
    load(loadNameAvgInfo_Rigid, 'DIYparams_Rigid', 'trlsResampID', 'trlsResampNum', 'frmsResampIDstt');

    DF.Rigid = DIYparams_Rigid.shifts(:, 1:2);

    loadNameAvgInfo_nonRigid = fullfile(savePath, 'frmsAvgInfo_nonRigid.mat');
    load(loadNameAvgInfo_nonRigid, 'trlsBinfrmsID', ...
        'XgridVec', 'YgridVec', 'nonRigidDFpadded');

    XYgridVecs = {XgridVec, YgridVec};
    XYZgridVecs = [XYgridVecs, {trlsBinfrmsID}];
    interpX = griddedInterpolant(XYZgridVecs, nonRigidDFpadded(:, :, :, 1), 'makima', 'nearest');
    interpY = griddedInterpolant(XYZgridVecs, nonRigidDFpadded(:, :, :, 2), 'makima', 'nearest');

    XYZgridVecsFull = [XYgridVecs, {1:frmsFullNum}];
    tic;
    DF.nonRigid = cat(4, interpX(XYZgridVecsFull), interpY(XYZgridVecsFull));
    toc;

    saveName = fullfile(savePath, 'frms4ROIdtct', 'plane0', ...
        ['data.', outputFormatFinal]);
    if ~exist(fullfile(savePath, 'frms4ROIdtct', 'plane0'), 'dir')
        mkdir(fullfile(savePath, 'frms4ROIdtct', 'plane0'));
    end
    loadName4rawData = fullfile(loadPathH5, ['rawData-', Date, '-', FOV, '.', outputFormat]);

    optionsTEMP = NoRMCorreSetParms(options4motCrct, ...
        'output_type', outputFormatFinal, 'output_filename', saveName, ...
        'init_batch', batchSize4load4th(1), 'bin_width', batchSize4load4th(2), ...
        'mem_batch_size', chunkSize4save2nd, 'print_msg', 10);
    optionsTEMP.final_chunk_size = 100;
    optionsTEMP.output_data_type = outputDataTypeFinal;
    [~, frmsAllAvgd, ~] = ...
        normcorre_batch_nonRigid_4GPU_applyShift(loadName4rawData, ...
        optionsTEMP, MaxCPUcoreNum4SPMD3, 0, 0, ...
        max(round(frmsFullNum/(MaxFrmsNum4uint16*MaxCPUcoreNum4SPMD3)), 1), ... % max(round(frmsFullNum/(MaxFrmsNum4dAND16*MaxCPUcoreNum4SPMD3)), 1)
        1, DF, XYgridVecs);

    frmsTempSum = NaN([FOVpixDim, trlsResampNum]);

    if strcmp(outputFormatFinal, 'bin')
        parfor n = 1:trlsResampNum
            frmsTemp = read_bin_file_DIY(saveName, frmsResampIDstt(n) + 1, frmsResampNumByTrl, ...
                FOVpixDim, bitSize, outputDataTypeFinal, inputDataType, 3);
            frmsTempSum(:, :, n) = sum(frmsTemp, 3);
        end
    elseif strcmp(outputFormatFinal, 'h5')
        parfor n = 1:trlsResampNum
            frmsTemp = h5read(saveName, ['/', h5GroupName], ...
                [1, 1, frmsResampIDstt(n) + 1], [FOVpixDim, frmsResampNumByTrl]);
            frmsTempSum(:, :, n) = sum(frmsTemp, 3);
        end
    end

    frmsTrlAvgd = sum(frmsTempSum, 3)./(trlsResampNum*frmsResampNumByTrl);
    frmsTempAvg = frmsTempSum./frmsResampNumByTrl;
    clear frmsTemp frmsTempSum;

    saveNameMOTcrctd_nonRigid = fullfile(savePath, ['frmsMOTcrctd-', ...
        num2str(frmsResampNumByTrl), 'avgd', suffix1, '.', outputFormat]);
    if exist(saveNameMOTcrctd_nonRigid, 'file')
        delete(saveNameMOTcrctd_nonRigid);
    end
    h5create(saveNameMOTcrctd_nonRigid, ['/', h5GroupName], [FOVpixDim, Inf], ...
        'ChunkSize', [FOVpixDim, trlsResampNum], 'Datatype', inputDataType);
    h5write(saveNameMOTcrctd_nonRigid, ['/', h5GroupName], ...
        cast(frmsTempAvg, inputDataType), [1, 1, 1], [FOVpixDim, trlsResampNum]);

    corrTrl2Trl = corr(reshape(frmsTempAvg, FOVpixNum, []));
    corrTrl2Avg = corr(reshape(frmsTempAvg, FOVpixNum, []), frmsTrlAvgd(:));
    clear frmsTempAvg;

    MaxCompNum = 8; % change to 8
    MaxRepeatNum = 100;
    optns = statset('Display', 'off', 'MaxIter', 1000, 'TolFun', fitgmdistTolFun(2));
    fitdGMM = cell(MaxRepeatNum, MaxCompNum);
    for n = 1:MaxCompNum
        fprintf(['------------------------------\n', ...
        'GMMs compNum = %d \n'], n);
        tic;
        parfor m = 1:MaxRepeatNum
            warning('off');
            fitdGMM{m, n} = fitgmdist(atanh(corrTrl2Avg), n, 'Options', optns, 'Start', 'plus', ...
                'RegularizationValue', fitgmdistRegVal(2));
            warning('on');
        end
        toc;
    end

    frmsTrlAvgd = frmsTrlAvgd.';
    frmsAllAvgd = frmsAllAvgd.';
    saveNameAvgInfo_Final = fullfile(savePath, 'frmsAvgInfo_Final.mat');
    save(saveNameAvgInfo_Final, 'corrTrl2Trl', 'corrTrl2Avg', ...
        'frmsTrlAvgd', 'frmsAllAvgd', 'fitdGMM', 'trlsResampID', 'trlsResampNum');

    if closeFigureON && max(stepNum) ~= 8
        close all;
    end
end


%% Step 9: Check the Results of the 4th Motion Correction
if motCrct4thCheckON
    loadNameAvgInfo_Rigid = fullfile(savePath, 'frmsAvgInfo_Rigid.mat');
    tempRigid = load(loadNameAvgInfo_Rigid, 'corrTrl2Avg', 'frmsTrlAvgd');

    loadNameAvgInfo_Final = fullfile(savePath, 'frmsAvgInfo_Final.mat');
    load(loadNameAvgInfo_Final, 'corrTrl2Trl', 'corrTrl2Avg', 'frmsTrlAvgd', 'trlsResampNum', 'fitdGMM');
    [MaxRepeatNum, MaxCompNum] = size(fitdGMM);

    corrTrl2Avg_atanh = atanh(corrTrl2Avg);
    x = linspace(min(corrTrl2Avg_atanh), max(corrTrl2Avg_atanh), 2^8).';
    ksd = ksdensity(corrTrl2Avg_atanh, x, 'Bandwidth', 0.03);

    alpha = 0.05;
    thresBIC = 1/5;

    k1 = round(0.5*MaxRepeatNum);
    k2 = ceil(0.1*k1);

    if strcmp(MSCmethod, 'RSS')
        MSC = cellfun(@(GMM) sum((pdf(GMM, x) - ksd).^2), fitdGMM);
    elseif strcmp(MSCmethod, 'BIC')
        MSC = cellfun(@(x) x.BIC, fitdGMM);
    end
    [MSC, minkMSCid] = mink(MSC, k1, 1);
    fitdGMM = fitdGMM(minkMSCid + (0:MaxCompNum - 1).*MaxRepeatNum);

    if isempty(CompNum2use0)
        if strcmp(MSCmethod, 'RSS')
%             CompNum2use = find((mean(MSC, 1) - min(mean(MSC, 1))) < ...
%                 thresBIC*(max(mean(MSC, 1)) - min(mean(MSC, 1))), 1);
            CompNum2use = find((median(MSC, 1) - min(median(MSC, 1))) < ...
                thresBIC*(max(median(MSC, 1)) - min(median(MSC, 1))), 1);
        elseif strcmp(MSCmethod, 'BIC')
            [~, CompNum2use] = min(mean(MSC, 1));
%             [~, CompNum2use] = min(median(MSC, 1));
        end
    else
        CompNum2use = CompNum2use0;
    end

    [~, ~, s] = anova1(MSC, [], 'off');
    figure;
    set(gcf, 'Position', ...
        [20 + 560, 50, ...
        315, 315].*scaleFactor);
    c = multcompare(s, 'Alpha', alpha);
    xLim = xlim;
    hold on;
    tempY = linspace(1, MaxCompNum, 2^8);
%     tempPlot(1) = plot(expfit(tempY).*temp, flip(tempY), 'Color', 'r');
    tempPlot = plot(xLim, MaxCompNum - CompNum2use + ones(1, 2), 'Color', '#77AC30');
    hold off;
    uistack(tempPlot, 'bottom');
    xlim(xLim);


%     [~, minBICid] = min(BIC(:, CompNum2use), [], 1);
%     fitdGMM2use = fitdGMM{minBICid, CompNum2use};
%     fitdParams2use = [fitdGMM2use.mu, fitdGMM2use.Sigma(:), fitdGMM2use.ComponentProportion.'];

    fitdParams = reshape(cell2mat(cellfun(@(x) sortrows( ...
        [x.mu, x.Sigma(:), x.ComponentProportion.'], 1), ...
        fitdGMM(1:k2, CompNum2use).', ...
        'UniformOutput', false)), CompNum2use, 3, []); % Top k2
%     fitdParams2use = median(fitdParams, 3);
    fitdParams2use = mean(fitdParams, 3);

    optns = statset('Display', 'off', 'MaxIter', 1000, 'TolFun', fitgmdistTolFun(2));
    initParams = struct('mu', fitdParams2use(:, 1), ...
        'Sigma', reshape(fitdParams2use(:, 2), 1, 1, []), ...
        'ComponentProportion', fitdParams2use(:, 3));
    fitdGMM2use = fitgmdist(corrTrl2Avg_atanh, CompNum2use, 'Options', optns, 'Start', initParams, ...
        'RegularizationValue', fitgmdistRegVal(2));
    fitdParams2use = sortrows([fitdGMM2use.mu, fitdGMM2use.Sigma(:), ...
        fitdGMM2use.ComponentProportion.'], 1);

%     fitdGMM2use = gmdistribution(fitdParams2use(:, 1), ...
%         reshape(fitdParams2use(:, 2), 1, 1, []), fitdParams2use(:, 3));

    fitdParams2use(:, 2) = sqrt(fitdParams2use(:, 2));
    tempID = fitdParams2use(:, 3) > 0.5*alpha;
    threshAll = tanh(min(fitdParams2use(tempID, 1) - 3*fitdParams2use(tempID, 2)));

    fig4save{1} = figure;
    set(gcf, 'Position', ...
        [30 + 560 + 315, 50, ...
        560, 315].*scaleFactor);
    histogram(corrTrl2Avg_atanh, 75, 'Normalization', 'pdf', 'DisplayStyle', 'bar', ...
        'EdgeColor', 'none', 'FaceColor', 0.75*[1 1 1]);
    hold on;
    clrID = round(linspace(1, clrN, CompNum2use));
    for n = 1:CompNum2use
        patch([x(1); x; x(end)], [0; fitdParams2use(n, 3).*pdf('normal', x, fitdParams2use(n, 1), fitdParams2use(n, 2)); 0], ...
            clr(clrID(n), :), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end

    tempPlot = plot(x, ksd, 'b', 'LineWidth', 6);
    tempPlot.Color(4) = 0.2;

    plot(x, pdf(fitdGMM2use, x), 'k');
    outlierID = find(corrTrl2Avg < threshAll);
    tempN = length(outlierID);
    save(loadNameAvgInfo_Final, 'outlierID', 'threshAll', 'fitdGMM2use', '-append');

    plot(repmat(atanh(threshAll), 1, 2), ylim, 'r');
    hold off;
    text(atanh(threshAll) + diff(xlim)*0.01, mean(ylim), {['\theta = ', num2str(threshAll, '%.4f')]; ...
        ['{\itN}(ZNCC < \theta) = ', num2str(tempN, '%d')]}, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
    ylabel('Probability Density');
    xlabel('ZNCC ( a.k.a. Pearson''s {\itr} )');
    xTick = linspace(min(corrTrl2Avg), max(corrTrl2Avg), 10);
    set(gca, 'XTick', atanh(xTick), 'XTickLabel', num2str(xTick.', '%.3f'), 'XTickLabelRotation', 45);

    fig4save{2} = figure;
    set(gcf, 'Position', ...
        [10, 50, ...
        560, 315].*scaleFactor);
    xLim = [0, trlsResampNum + 1];
    t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    plot(tempRigid.corrTrl2Avg, 'Color', '#0072BD');
    hold on;
    plot(corrTrl2Avg, 'Color', '#77AC30');
    hold off;

    hold on;
    plot(xlim, repmat(threshAll, 1, 2), 'r');
    hold off;
    title('Chronological Order');
    xlim(xLim);

    nexttile;
    [tempPlot, tempID] = sort(corrTrl2Avg);
    lgd(1) = semilogx(tempRigid.corrTrl2Avg(tempID), 'Color', '#0072BD');
    hold on;
    lgd(2) = semilogx(tempPlot, 'Color', '#77AC30');
    yLim = ylim;
    lgd(3) = semilogx([xLim(1) + 1, repmat(tempN + 0.5, 1, 2)], [repmat(threshAll, 1, 2), yLim(1)], 'r');
    hold off;
    legend(lgd, {'Rigid'; 'nonRigid'; 'threshold'}, 'box', 'off', ...
        'Location', 'southeast');
    title('ZNCC Ascending Order');
    xlim(xLim + [1 0]);
    xLim = xlim;
    yLim = ylim;
    text(10^(0.025*diff(log10(xLim))), yLim(1) + 0.9*diff(yLim), {['\theta = ', num2str(threshAll, '%.4f')]; ...
        ['{\itN}(ZNCC < \theta) = ', num2str(tempN, '%d')]}, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    xlabel(t, 'trial number');
    ylabel(t, 'ZNCC ( a.k.a. Pearson''s {\itr} )');


    fig4save{3} = figure;
    set(gcf, 'Position', ...
        [40 + 2*560 + 315, 50, ...
        400, 315].*scaleFactor);
    imagesc(corrTrl2Trl);
    cb = colorbar;
    title(cb, 'ZNCC', 'FontAngle', 'normal', 'VerticalAlignment', 'baseline');
    colormap(flip(gray));
    axis tight;
    axis equal;

    for n = 1:length(fig4save)
        fig4save{n}.Renderer = 'painters';
    end
    if saveImageON
        saveas(fig4save{1}, fullfile(savePath, 'ZNCC_trlPDF_nonRigid.svg'));
        saveas(fig4save{2}, fullfile(savePath, 'ZNCC_trl_nonRigid.svg'));
        saveas(fig4save{3}, fullfile(savePath, 'ZNCC_trl2trl_nonRigid.svg'));
        imwrite(ORGNLgreen(frmsTrlAvgd), ...
            fullfile(savePath, 'frmsTrlAvgd_ORGNLgreen_uint8_nonRigid.tif'), ...
            'Compression', 'none');
        imwrite(CLAHEslight(frmsTrlAvgd), ...
            fullfile(savePath, 'frmsTrlAvgd_CLAHEslight_uint8_nonRigid.tif'), ...
            'Compression', 'none');
        imwrite(CLAHEstrong(frmsTrlAvgd), ...
            fullfile(savePath, 'frmsTrlAvgd_CLAHEstrong_uint8_nonRigid.tif'), ...
            'Compression', 'none');
    end

    figure;
    set(gca, 'Position', [0, 0, 1, 1]);
    set(gcf, 'Position', ...
        [10, basePixSize - FOVpixDim(2) - 110, ...
        FOVpixDim].*scaleFactor);
    imshow(CLAHEslight(tempRigid.frmsTrlAvgd));

    figure;
    set(gca, 'Position', [0, 0, 1, 1]);
    set(gcf, 'Position', ...
        [20 + FOVpixDim(1), basePixSize - FOVpixDim(2) - 110, ...
        FOVpixDim].*scaleFactor);
    imshow(CLAHEslight(frmsTrlAvgd));

    if playMovieON(2)
        loadNameMOTcrctd_nonRigid = fullfile(savePath, ['frmsMOTcrctd-', ...
            num2str(frmsResampNumByTrl), 'avgd', suffix1, '.', outputFormat]);
        trlsResamp = permute(h5read(loadNameMOTcrctd_nonRigid, ['/', h5GroupName]), ...
            [2, 1, 3]);
        movieTEMP = implay(single(rescale(trlsResamp)), 10);
        movieTEMP.Parent.Position = ...
            [30 + 2*FOVpixDim(1), basePixSize - FOVpixDim(2) - 134, ...
            FOVpixDim + [0, 24]].*scaleFactor;
        clear frms4resamp;
    end

end

%% Step 10: Prepare `ops.npy` for suite2p
if opsPrpr4suite2pON
    loadName = dir(fullfile(savePath, 'frms4ROIdtct', 'plane0', ...
        ['data.', outputFormatFinal]));
    loadPath = loadName.folder;

    loadNameAvgInfo_Final = fullfile(savePath, 'frmsAvgInfo_Final.mat');
    load(loadNameAvgInfo_Final, 'frmsAllAvgd');

    if ne(pyenv().Status, 'Loaded')
        pyenv('Version', pythonPath);
    end

    savePathFull = dir(savePath);
    savePathFull = savePathFull(1).folder;
    
    ops = pyrunfile(['default_ops_', suite2pVersion, '.py'], 'default');
    ops = ops();
    ops_add = py.dict(pyargs( ...
        'nframes', int32(frmsFullNum), ...
        'Ly', int32(FOVpixDim(2)), ...
        'Lx', int32(FOVpixDim(1)), ...
        'yrange', py.numpy.array(int32([0, FOVpixDim(2)])), ...
        'xrange', py.numpy.array(int32([0, FOVpixDim(1)])), ...
        'meanImg', py.numpy.array(min(frmsAllAvgd, peakBRT)), ...
        'meanImgE', py.numpy.array(CLAHEstrong(frmsAllAvgd)), ...
        'data_path', py.list({loadPath}), ...
        'input_format', outputFormatFinal, ...
        'save_path0', savePathFull, ...
        'save_folder', 'frms4ROIdtct', ...
        'tau', 0.7, ... % GCaMP6f/m/s: 0.7/1.0/1.25~1.5
        'fs', FPS, ...
        'save_mat', true, ...
        'combined', false, ...
        'batch_size', int32(FPS*100), ...
        'do_registration', int32(0), ...
        'smooth_sigma_time', 0.5, ...
        'nonrigid', false, ...
        'sparse_mode', true, ... % false: sourcery mode
        'spatial_scale', int32(2), ... % If you are not sure, then just set it to `0`.
        'connected', false, ... % only for sourcery mode
        'smooth_masks', true, ... % only for sourcery mode
        'nbinned', int32(floor(frmsFullNum/FPS)), ...
        'threshold_scaling', 0.25, ... % 0.315 | 0.335 | 0.375 | 0.425 | 0.5
        'active_percentile', 95, ... % *** 0~100, affect a lot of things about ROI detection when > 0
        'splitThres', 1.2, ... % EVafterSplit/EVb4Split
        'max_iterations', 8, ... % *** max ROI number before remove overlaps (actually it was `250 * max_iterations`)
        'max_overlap', 0.75, ...
        'high_pass', int32(8*FPS), ... % <10: gaussian_filter; >=10: rolling_mean_filter
        'spatial_hp_detect', int32(32), ... % uniform filter, recommended value: 2^n
        'denoise', false, ... % open it if data is very noisy
        'block_size', py.list(int32(flip(FOVpixDim)*2)), ... % also affects PCA denoise
        'PCAdenoiseThres', 5, ...% should be larger with larger 'block_size' (affect how many components to keep)
        'anatomical_only', int32(0), ...
        'diameter', int32(0), ... % If you are not sure, then just set it to `0`.
        'neuropil_extract', true, ...
        'circular_neuropil', true, ... % otherwise the neuropil mask is a box
        'inner_neuropil_radius', int32(2), ... % number of pixels to keep between ROI and neuropil donut
        'min_neuropil_pixels', int32(3), ... % thickness of neuropil donut (if 'circular_neuropil' is true)
        'lam_percentile', 50, ...
        'allow_overlap', false, ...
        'baseline', 'maximin', ...
        'win_baseline', 60, ... % win = int(win_baseline*fs)
        'sig_baseline', 10, ...
        'prctile_baseline', 8, ...
        'neucoeff', 0.7));
    update(ops, ops_add);
    py.numpy.save(fullfile(loadPath, 'ops.npy'), ops);

end

%% Shut Down the Parallel Pool
if deletePoolON
    delete(gcp('nocreate'));
end

if moveOriH5ON(2)
    delete(fullfile(loadPathH5, ['rawData-', Date, '-', FOV, '.', outputFormat]));
end


toc(tStart);

