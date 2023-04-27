function [DIYparams, shifts, templateUpdating, options] = ...
    normcorre_batch_Rigid_4GPU_outputDIY(Y, options, template, SPMDcoreNum, DIYparamsON, clearGPU, iterNum, saveON, plotON, avgON, frmsID)
% ONLY FOR `options.output_type == 'h5'`
% `SPMDcoreNum` will be actually `SPMDcoreNum + 1`

% Based on the dftregistration.m function from Manuel Guizar and Jim Fienup

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction (optional, rigid registration is performed if not provided)
% template:         provide template

% OUTPUTS
% DIYparams:        DIY parameters
% shifts:           calculated shifts
% templateUpdating: calculated template

if plotON
    screenPixSize = get(0, 'ScreenSize');
    basePixSize = 1080; % never change it
    scaleFactor = screenPixSize(4)/basePixSize;

    CLAHEstrong = @(I) adapthisteq(rescale(I), 'NumTiles', repmat(2^5, 1, 2), ...
        'Distribution', 'rayleigh', 'Alpha', 0.32, 'NBins', 2^16, 'Range', 'original');
end

%% first determine filetype

if isa(Y,'char')
    [~,~,ext] = fileparts(Y);
    ext = ext(2:end);
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff')
        tiffInfo = imfinfo(Y);
        filetype = 'tif';
        T = length(tiffInfo);
        sizY = [tiffInfo(1).Height,tiffInfo(1).Width,T];
    elseif strcmpi(ext,'mat')
        filetype = 'mem';
        Y = matfile(Y,'Writable',true);
        details = whos(Y);
        var_sizes = [details.bytes];
        [~,var_ind] = max(var_sizes);
        var_name = details(var_ind).name;
        sizY = size(Y,var_name);
        T = sizY(end);
    elseif strcmpi(ext,'hdf5') || strcmpi(ext,'h5')
        filetype = 'hdf5';
        fileinfo = hdf5info(Y);
        data_name = fileinfo.GroupHierarchy.Datasets.Name;
        sizY = fileinfo.GroupHierarchy.Datasets.Dims;
        T = sizY(end);
    elseif strcmpi(ext,'raw')
        filetype = 'raw';
        fid = fopen(Y);
        FOV = [options.d1,options.d2];
        bitsize = options.bitsize;
        imsize = FOV(1)*FOV(2)*bitsize;                                                   % Bit size of single frame
        current_seek = ftell(fid);
        fseek(fid, 0, 1);
        file_length = ftell(fid);
        fseek(fid, current_seek, -1);
        T = file_length/imsize;
        sizY = [FOV,T];
        fclose(fid);        
    elseif strcmpi(ext,'avi')
        filetype = 'avi';
        sizY = size(read_file(Y));
        FOV = sizY(1:2);
        T = sizY(end);
    end    
elseif isobject(Y)
    filetype = 'mem';
    var_name = 'Y';
    sizY = size(Y,var_name);
    T = sizY(end);
else % array loaded in memory
    filetype = 'mat';
    sizY = size(Y);
    T = sizY(end);
end

nd = length(sizY)-1; % determine whether imaging is 2d or 3d
sizY = sizY(1:nd);
otherdims = repmat({':'}, 1, nd);
%% set default parameters if not present

if ~exist('options','var') || isempty(options)
    options = NoRMCorreSetParms('d1', sizY(1), 'd2', sizY(2));
    if nd > 2; options.d3 = sizY(3); end
end

bin_width = options.bin_width;
init_batch = options.init_batch;
print_msg = options.print_msg;


%% read initial batch and compute template

init_batch = min(T,init_batch);
interval = ceil((T - init_batch + 1)*0.5):floor((T + init_batch)*0.5);
switch filetype
    case 'tif'
        Y_temp = read_file(Y,interval(1),init_batch,[],tiffInfo);
    case 'hdf5'
        Y_temp = read_file(Y,interval(1),init_batch);        
    case 'avi'
        Y_temp = read_file(Y,interval(1),init_batch);
    case 'mem'
        Y_temp = Y.(var_name)(otherdims{:},interval);
    case 'mat'
        Y_temp = Y(otherdims{:},interval);
    case 'raw'
        Y_temp = read_raw_file(Y,interval(1),init_batch,FOV,bitsize);
end
data_type = class(Y_temp);
% Y_temp = single(Y_temp);
% use_proj = true;

[d1,d2,d3,~] = size(Y_temp);
clear Y_temp;
if nd == 2; d3 = 1; end
%% setup grids for patches
dim = [d1, d2, d3];
T2 = length(frmsID);

if mod(options.mem_batch_size, bin_width)
    warning('error: mod(options.mem_batch_size, bin_width) ~= 0');
    if isempty(frmsID)
        return;
    end
end

mem_buffer = squeeze(zeros([dim, options.mem_batch_size], data_type));

if saveON
    if exist(options.h5_filename,'file')
        [pathstr,fname,ext] = fileparts(options.h5_filename);
        new_filename = fullfile(pathstr,[fname,'_',datestr(now,30),ext]);
        warning_msg = ['File ',options.h5_filename,'already exists. Saving motion corrected file as',new_filename];
        warning('%s',warning_msg);
        options.h5_filename = new_filename;
    end
    if avgON
        h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,T2./options.mem_batch_size],'Datatype',data_type);
    else
        h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,options.mem_batch_size],'Datatype',data_type);
    end
end


if print_msg; fprintf('Template initialization complete. Now registering all the frames with new template. \n'); end
%% Set Options
options.pixDim = [options.d1, options.d2];
options.pixNum = options.d1.*options.d2;
options.pixDim2 = 2.*options.pixDim;
options.pixNum2 = 4.*options.pixNum;

options.padID = zeros(options.pixDim2, 'logical');
options.padID([1:0.25*end, 0.75*end+1:end], [1:0.25*end, 0.75*end+1:end]) = 1;

options.kernR = (1i*2*pi/(options.d1*options.us_fac))*(ifftshift(0:options.d1-1)   - floor(options.d1/2));
options.kernC = (1i*2*pi/(options.d2*options.us_fac))*(ifftshift(0:options.d2-1).' - floor(options.d2/2));
options.interpGridR = (1:options.d1).';
options.interpGridC = 1:options.d2;

options.dftNum = ceil(options.us_fac*1.5);
options.dftDim = repmat(options.dftNum, 1, 2);
options.dftNum2 = options.dftNum.^2;
options.dftShift = fix(options.dftNum/2);

options.Nr = ifftshift(-fix(options.d1/2):ceil(options.d1/2)-1).';
options.Nc = ifftshift(-fix(options.d2/2):ceil(options.d2/2)-1).';
options.Nr2 = ifftshift(-fix(options.d1):ceil(options.d1)-1).'./2;
options.Nc2 = ifftshift(-fix(options.d2):ceil(options.d2)-1).'./2;


us_fac = gpuArray(options.us_fac);
pixNum = gpuArray(options.pixNum);
pixDim = gpuArray(options.pixDim);
pixNum2 = gpuArray(options.pixNum2);
pixDim2 = gpuArray(options.pixDim2);
Nr = gpuArray(options.Nr);
Nc = gpuArray(options.Nc);
Nr2 = gpuArray(options.Nr2);
Nc2 = gpuArray(options.Nc2);
padID = gpuArray(options.padID);
dftNum = gpuArray(options.dftNum);
dftNum2 = gpuArray(options.dftNum2);
dftDim = gpuArray(options.dftDim);
dftShift = gpuArray(options.dftShift);
kernR = gpuArray(options.kernR);
kernC = gpuArray(options.kernC);
interpGridR = gpuArray(options.interpGridR);
interpGridC = gpuArray(options.interpGridC);

fftTemp_conj = gpuArray(conj(fft(fft(single(template), [], 1), [], 2)));

% varList = {'add_value', 'us_fac', 'pixNum', 'pixDim', 'pixNum2', 'pixDim2', ...
%     'Nr', 'Nc', 'Nr2', 'Nc2', 'padID', 'dftNum', 'dftNum2', 'dftDim', 'dftShift', ...
%     'kernR', 'kernC'};
% for i = 1:length(varList)
%     eval([varList{i}, ' = gpuArray(options.', varList{i}, ');']);
% end

%%

bin_num = ceil(T/bin_width);
bin_num2 = ceil(bin_num/SPMDcoreNum);
tID = NaN(SPMDcoreNum, bin_num2);
tID(1:bin_num) = 1;
tID = tID.';
tIDnum = sum(tID, 1, 'omitnan');

tID2 = permute(repmat(tID, 1, 1, bin_width), [3, 1, 2]);
tempID = find(~isnan(tID2));
tID2(tempID(T+1:end)) = NaN;
IfSave = NaN(size(tID2));
IfSave(tID2 == 1) = ismember(1:T, frmsID);
IfSave = permute(IfSave, [2 3 1]);
IfSaveDiff = sum(IfSave, 3, 'omitnan');
IfSaveDiff(isnan(tID)) = NaN;
IfSaveDiffCumsum = cumsum(IfSaveDiff, 1, 'omitnan');
IfSaveDiffCumsum(isnan(tID)) = NaN;
[~, IfSaveID]= max(IfSaveDiffCumsum, [], 1, 'omitnan');

binIDstt = ceil(frmsID([0, find(diff(frmsID)>1)] + 1)/bin_width);
binIDend = ceil(frmsID([find(diff(frmsID)>1), end])/bin_width);
addID = cumsum([0, tIDnum(1:end-1)]);
binIDsave = IfSaveID + addID;

tempID = binIDstt - binIDsave.';
tempID(tempID > 0) = NaN;
[~, tempID] = max(tempID, [], 2);
edgeIDstt = binIDstt(tempID) - addID;

tempID = binIDend - binIDsave.';
tempID(tempID < 0) = NaN;
[~, tempID] = min(tempID, [], 2);
edgeIDend = binIDend(tempID) - addID;

coreIDmove = find((edgeIDend - tIDnum) > 0);
edgeIDlen = tIDnum(coreIDmove) - edgeIDstt(coreIDmove) + 1;
tIDnum(coreIDmove) = tIDnum(coreIDmove) - edgeIDlen;
tIDnum(coreIDmove + 1) = tIDnum(coreIDmove + 1) + edgeIDlen;

tID = NaN(max(tIDnum), SPMDcoreNum);
for c = 1:SPMDcoreNum
    tID(1:tIDnum(c), c) = 1;
end

tID2 = permute(repmat(tID, 1, 1, bin_width), [3, 1, 2]);
tempID = find(~isnan(tID2));
tID2(tempID(T+1:end)) = NaN;
IfSave = NaN(size(tID2));
IfSave(tID2 == 1) = ismember(1:T, frmsID);
IfSave = permute(IfSave, [2 3 1]);
IfSaveDiff = sum(IfSave, 3, 'omitnan');
IfSaveDiff(isnan(tID)) = NaN;
IfSaveDiffCumsum = cumsum(IfSaveDiff, 1, 'omitnan');
IfSaveDiffCumsum(isnan(tID)) = NaN;
[IfSaveDiffSum, IfSaveID] = max(IfSaveDiffCumsum, [], 1, 'omitnan');
IfSaveSTT = cumsum([1, IfSaveDiffSum(1:end-1)]);

if avgON
    IfSaveDiffCumsum2 = IfSaveDiffCumsum./options.mem_batch_size;
    IfSaveDiffSum = max(IfSaveDiffCumsum2, [], 1, 'omitnan');
    IfSaveSTT = cumsum([1, IfSaveDiffSum(1:end-1)]);
end


tIDstt = tID;
tIDend = tID;
tIDstt(tID == 1) = 1:bin_width:T;
tIDend(tID == 1) = unique([bin_width:bin_width:T, T], 'stable');

% tIDnum = sum(tID, 1, 'omitnan');
% sCRnum = sum(tID, 2, 'omitnan');

if saveON
    tempH5name = cell(SPMDcoreNum, 1);
    for c = 1:SPMDcoreNum
        tempH5name{c} = [options.h5_filename(1:end-3), '_temp', num2str(c, '%03d'), '.h5'];
        if avgON
            h5create(tempH5name{c},['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,IfSaveDiffSum(c)],'Datatype',data_type);
        else
            h5create(tempH5name{c},['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,options.mem_batch_size],'Datatype',data_type);
        end
    end
end

tic
for iter = 1:iterNum
    if print_msg; fprintf(['>>>Start iteration ', num2str(iter), ' out of ', num2str(iterNum), '... \n']); end
    if iter > 1
        fftTemp_conj = gpuArray(conj(fft(fft(single(templateUpdating), [], 1), [], 2)));
    end
    spmd(SPMDcoreNum + 1)
        switch labindex
            case SPMDcoreNum + 1
                tic;
                tNOW = 0;
                Ttxt = num2str(T, '%d');
                for m = 1:max(tIDnum)
                    for n = find(tID(m, :) == 1)
                        tNOW = labReceive(n, m) + tNOW;
                    end
                    tempTIME = toc;
                    if ~mod(m, print_msg) || m == max(tIDnum)
                        tNOWtxt = pad(num2str(tNOW, '%d'), length(Ttxt), 'left');
                        tempTIMEtxt = pad(num2str(tempTIME, '%.3f'), 8, 'left');
                        fprintf([tNOWtxt ' out of ' Ttxt ' frames registered... ' tempTIMEtxt ' seconds passed...\n']);
                    end
                end
            otherwise
                if iter == iterNum
                    shifts = cell(tIDnum(labindex), 1);
                end
                templateUpdating = zeros(sizY(1:nd));
                for j = 1:tIDnum(labindex)
                    tSTT = tIDstt(j, labindex);
                    tEND = tIDend(j, labindex);
                    if labindex == SPMDcoreNum && j == tIDnum(SPMDcoreNum)
                        nf = tEND - tSTT + 1;
                    else
                        nf = bin_width;
                    end
                    nf2 = IfSaveDiff(j, labindex);
                    switch filetype
                        case 'tif'
                            Ytm = single(read_file(Y, tSTT, nf, [], tiffInfo));
                        case 'avi'
                            Ytm = single(read_file(Y, tSTT, nf));
                        case 'hdf5'
                            Ytm = single(h5read(Y,data_name,[ones(1,nd),tSTT],[sizY(1:nd),nf]));
                        case 'mem'
                            Ytm = single(Y.(var_name)(otherdims{:},tSTT:tEND));
                        case 'mat'
                            Ytm = single(Y(otherdims{:},tSTT:tEND));
                        case 'raw'
                            Ytm = single(read_raw_file(Y,tSTT,nf,FOV,bitsize));
                    end

                    if iter == iterNum
                        [shifts{j}, Mf] = register_frame_Rigid_4GPU_outputDIY(Ytm, fftTemp_conj, nf, data_type, ...
                            us_fac, pixNum, pixDim, pixNum2, pixDim2, Nr, Nc, Nr2, Nc2, padID, ...
                            dftNum, dftNum2, dftDim, dftShift, kernR, kernC, interpGridR, interpGridC, ...
                            nf2, find(IfSave(j, labindex, :)));
                    else
                        [~, Mf] = register_frame_Rigid_4GPU_outputDIY(Ytm, fftTemp_conj, nf, data_type, ...
                            us_fac, pixNum, pixDim, pixNum2, pixDim2, Nr, Nc, Nr2, Nc2, padID, ...
                            dftNum, dftNum2, dftDim, dftShift, kernR, kernC, interpGridR, interpGridC, ...
                            nf2, find(IfSave(j, labindex, :)));
                    end

                    if nf2
                        templateUpdating = templateUpdating + sum(Mf, nd + 1);
                    end
                    if saveON && iter == iterNum && nf2
                        rem_mem = rem(IfSaveDiffCumsum(j, labindex), options.mem_batch_size);
                        if rem_mem == 0; rem_mem = options.mem_batch_size; end
                        mem_buffer(otherdims{:}, rem_mem-nf2+1:rem_mem) = Mf;
                        if rem_mem == options.mem_batch_size || j == IfSaveID(labindex)
                            if avgON
                                h5write(tempH5name{labindex}, ['/', options.h5_groupname], ...
                                    cast(sum(mem_buffer, 3)./options.mem_batch_size, data_type), ...
                                    [ones(1, nd), IfSaveDiffCumsum2(j, labindex)], ...
                                    [sizY(1:nd), 1]);
                            else
                                h5write(tempH5name{labindex}, ['/', options.h5_groupname], ...
                                    mem_buffer(otherdims{:},1:rem_mem), ...
                                    [ones(1, nd), IfSaveDiffCumsum(j, labindex) - rem_mem + 1], ...
                                    [sizY(1:nd), rem_mem]);
                            end
                        end
                    end

                    labSend(nf, SPMDcoreNum + 1, j);
                end
        end
    end
    templateUpdating = sum(cat(nd + 1, templateUpdating{1:SPMDcoreNum}), nd + 1)./T2;
    
    if plotON
        figure;
        set(gca, 'Position', [0, 0, 1, 1]);
        set(gcf, 'Position', ...
            [(10 + options.pixDim(1))*2 + (10 + options.pixDim(1)/iterNum)*(iter - 1), basePixSize - options.pixDim(2) - 110, ...
            options.pixDim].*scaleFactor);
        imshow(CLAHEstrong(templateUpdating.'));
    end
end
toc


shifts = cell2mat(cat(1, shifts{1:SPMDcoreNum}));

if clearGPU
    delete(gcp('nocreate'));
else
    maxNumCompThreads('automatic');
end

if saveON
    tic
    if avgON
        ZNCC = NaN(T2./options.mem_batch_size, 1);
    else
        ZNCC = NaN(T2, 1);
    end
    for c = 1:SPMDcoreNum
        [N, tempID, tempID2useNum] = autoGroup(IfSaveDiffSum(c), ceil(IfSaveDiffSum(c)/saveON));
        for n = 1:N
            tempH5 = h5read(tempH5name{c}, ['/', options.h5_groupname], ...
                [ones(1, nd), tempID(1, n)], [sizY(1:nd), tempID2useNum(n)]);
            h5write(options.h5_filename, ['/', options.h5_groupname], tempH5, ...
                [ones(1, nd), IfSaveSTT(c) + tempID(1, n) - 1], [sizY(1:nd), tempID2useNum(n)]);
            ZNCC(IfSaveSTT(c) + tempID(1, n) - 2 + (1:tempID2useNum(n))) = corr(double(reshape(tempH5, options.pixNum, [])), templateUpdating(:));
        end
        delete(tempH5name{c});
    end
    clear tempH5;
    toc
end

if DIYparamsON
    DIYparams.shifts = shifts;
    DIYparams.templateInput = template;
    DIYparams.templateOutput = templateUpdating;
    if exist('ZNCC', 'var')
        DIYparams.ZNCC = ZNCC;
    end
else
    DIYparams = [];
end

if print_msg; fprintf('\n'); end


