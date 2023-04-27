function A = read_bin_file_DIY(filename, framenum, window, FOV, bitsize, dataType, dataTypeTrans)

%% Starting from frame number: 'framenum', read 'window'
fid = fopen(filename);
imsize = FOV(1)*FOV(2)*bitsize;                                                   % Bit size of single frame

current_seek = ftell(fid);
fseek(fid, 0, 1);
file_length = ftell(fid);
fseek(fid, current_seek, -1);
frame_length = file_length/imsize;

%%
window = min(window, frame_length);
A = zeros(FOV(1), FOV(2), window, dataType);

for w = 0:window-1                                                          % For each frame inside the window
    fseek(fid, (framenum-1+w)*imsize, 'bof');                                 % Position the file-indicator at the beginning of the frame
    A(:, :, w+1) = fread(fid, [FOV(1), FOV(2)], dataTypeTrans, 0, 'n');                              % Read the frame                                                                  
end
fclose(fid);