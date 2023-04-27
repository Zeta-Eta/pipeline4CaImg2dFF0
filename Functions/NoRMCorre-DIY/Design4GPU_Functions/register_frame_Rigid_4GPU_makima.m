function [shifts, Mf] = register_frame_Rigid_4GPU_makima(Yt, fftTemp_conj, nf, data_type, ...
                    us_fac, pixNum, pixDim, pixNum2, pixDim2, Nr, Nc, Nr2, Nc2, padID, ...
                    dftNum, dftNum2, dftDim, dftShift, kernR, kernC, interpGridR, interpGridC)

Yt = gpuArray(Yt);
nf = gpuArray(nf);

shifts = dftregistration_DIY_4GPU(fft(fft(Yt, [], 1), [], 2), fftTemp_conj, nf, ...
                    us_fac, pixNum, pixDim, pixNum2, pixDim2, Nr, Nc, Nr2, Nc2, padID, ...
                    dftNum, dftNum2, dftDim, dftShift, kernR, kernC);

shifts = gather(shifts);
Yt = gather(Yt);
interpGridR = gather(interpGridR);
interpGridC = gather(interpGridC);
nf = gather(nf);

Mf = zeros([pixDim, nf], data_type);
for i = 1:nf
%     tempYt = Yt(:, :, i);
%     F = griddedInterpolant(Yt(:, :, i), 'makima', 'none');
%     tempMf = F({interpGridR + shifts(i, 1), interpGridC + shifts(i, 2)});
%     tempMf(isnan(tempMf)) = 0;
    Mf(:, :, i) = min(1, max(0, interp2(interpGridC, interpGridR, Yt(:, :, i), interpGridC + shifts(i, 2), interpGridR + shifts(i, 1), 'makima', 0)));

%     tempMf = imwarp(tempYt, cat(3, repmat(shifts(i, 2), pixDim), repmat(shifts(i, 1), pixDim)), 'linear', 'FillValues', add_value);

%     minY = min(tempYt(:));
%     maxY = max(tempYt(:));
%     tempMf(tempMf < minY) = minY;
%     tempMf(tempMf > maxY) = maxY;

%     Mf(:, :, i) = tempMf;
end

end