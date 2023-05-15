function [shifts, Mf] = register_frame_Rigid_4GPU(Yt, fftTemp_conj, nf, data_type, ...
                    us_fac, pixNum, pixDim, pixNum2, pixDim2, Nr, Nc, Nr2, Nc2, padID, ...
                    dftNum, dftNum2, dftDim, dftShift, kernR, kernC, interpGridR, interpGridC)

Yt = gpuArray(Yt);
nf = gpuArray(nf);

shifts = dftregistration_DIY_4GPU(fft(fft(Yt, [], 1), [], 2), fftTemp_conj, nf, ...
                    us_fac, pixNum, pixDim, pixNum2, pixDim2, Nr, Nc, Nr2, Nc2, padID, ...
                    dftNum, dftNum2, dftDim, dftShift, kernR, kernC);

Mf = zeros([pixDim, nf], data_type, 'gpuArray');
for i = 1:nf
    Mf(:, :, i) = interp2(interpGridC, interpGridR, Yt(:, :, i), interpGridC + shifts(i, 2), interpGridR + shifts(i, 1), 'linear', 0);
end

shifts = gather(shifts);
Mf = gather(Mf);
end