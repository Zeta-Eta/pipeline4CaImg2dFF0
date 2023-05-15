function [shifts, Mf] = register_frame_Rigid_4GPU_outputDIY(Yt, fftTemp_conj, nf, data_type, ...
                    us_fac, pixNum, pixDim, pixNum2, pixDim2, Nr, Nc, Nr2, Nc2, padID, ...
                    dftNum, dftNum2, dftDim, dftShift, kernR, kernC, interpGridR, interpGridC, ...
                    nf2, frmsID)

Yt = gpuArray(Yt);
nf = gpuArray(nf);

shifts = dftregistration_DIY_4GPU(fft(fft(Yt, [], 1), [], 2), fftTemp_conj, nf, ...
                    us_fac, pixNum, pixDim, pixNum2, pixDim2, Nr, Nc, Nr2, Nc2, padID, ...
                    dftNum, dftNum2, dftDim, dftShift, kernR, kernC);
if nf2
    Mf = zeros([pixDim, nf2], data_type, 'gpuArray');
    for i = 1:nf2
        Mf(:, :, i) = interp2(interpGridC, interpGridR, Yt(:, :, frmsID(i)), ...
            interpGridC + shifts(frmsID(i), 2), interpGridR + shifts(frmsID(i), 1), 'linear', 0);
    end
    Mf = gather(Mf);
else
    Mf = [];
end

shifts = gather(shifts);

end