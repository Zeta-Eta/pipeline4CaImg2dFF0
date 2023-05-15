function [shifts, CCmax] = dftregistration_DIY_4GPU(buf1ft, buf2ft_conj, nf, us_fac, ...
    pixNum, pixDim, pixNum2, pixDim2, Nr, Nc, Nr2, Nc2, padID, ...
    dftNum, dftNum2, dftDim, dftShift, kernR, kernC)
%% !!! `buf2ft_conj` is the conj FFT of the template

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

if us_fac == 0
    % without registration
    row_shift = zeros(nf, 1);
    col_shift = row_shift;
    shifts = [row_shift, col_shift];
elseif us_fac == 1
    % Single pixel registration
    CC = ifft(ifft(buf1ft.*buf2ft_conj, [], 2), [], 1, 'symmetric'); % <- cross-correlation
    % CC = ifft(ifft(sign(buf1ft.*buf2ft_conj), [], 2), [], 1, 'symmetric'); % <- phase correlation
    [shifts, CC] = max(reshape(CC, pixNum, nf), [], 1); % row vector
    [row_shift, col_shift] = ind2sub(pixDim, CC); % row vector
    % Now change shifts so that they represent relative shifts and not indices
    row_shift = Nr(row_shift);
    col_shift = Nc(col_shift);
    shifts = [row_shift, col_shift, shifts.'];
elseif us_fac > 1
    % Start with us_fac == 2
    CC = buf1ft.*buf2ft_conj;
    CC2 = zeros([pixDim2, nf], 'single', 'gpuArray');
    CC2(repmat(padID, 1, 1, nf)) = buf1ft.*buf2ft_conj;
    CC2 = ifft(ifft(CC2, [], 2), [], 1, 'symmetric'); % <- cross-correlation
    % CC2 = ifft(ifft(sign(CC2), [], 2), [], 1, 'symmetric'); % <- phase correlation
    [shifts, CC2] = max(reshape(CC2, pixNum2, nf), [], 1); % row vector
    [row_shift, col_shift] = ind2sub(pixDim2, CC2); % row vector
    % Now change shifts so that they represent relative shifts and not indices
    row_shift = Nr2(row_shift); % col vector
    col_shift = Nc2(col_shift); % col vector
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if us_fac > 2
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*us_fac);
        col_shift = round(col_shift*us_fac);
        % Matrix multiply DFT around the current shift estimate
        CC = abs(pagemtimes( ...
            pagemtimes( ...
            exp(pagemtimes(((0:dftNum-1).' - dftShift + reshape(row_shift, 1, 1, nf)), kernR)), ...
            CC), ...
            exp(pagemtimes(kernC, ((0:dftNum-1) - dftShift + reshape(col_shift, 1, 1, nf))) ...
            ) ...
            ));
        % Locate maximum and map back to original pixel grid
        CC = reshape(CC, dftNum2, nf);
        shifts = [shifts; mean(CC, 1); median(CC, 1)];
        [CCmax, CC] = max(CC, [], 1); % row vector
        shifts = [shifts; CCmax];
        [rloc, cloc] = ind2sub(dftDim, CC.'); % col vector
        row_shift = (row_shift + rloc - dftShift - 1)/us_fac;% dftShift: Center of output array at dftshift+1
        col_shift = (col_shift + cloc - dftShift - 1)/us_fac;
    end
    shifts = [row_shift, col_shift, shifts.'];
end

end