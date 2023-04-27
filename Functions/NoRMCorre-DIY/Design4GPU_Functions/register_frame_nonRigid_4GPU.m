function Mf = register_frame_nonRigid_4GPU(Yt, nf, RigidDF, nonRigidDF, data_type, pixDim, XgridVec, YgridVec, interpGridR, interpGridC)

Yt = gpuArray(Yt);
nf = gpuArray(nf);
RigidDF = gpuArray(RigidDF);
nonRigidDF = gpuArray(nonRigidDF);

Mf = zeros([pixDim, nf], data_type, 'gpuArray');
for i = 1:nf
    Mf(:, :, i) = interp2(interpGridC, interpGridR, Yt(:, :, i), ...
        interp2(YgridVec, XgridVec, nonRigidDF(:, :, i, 2), interpGridC, interpGridR, 'cubic', 0) + RigidDF(i, 2) + interpGridC, ...
        interp2(YgridVec, XgridVec, nonRigidDF(:, :, i, 1), interpGridC, interpGridR, 'cubic', 0) + RigidDF(i, 1) + interpGridR, ...
        'linear', 0);
end

Mf = gather(Mf);
end