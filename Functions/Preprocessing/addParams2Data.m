function data = addParams2Data(data, kappa, rotA, sqnsN)
%add params to data

tempI = floor((data.sqnsAll - 1)./sqnsN) + 1;

kappa = [1; kappa];
data.kappa = kappa(tempI);

rotA = [NaN; rotA];
data.rotA = rotA(tempI);

end

