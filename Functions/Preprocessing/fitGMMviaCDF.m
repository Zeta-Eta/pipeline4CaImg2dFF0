function [modelParams, Rsqrd, RNG] = fitGMMviaCDF(data, k, initMethod)
%% fit Gaussian Mixture Model via CDF
dataNum = length(data);

Xdata = reshape(sort(data), 1, dataNum);
Ydata = ((1:dataNum) - 0.5)./dataNum;

dataMin = min(data);
dataMax = max(data);
RNG = [dataMin, dataMax];
MAXabsRNG = max(abs(RNG));
dRNG = diff(RNG, 1, 2);

if strcmp(initMethod, 'k-medoids')
    clstrID = kmedoids(data, k);
else
    clstrID = kmeans(data, k);
end
for i = 1:k
    tempID = clstrID == i;
    tempData = data(tempID);
    if i ~= k
        init.w(i) = mean(tempID);
    end
    init.mu(i) = mean(tempData);
    init.sigma(i) = std(tempData);
end
initParams = [log(init.w), atanh(init.mu./MAXabsRNG), log(init.sigma)]; %

% options = optimset('fminsearch');
% options.TolFun = 1e-7;
% options.TolX = 1e-7;
% options.MaxFunEvals = 1e5;
% options.MaxIter = 1e5;
% options.Display = 'notify';

options = optimoptions('fminunc');
options.TolFun = 1e-7;
options.TolX = 1e-7;
options.MaxFunEvals = 1e5;
options.MaxIter = 1e5;
options.Display = 'notify-detailed';

tic;

lossFunc = @(params) LF(params, k, MAXabsRNG, dRNG, Xdata, Ydata);

% fitdParams = fminsearch(lossFunc, initParams, options);
fitdParams = fminunc(lossFunc, initParams, options);

[SSres, ~, w, mu, sigma] = LF(fitdParams, k, MAXabsRNG, dRNG, Xdata, Ydata);
toc;
modelParams = table(w, mu, sigma);
SSres = SSres./(dataNum - length(initParams));
SStot = sum((Ydata - mean(Ydata, 'all', 'omitnan')).^2, 'all', 'omitnan')./(dataNum - 1);
Rsqrd = 1 - SSres./SStot;

end

function [L, Ymodel, w, mu, sigma] = LF(params, k, MAXabsRNG, dRNG, Xdata, Ydata)

w = exp(-abs(params(1:k - 1)));
w = [w, 1 - sum(w)];
if w(k) < 0
    w = abs(w);
    w = w./sum(w);
end
mu = MAXabsRNG.*tanh(params(k:2*k - 1));
sigma = dRNG.*exp(-abs(params(2*k:end)));

CDF = erfc((mu.' - Xdata)./sigma.'./sqrt(2));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
Ymodel = (0.5.*w) * CDF;

L = sum((Ymodel - Ydata).^2, 'all', 'omitnan');

end
