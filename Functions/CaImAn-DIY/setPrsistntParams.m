function prsistntParams = setPrsistntParams(n, len)
% set persistent params
prsistntParams.n = 2^ceil(log2(n));
prsistntParams.I = (1:prsistntParams.n-1).'.^2;
prsistntParams.len = len; % default: 7
prsistntParams.L = 2 .* (pi.^2 .* prsistntParams.I).^(1:prsistntParams.len);
prsistntParams.R = -prsistntParams.I.*pi^2;
prsistntParams.tBase = (1 + 0.5.^((1:prsistntParams.len) + 0.5)).*cumprod(1:2:2*prsistntParams.len - 1)./(3*n*sqrt(0.5*pi));
prsistntParams.tPower = 2./(3 + 2.*(1:prsistntParams.len));
prsistntParams.tBaseMin = 2*n*sqrt(pi);

end

