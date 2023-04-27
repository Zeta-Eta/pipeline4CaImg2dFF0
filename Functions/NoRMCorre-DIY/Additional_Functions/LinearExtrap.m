function [diffArange, x] = LinearExtrap(XYin, XYout, cntrPos)
% Linear Extrapolation
XYin = XYin - cntrPos;
XYout = XYout - cntrPos;

[Ain, ~] = cart2pol(XYin(:, 1), XYin(:, 2)); % theta: rot Angle; rho: mov Distance
[Aout, ~] = cart2pol(XYout(:, 1), XYout(:, 2));

diffA = mod(Aout - Ain.', 2*pi);
[~, diffArange(:, 1)] = min(diffA, [], 2);
[~, diffArange(:, 2)] = max(diffA, [], 2);
A = reshape(XYin(diffArange', :).', 2, 2, []);
b = reshape(XYout.', 2, 1, []);
x = squeeze(pagemldivide(A, b)).';

end

