function [B2, I] = autoMatch(A1, A2, theta)
%% Auto Matching
% find the element from A2, which is closest to each element in A1. 

n = length(A1);
dA = A1' - A2;
dAlogic = [theta(1) <= dA(:, 1:end - 1) & dA(:, 1:end - 1) <= theta(2), ...
    abs(dA(:, end)) <= theta(end)];
[I, J] = find(dAlogic);

disp([J, I, dA(dAlogic), A2(I)]);
if length(J) ~= n || any(J' ~= 1:n)
    warning('Need artificial judging!');
end

B2 = A2(I);

end

