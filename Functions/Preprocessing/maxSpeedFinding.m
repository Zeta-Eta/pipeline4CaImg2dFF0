function tID_M = maxSpeedFinding(Response, eventTime, eyeMov4mula, layoutID, layoutInf, dThres, vThres, respFixT, odrN)
% max(v) finding
z0 = layoutInf.r(layoutID, :).*exp(1i.*layoutInf.theta(layoutID, :));% Layout
% tNode = eventTime{1, [6, 7 + (1:odrN)]};
tNode = eventTime{1, [end-odrN-2, end-odrN:end-1]};
T = tNode(1):tNode(end);
X = ppval(eyeMov4mula.x, T);
Y = ppval(eyeMov4mula.y, T);
Z = X + Y.*1i;
V_X = ppval(eyeMov4mula.v_x, T);
V_Y = ppval(eyeMov4mula.v_y, T);
V = vecnorm([V_X; V_Y], 2, 1);
[~, ~, vThres2, ~] = isoutlier(V);
vThres = vThres*vThres2; % vThres = 2;

tID_M = NaN(1, odrN);
for odr = 1:odrN
    t = (tNode(odr):tNode(odr + 1)) - T(1) + 1;
    z = Z(t);
    d = abs(z - z0(Response(odr)));
    v = V(t);
    
    % time index when d < dThres
    tIDd = [find(diff(d < dThres) == -1), length(t) - respFixT];
    tIDd(tIDd > tIDd(end)) = [];
    tIDd0 = find(diff(d < dThres) == 1);
    tIDd0(tIDd0 > tIDd(end)) = [];
    if length(tIDd0) < length(tIDd)
        tIDd0 = [1, tIDd0];
    end
    tIDd = tIDd(find(tIDd - tIDd0 >= 15, 1, 'first'));
    
    % time index when v > vThres
    vThresUSE = vThres;
    while 1
        diff_tIDv = diff(v(1:tIDd) > vThresUSE);
        tIDv = [find(diff_tIDv == 1, 1, 'last'), ...
            find(diff_tIDv == -1, 1, 'last')];
        if isempty(tIDv)
            if any(v(1:tIDd) > vThresUSE)
                v(1) = 0; % set v(1) to 0
            else
                vThresUSE = 0.9*vThresUSE;
            end
        elseif length(tIDv) == 1
            temp = find(diff(v > vThresUSE) == -1, 1);
            if isempty(temp)
                temp = length(t);
            end
            tIDv = [1, temp];
            break;
        elseif length(tIDv) == 2
            if tIDv(1) > tIDv(2)
                temp = find(diff_tIDv(1:tIDv(2)) == 1, 1, 'last');
                if isempty(temp)
                    temp = 1;
                end
                tIDv(1) = temp;
            end
            timeThres = 15;
            while 1
                if tIDv(1) == 1
                    break;
                elseif tIDv(2) - tIDv(1) < timeThres
                    tIDv2 = [find(diff_tIDv(1:tIDv(1) - 1) == 1, 1, 'last'), ...
                        find(diff_tIDv(1:tIDv(1) - 1) == -1, 1, 'last')];
                    if isempty(tIDv2)
                        tIDv2 = [1, tIDv(1)];
                    elseif length(tIDv2) == 1
                        tIDv2 = [1, tIDv2(1)];
                    end
                    timeThres = timeThres + tIDv(2) - tIDv2(2);
                    tIDv(1) = tIDv2(1);
                else
                    break;
                end
            end
            break;
        end
    end
    k = 0;
    [~, tID_M(odr)] = max(v((tIDv(1):tIDv(2)) + k));
    
    tID_M(odr) = tID_M(odr) + tIDv(1) + k - 1;
    
    
end

