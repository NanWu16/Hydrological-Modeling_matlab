function [snew, fnew] = cceua_N(fn, s, sf, bl, bu, x_obs, y_obs, UL, CON)

    [nps,nopt]=size(s);
    n = nps;
    % m = nopt;
    alpha = 1.0;
    beta = 0.5;

    % Assign the best and worst points:
    % sb=s(1,:); fb=sf(1);
    sw=s(n,:); fw=sf(n);

    % Compute the centroid of the simplex excluding the worst point:
    ce=mean(s(1:n-1,:));

    % Attempt a reflection point
    snew = ce + alpha*(ce-sw);

    % Check if is outside the bounds:
    ibound=0;
    s1=snew-bl; idx=find(s1<0, 1); if ~isempty(idx); ibound=1; end
    s1=bu-snew; idx=find(s1<0, 1); if ~isempty(idx); ibound=2; end

    if ibound >=1
        snew = bl + rand(1,nopt).*(bu-bl);
    end
    fnew = fn(x_obs, snew, y_obs, UL, CON);

    % Reflection failed; now attempt a contraction point:
    if fnew > fw
        snew = sw + beta*(ce-sw);
        fnew = fn(x_obs, snew, y_obs, UL, CON);

        % Both reflection and contraction have failed, attempt a random point;
        if fnew > fw
            snew = bl + rand(1,nopt).*(bu-bl);
            fnew = fn(x_obs, snew, y_obs, UL, CON);
        end
    end
end