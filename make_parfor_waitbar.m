function [h, q] = make_parfor_waitbar(N, titleStr)
% Waitbar that updates from parfor workers via a DataQueue
h = waitbar(0, sprintf('%s 0%%', titleStr));
q = parallel.pool.DataQueue;
cnt = 0;  % lives on client
afterEach(q, @update);
    function update(~)
        cnt = cnt + 1;
        frac = min(cnt / N, 1);
        if isvalid(h)
            waitbar(frac, h, sprintf('%s %2.0f%%', titleStr, 100*frac));
        end
    end
end