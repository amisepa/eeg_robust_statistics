function [tfce_score, thresholded_maps] = limo_tfce(varargin)
% Threshold-free cluster enhancement (TFCE)
% tfce = sum( extent(h)^E * h^H * dh )
% 
% developped for fMRI by Smith & Nichols, NeuroImage 44(2009), 83-98
% tfce = sum(extent(h)^E*height^H*dh)
% 
% Inputs
%   tfce_score = limo_tfce(type, data, channeighbstructmat)
%   tfce_score = limo_tfce(type, data, channeighbstructmat, updatebar, E, H, dh)
%
%   type: 1 = 1D, 2 = 2D, 3 = 3D. First dim must be channels when using a neighbor graph
%   data: t/F map or stack under H0 in last dim (boot/permutation)
%   channeighbstructmat: [] to use grid connectivity; otherwise neighbor graph
%   updatebar: 0/1 lightweight console prints every ~5 percent
%   E, H, dh: TFCE parameters (defaults E=0.5, H=2, dh=0.1)
%
% Outputs
%   tfce_score: TFCE transformed map(s) same shape as input map(s) excluding the H0 dim
%   thresholded_maps: per-height TFCE contributions (optional)
% 
% References
%
% Pernet, C., Latinus, M., Nichols, T.E., & Rousselet, G.A. (2015)
% Cluster-based computational methods for mass univariate analyses
% of event-related brain potentials/fields: a simulation study
% Journal Of Neuroscience Method 250, Pages 85–93
% <10.1016/j.jneumeth.2014.08.003>
%
% Pernet, Cyril; Rousselet, Guillaume (2014): Type 1 error rate using TFCE for ERP.
% figshare. http://dx.doi.org/10.6084/m9.figshare.1008325
%
% Cyril Pernet V5 20-08-2015
% use limo_findcluster which is faster (clustering speed x60)
% changed the integration from a loop to hist - thx to Bruno Giordano
% ------------------------------
%  Copyright (C) LIMO Team 2019
% 
% 2025 update: performance and robustness pass

%% parse args
if nargin < 3
    error('limo_tfce: not enough arguments')
elseif nargin == 3 || nargin == 4
    if nargin == 3, updatebar = 1; else, updatebar = varargin{4}; end
    E = 0.5; H = 2; dh = 0.1;
elseif nargin == 7
    updatebar = varargin{4}; E = varargin{5}; H = varargin{6}; dh = varargin{7};
else
    error('limo_tfce: wrong number of inputs')
end

type  = varargin{1};
data  = varargin{2};
G     = varargin{3};           % neighbor graph or []

thresholded_maps = [];

% sanitize neighbor graph once
if ~isempty(G)
    G = logical(G);
    G = G | G.';
    nd = size(G,1);
    G(1:nd+1:end) = false;
end

% use like-typed outputs
outClass = class(data);

%% dispatch by dimensionality and presence of H0 stack
switch type
case 1
    % data: [T] or [T x B]
    sz = size(data);
    if numel(sz) == 2 && sz(2) > 1
        subtype = 2;  % H0 stack
        T = sz(1); B = sz(2);
    else
        subtype = 1;  % single map
        T = numel(data); B = 1;
    end

    if subtype == 1
        tfce_score = tfce_one_map_1d(reshape(data, [T 1]), G, E, H, dh, updatebar);
    else
        tfce_score = zeros(T, B, outClass);
        if updatebar, fprintf('TFCE 1D %d maps: ', B); end
        for b = 1:B
            tfce_score(:,b) = tfce_one_map_1d(data(:,b), G, E, H, dh, 0);
            if updatebar && mod(b, max(1, floor(B/20))) == 0, fprintf('.'); end
        end
        if updatebar, fprintf(' %d/%d\n', B, B); end
    end

case 2
    % data: [X Y] or [X Y B]
    [X, Y, B] = size(data);
    if B == 1
        tfce_score = tfce_one_map_nd(data, G, E, H, dh, updatebar, 2);
    else
        tfce_score = zeros(X, Y, B, outClass);
        if updatebar, fprintf('TFCE 2D %d maps: ', B); end
        for b = 1:B
            tfce_score(:,:,b) = tfce_one_map_nd(data(:,:,b), G, E, H, dh, 0, 2);
            if updatebar && mod(b, max(1, floor(B/20))) == 0, fprintf('.'); end
        end
        if updatebar, fprintf(' %d/%d\n', B, B); end
    end

case 3
    % data: [X Y Z] or [X Y Z B]
    sz = size(data);
    if numel(sz) == 3
        tfce_score = tfce_one_map_nd(data, G, E, H, dh, updatebar, 3);
    else
        [X, Y, Z, B] = size(data);
        tfce_score = zeros(X, Y, Z, B, outClass);
        if updatebar, fprintf('TFCE 3D %d maps: ', B); end
        for b = 1:B
            tfce_score(:,:,:,b) = tfce_one_map_nd(data(:,:,:,b), G, E, H, dh, 0, 3);
            if updatebar && mod(b, max(1, floor(B/20))) == 0, fprintf('.'); end
        end
        if updatebar, fprintf(' %d/%d\n', B, B); end
    end

otherwise
    error('limo_tfce: unsupported type %d', type)
end

% thresholded_maps remains unset for speed; can be added back if needed

end


% ======================== helpers ========================

function tfce = tfce_one_map_1d(map, G, E, H, dh, updatebar)
map = double(map);
pos = max(map, 0);
neg = max(-map, 0);
if ~any(pos(:)) && ~any(neg(:)), tfce = zeros(size(map)); return, end

hpos = build_heights_one_sided(max(pos(:)), dh, 200);
hneg = build_heights_one_sided(max(neg(:)), dh, 200);

tfce_pos = tfce_accumulate_1d(pos, G, hpos, E, H, updatebar);
tfce_neg = tfce_accumulate_1d(neg, G, hneg, E, H, updatebar);
tfce = cast(tfce_pos + tfce_neg, 'like', map);
end

function tfce = tfce_one_map_nd(map, G, E, H, dh, updatebar, nd)
map = double(map);
pos = max(map, 0);
neg = max(-map, 0);
if ~any(pos(:)) && ~any(neg(:)), tfce = zeros(size(map)); return, end

hpos = build_heights_one_sided(max(pos(:)), dh, 200);
hneg = build_heights_one_sided(max(neg(:)), dh, 200);

tfce_pos = tfce_accumulate_nd(pos, G, hpos, E, H, updatebar, nd);
tfce_neg = tfce_accumulate_nd(neg, G, hneg, E, H, updatebar, nd);
tfce = cast(tfce_pos + tfce_neg, 'like', map);
end

function hvec = build_heights_one_sided(maxv, dh, cap)
% Heights start at 0 and go up to the positive maximum ONLY
if ~(isfinite(maxv) && maxv > 0), hvec = 0; return, end
if maxv > 1
    precision = min(cap, max(2, round(maxv / dh)));
else
    precision = min(cap, max(2, round(1 / dh)));  % finer grid for small ranges
end
hvec = linspace(0, maxv, precision);
end


function tfce = tfce_accumulate_1d(map, G, hvec, E, H, updatebar)
T = size(map,1);
tfce = zeros(T,1);
if all(map(:) == 0), return, end

print_every = max(1, floor(numel(hvec)/20));
for i = 1:numel(hvec)
    h = hvec(i);
    bw = map > h;                           % logical [T x 1]
    if isempty(G)
        % connectivity along the vector
        CC = bwconncomp(bw, 2);
        ext = zeros(T,1);
        for k = 1:CC.NumObjects
            idx = CC.PixelIdxList{k};
            ext(idx) = numel(idx);
        end
    else
        % neighbor graph in channels space; here T must match size(G,1)
        if size(G,1) ~= T
            error('limo_tfce: neighbor graph size does not match 1D length')
        end
        % build components via graph traversal
        ext = components_from_graph(bw, G);
    end
    tfce = tfce + (ext.^E) * (h^H) * step_size(hvec, i);
    if updatebar && mod(i, print_every) == 0
        fprintf('.');
    end
end
if updatebar, fprintf(' %d/%d\n', numel(hvec), numel(hvec)); end
end

function tfce = tfce_accumulate_nd(map, G, hvec, E, H, updatebar, nd)
sz = size(map);
tfce = zeros(sz);

if all(map(:) == 0), return, end

% pick grid connectivity when no graph
if isempty(G)
    if nd == 2
        conn = 4;
    else
        conn = 6;
    end
end

print_every = max(1, floor(numel(hvec)/20));
for i = 1:numel(hvec)
    h = hvec(i);
    bw = map > h;

    if isempty(G)
        CC = bwconncomp(bw, conn);
        ext = zeros(sz);
        for k = 1:CC.NumObjects
            idx = CC.PixelIdxList{k};
            ext(idx) = numel(idx);
        end
    else
        % When a channel neighbor graph is provided, we cluster across channels
        % at each time or time x freq slice. This assumes dim 1 = channels.
        if nd == 2
            % bw: [chan x time]
            [nChan, nTime] = size(bw);
            if size(G,1) ~= nChan
                error('limo_tfce: neighbor graph rows do not match channels')
            end
            ext = zeros(nChan, nTime);
            for t = 1:nTime
                ext(:,t) = components_from_graph(bw(:,t), G);
            end
        else
            % nd == 3: [chan x freq x time]
            [nChan, nF, nT] = size(bw);
            if size(G,1) ~= nChan
                error('limo_tfce: neighbor graph rows do not match channels')
            end
            ext = zeros(nChan, nF, nT);
            for f = 1:nF
                for t = 1:nT
                    ext(:,f,t) = components_from_graph(bw(:,f,t), G);
                end
            end
        end
    end

    tfce = tfce + (ext.^E) * (h^H) * step_size(hvec, i);
    if updatebar && mod(i, print_every) == 0
        fprintf('.');
    end
end
if updatebar, fprintf(' %d/%d\n', numel(hvec), numel(hvec)); end
end

function ds = step_size(hvec, i)
% piecewise constant step size between heights
if numel(hvec) == 1
    ds = 1;
elseif i == 1
    ds = hvec(2) - hvec(1);
else
    ds = hvec(i) - hvec(i-1);
end
if ~isfinite(ds) || ds <= 0
    ds = 0;
end
end

function extent = components_from_graph(active, G)
% active: logical [n x 1], G: logical [n x n] symmetric, zero diag
n = numel(active);
extent = zeros(n,1);
if ~any(active), return, end

% restrict to active nodes
idx = find(active);
sub = G(idx, idx);
% find components by BFS/DFS
visited = false(numel(idx),1);
for k = 1:numel(idx)
    if ~visited(k)
        comp = bfs(sub, k);
        visited(comp) = true;
        nodes = idx(comp);
        extent(nodes) = numel(nodes);
    end
end
end

function comp = bfs(A, s)
% breadth-first search on adjacency A starting at node s
n = size(A,1);
comp = false(n,1);
q = s;
comp(s) = true;
head = 1;
while head <= numel(q)
    v = q(head); head = head + 1;
    nbr = find(A(v,:));
    nbr = nbr(~comp(nbr));
    comp(nbr) = true;
    q = [q nbr]; %#ok<AGROW>
end
comp = find(comp);
end

