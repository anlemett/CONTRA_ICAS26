function [sector_ab, a_band, flows_j] = contra_function_flows_sector_k(k, main_sectors, adjacent_sectors)
% Build sector polygons per altitude band and admissible flows (triplets)
% for sector k based on shared borders with other (main + adjacent) sectors.
%
%  1) Adjacency detection uses "near boundary distance".
%  2) Each neighbor adjacency mask is reduced to ONE contiguous run (longest),
%     so masks cannot cover the whole boundary accidentally.
%  3) Top/Bottom construction is robust:
%     - primary: complement of (S U D) split into 2 components
%     - fallback A (when complement empty): T/B arcs between S and D runs
%     - fallback B (single complement component): split deterministically
%
% Diagnostics:
%   Set DEBUG_FLOWS=true to see skip reasons and counts.

%% ------------------ Tuning / Diagnostics ------------------
DEBUG_FLOWS = false;     % set true to print skip reasons + counts

% Adjacency tolerance in degrees (lon/lat).
% Typical: 1e-4..5e-4 depending on vertex density.
ADJ_TOL_DEG = 2e-4;

% Small gap closing and minimum run length (vertex units)
GAP_MAX_VERTICES = 1;
MIN_RUN_VERTICES = 1;

% Minimum points per edge polyline
MIN_EDGE_POINTS_SD = 2;  % S and D must have >=2 points
MIN_EDGE_POINTS_TB = 2;  % T and B must have >=2 points (for polyxpoly/mincut)

%% ------------------------------------------------------------------------
nSec = numel(main_sectors);       % number of main sectors
nAdj = numel(adjacent_sectors);   % number of adjacent sectors

sector_k = main_sectors{k};
nk = numel(sector_k);

% Determine altitude limits for sector k (FL)
low_lim = inf;
upp_lim = -inf;
for i = 1:nk
    low_lim = min(low_lim, sector_k(i).properties.LOWER_LIMIT_VALUE);
    upp_lim = max(upp_lim, sector_k(i).properties.UPPER_LIMIT_VALUE);
end
if ~isfinite(low_lim) || ~isfinite(upp_lim) || upp_lim <= low_lim
    error('Invalid altitude limits for sector k=%d.', k);
end

% Altitude bands within the sector
a_band = (0:10:600)'; % FL
a_band = a_band((a_band > low_lim) & (a_band < upp_lim));
nab = numel(a_band);

sector_ab = cell(nab, 1);
flows_j   = cell(nab, 1);

for i = 1:nab

    [sectors_pgon, adj_pgon] = contra_function_sector_adjacent_pgon_a_band( ...
        main_sectors, adjacent_sectors, a_band(i));

    if k > numel(sectors_pgon) || isempty(sectors_pgon{k}) || sectors_pgon{k}.NumRegions == 0
        sector_ab{i} = polyshape.empty(0,1);
        flows_j{i}   = [];
        continue
    end

    sector_ab{i} = sectors_pgon{k};

    Vk = sectors_pgon{k}.Vertices;
    if isempty(Vk) || size(Vk,1) < 3
        flows_j{i} = [];
        continue
    end

    % Remove NaN separators if any
    Vk = Vk(all(isfinite(Vk),2),:);
    nV = size(Vk,1);
    if nV < 3
        flows_j{i} = [];
        continue
    end

    edge_data = false(nV, nSec + nAdj);

    % ----------------- Main sectors adjacency (near boundary) -----------------
    n_vec = 1:nSec;
    n_vec(n_vec == k) = [];
    for n = n_vec
        if n <= numel(sectors_pgon) && ~isempty(sectors_pgon{n}) && sectors_pgon{n}.NumRegions > 0
            Vn = sectors_pgon{n}.Vertices;
            Vn = Vn(all(isfinite(Vn),2),:);
            if size(Vn,1) < 2, continue; end
            m = local_near_poly_boundary(Vk, Vn, ADJ_TOL_DEG);
            m = local_keep_longest_run(m);     % IMPORTANT: prevent "almost all boundary"
            edge_data(:, n) = m;
        end
    end

    % ---------------- Adjacent sectors adjacency (near boundary) --------------
    for n = 1:nAdj
        if n <= numel(adj_pgon) && ~isempty(adj_pgon{n}) && adj_pgon{n}.NumRegions > 0
            Va = adj_pgon{n}.Vertices;
            Va = Va(all(isfinite(Va),2),:);
            if size(Va,1) < 2, continue; end
            m = local_near_poly_boundary(Vk, Va, ADJ_TOL_DEG);
            m = local_keep_longest_run(m);     % IMPORTANT: prevent "almost all boundary"
            edge_data(:, n + nSec) = m;
        end
    end

    % All candidate edges that exist (touch boundary somewhere)
    all_edges = 1:(nSec + nAdj);
    all_edges = all_edges(any(edge_data, 1));

    if DEBUG_FLOWS
        fprintf('[FL=%g] nnz(edge_data)=%d, cols_with_any=%d\n', ...
            a_band(i), nnz(edge_data), nnz(any(edge_data,1)));
        fprintf('[FL=%g] all_edges count = %d\n', a_band(i), numel(all_edges));
    end

    if numel(all_edges) <= 1
        flows_j{i} = [];
        continue
    end

    pairs  = nchoosek(all_edges, 2);
    npairs = size(pairs, 1);

    flows = repmat(struct('triplet', [], 'S', [], 'D', [], 'T', [], 'B', [], 'Omincut', []), ...
                   2*npairs, 1);

    nValid = 0;

    for p = 1:npairs
        e1 = pairs(p,1);
        e2 = pairs(p,2);

        maskS = logical(edge_data(:, e1));
        maskD = logical(edge_data(:, e2));

        [ok, S_edge, D_edge, T_edge, B_edge, reason] = local_build_SDTB_from_masks( ...
            Vk, maskS, maskD, GAP_MAX_VERTICES, MIN_RUN_VERTICES, ...
            MIN_EDGE_POINTS_SD, MIN_EDGE_POINTS_TB);

        if ~ok
            if DEBUG_FLOWS
                fprintf('  skip pair (%d,%d): %s\n', e1, e2, reason);
            end
            continue
        end

        % Forward flow
        nValid = nValid + 1;
        flows(nValid).triplet = [k, e1, e2];
        flows(nValid).S = S_edge;
        flows(nValid).D = D_edge;
        flows(nValid).T = T_edge;
        flows(nValid).B = B_edge;
        flows(nValid).Omincut = contra_function_Omincut(sectors_pgon{k}, T_edge, B_edge);

        % Inverse flow
        nValid = nValid + 1;
        flows(nValid).triplet = [k, e2, e1];
        flows(nValid).S = D_edge;
        flows(nValid).D = S_edge;
        flows(nValid).T = B_edge;
        flows(nValid).B = T_edge;
        flows(nValid).Omincut = flows(nValid-1).Omincut;
    end

    if DEBUG_FLOWS
        fprintf('[FL=%g] nValid flows = %d\n', a_band(i), nValid);
    end

    if nValid == 0
        flows_j{i} = [];
    else
        flows_j{i} = flows(1:nValid);
    end
end

end

%% ==================== Robust boundary partition ====================

function [ok, S_edge, D_edge, T_edge, B_edge, reason] = local_build_SDTB_from_masks( ...
    V, maskS, maskD, gapMax, minRun, minEdgePtsSD, minEdgePtsTB)

ok = false;
reason = '';

S_edge = zeros(0,2);
D_edge = zeros(0,2);
T_edge = zeros(0,2);
B_edge = zeros(0,2);

n = size(V,1);
if n < 3
    reason = 'boundary has <3 vertices';
    return
end

maskS = logical(maskS(:));
maskD = logical(maskD(:));
if numel(maskS) ~= n || numel(maskD) ~= n
    reason = 'mask length mismatch';
    return
end
if ~any(maskS) || ~any(maskD)
    reason = 'one of masks is empty';
    return
end

% Light denoise
maskS2 = local_close_gaps_circular(maskS, gapMax);
maskD2 = local_close_gaps_circular(maskD, gapMax);

maskS2 = local_remove_short_runs_circular(maskS2, minRun);
maskD2 = local_remove_short_runs_circular(maskD2, minRun);

if ~any(maskS2), maskS2 = maskS; end
if ~any(maskD2), maskD2 = maskD; end

S_runs = local_runs_circular(maskS2);
D_runs = local_runs_circular(maskD2);
if isempty(S_runs) || isempty(D_runs)
    reason = 'no runs found';
    return
end

% Choose run pair
[selS, selD] = local_select_best_run_pair(S_runs, D_runs, n);
if isempty(selS) || isempty(selD)
    [selS, selD] = local_select_longest_pair(S_runs, D_runs);
end

sS = selS(1); eS = selS(2);
sD = selD(1); eD = selD(2);

% Inclusive S/D masks
idxS = false(n,1); idxS(local_indices_cyclic_arc(n, sS, eS)) = true;
idxD = false(n,1); idxD(local_indices_cyclic_arc(n, sD, eD)) = true;

S_edge = local_clean_polyline(V(idxS,:));
D_edge = local_clean_polyline(V(idxD,:));

if size(S_edge,1) < minEdgePtsSD || size(D_edge,1) < minEdgePtsSD
    reason = 'S or D too short';
    return
end

% ---------------- TOP/BOTTOM primary: complement of (S U D) ----------------
idxU = idxS | idxD;
idxC = ~idxU;

if any(idxC)
    C_runs = local_runs_circular(idxC);
    if isempty(C_runs)
        reason = 'no complement runs';
        return
    end

    % Sort complement runs by length desc
    [~, ord] = sort(C_runs(:,3), 'descend');
    C_runs = C_runs(ord,:);

    if size(C_runs,1) >= 2
        idxT = false(n,1); idxT(local_indices_cyclic_arc(n, C_runs(1,1), C_runs(1,2))) = true;
        idxB = false(n,1); idxB(local_indices_cyclic_arc(n, C_runs(2,1), C_runs(2,2))) = true;
    else
        % Only one complement component -> split deterministically
        sC = C_runs(1,1); eC = C_runs(1,2);
        idxList = local_indices_cyclic_arc(n, sC, eC);
        L = numel(idxList);
        if L < 2
            reason = 'complement too small';
            return
        end
        mid = floor(L/2);
        idxT = false(n,1); idxT(idxList(1:mid)) = true;
        idxB = false(n,1); idxB(idxList(mid+1:end)) = true;
    end

    T_edge = local_clean_polyline(V(idxT,:));
    B_edge = local_clean_polyline(V(idxB,:));

    if size(T_edge,1) < minEdgePtsTB || size(B_edge,1) < minEdgePtsTB
        reason = 'T or B too short after complement split';
        return
    end

    ok = true;
    reason = 'ok';
    return
end

% ---------------- FALLBACK A: complement empty => arcs between S and D -----
% This handles the case that previously caused "S U D covers entire boundary".
% Define:
%   Top    = arc from end(S) to start(D)
%   Bottom = arc from end(D) to start(S)
idxT = false(n,1); idxT(local_indices_cyclic_arc(n, eS, sD)) = true;
idxB = false(n,1); idxB(local_indices_cyclic_arc(n, eD, sS)) = true;

% Ensure we do not accidentally make T/B identical or degenerate
T_edge = local_clean_polyline(V(idxT,:));
B_edge = local_clean_polyline(V(idxB,:));

if size(T_edge,1) < minEdgePtsTB || size(B_edge,1) < minEdgePtsTB
    % FALLBACK B: deterministic split of full boundary into two halves
    idxList = (1:n).';
    mid = floor(n/2);
    idxT = false(n,1); idxT(idxList(1:mid)) = true;
    idxB = false(n,1); idxB(idxList(mid+1:end)) = true;

    T_edge = local_clean_polyline(V(idxT,:));
    B_edge = local_clean_polyline(V(idxB,:));

    if size(T_edge,1) < minEdgePtsTB || size(B_edge,1) < minEdgePtsTB
        reason = 'T/B degenerate even after fallbacks';
        return
    end
end

ok = true;
reason = 'ok';

end

%% ==================== Adjacency by near-boundary distance ====================

function tf = local_near_poly_boundary(P, Vpoly, tol)
% tf(i)=true if point P(i,:) is within tol (lon/lat degrees) of polyline Vpoly.
P = P(all(isfinite(P),2),:);
Vpoly = Vpoly(all(isfinite(Vpoly),2),:);

nP = size(P,1);
tf = false(nP,1);

if nP == 0 || size(Vpoly,1) < 2
    return
end

% Ensure closed boundary polyline
if any(Vpoly(1,:) ~= Vpoly(end,:))
    Vpoly(end+1,:) = Vpoly(1,:);
end

for ii = 1:nP
    px = P(ii,1); py = P(ii,2);
    dmin = inf;

    for jj = 1:(size(Vpoly,1)-1)
        x1 = Vpoly(jj,1);   y1 = Vpoly(jj,2);
        x2 = Vpoly(jj+1,1); y2 = Vpoly(jj+1,2);

        d = local_point_to_segment_dist(px, py, x1, y1, x2, y2);
        if d < dmin
            dmin = d;
        end
        if dmin <= tol
            break
        end
    end

    tf(ii) = (dmin <= tol);
end
end

function d = local_point_to_segment_dist(px, py, x1, y1, x2, y2)
% Euclidean distance in lon/lat degrees (adequate for small regions).
vx = x2 - x1; vy = y2 - y1;
wx = px - x1; wy = py - y1;

c1 = wx*vx + wy*vy;
if c1 <= 0
    d = hypot(px - x1, py - y1);
    return
end

c2 = vx*vx + vy*vy;
if c2 <= eps
    d = hypot(px - x1, py - y1);
    return
end

t = c1 / c2;
if t >= 1
    d = hypot(px - x2, py - y2);
else
    projx = x1 + t*vx;
    projy = y1 + t*vy;
    d = hypot(px - projx, py - projy);
end
end

%% ==================== Mask compression: keep only longest run ====================

function tf = local_keep_longest_run(tf)
% Reduce a circular logical vector to a single (longest) TRUE run.
tf = logical(tf(:));
n = numel(tf);
if n == 0 || ~any(tf)
    return
end
if all(tf)
    % This is pathological for adjacency; keep a conservative half instead
    % to avoid "touches everywhere" degeneracy.
    out = false(n,1);
    out(1:floor(n/2)) = true;
    tf = out;
    return
end

runs = local_runs_circular(tf);
if isempty(runs)
    tf(:) = false;
    return
end

[~, idx] = max(runs(:,3));
s = runs(idx,1);
e = runs(idx,2);

out = false(n,1);
out(local_indices_cyclic_arc(n, s, e)) = true;
tf = out;
end

%% ==================== Circular run utilities ====================

function runs = local_runs_circular(tf)
% Return all contiguous TRUE runs in a circular logical vector.
% runs: [K x 3] rows [start, end, length] inclusive.
tf = logical(tf(:));
n = numel(tf);
runs = zeros(0,3);

if n == 0 || ~any(tf)
    return
end
if all(tf)
    runs = [1, n, n];
    return
end

d = diff([tf; tf(1)]);
starts = find(d == 1) + 1;
ends   = find(d == -1);

starts(starts == n+1) = 1;

starts = starts(:);
ends   = ends(:);

if tf(1) && (isempty(starts) || starts(1) ~= 1)
    starts = [1; starts];
end

% If mismatch, fallback to a single longest run
if numel(starts) ~= numel(ends)
    [s,e] = local_longest_true_run_fallback(tf);
    if ~isempty(s)
        runs = [s,e, local_cyclic_len(n, s, e)];
    end
    return
end

K = numel(starts);
runs = zeros(K,3);
for k = 1:K
    s = starts(k);
    e = ends(k);
    L = local_cyclic_len(n, s, e);
    runs(k,:) = [s,e,L];
end
end

function [selS, selD] = local_select_best_run_pair(S_runs, D_runs, n)
% Pick run pair that yields stable partition.
% Objective: maximize min(lenT,lenB), tie-break minimize max(lenT,lenB).
selS = [];
selD = [];

bestMin = -inf;
bestMax = inf;

for i = 1:size(S_runs,1)
    sS = S_runs(i,1); eS = S_runs(i,2);
    for j = 1:size(D_runs,1)
        sD = D_runs(j,1); eD = D_runs(j,2);

        lenT = local_cyclic_len(n, eS, sD);
        lenB = local_cyclic_len(n, eD, sS);

        m = min(lenT, lenB);
        M = max(lenT, lenB);

        if (m > bestMin) || (m == bestMin && M < bestMax)
            bestMin = m;
            bestMax = M;
            selS = [sS, eS];
            selD = [sD, eD];
        end
    end
end
end

function [selS, selD] = local_select_longest_pair(S_runs, D_runs)
[~, iS] = max(S_runs(:,3));
[~, iD] = max(D_runs(:,3));
selS = S_runs(iS,1:2);
selD = D_runs(iD,1:2);
end

function tf = local_close_gaps_circular(tf, gapMax)
% Close FALSE gaps of length <= gapMax between TRUE runs on a circular vector.
tf = logical(tf(:));
n = numel(tf);
if gapMax <= 0 || n == 0 || ~any(tf)
    return
end

f = ~tf;
runsF = local_runs_circular(f);
if isempty(runsF)
    return
end

tf2 = tf;
for k = 1:size(runsF,1)
    if runsF(k,3) <= gapMax
        idx = local_indices_cyclic_arc(n, runsF(k,1), runsF(k,2));
        tf2(idx) = true;
    end
end
tf = tf2;
end

function tf = local_remove_short_runs_circular(tf, minRun)
% Remove TRUE runs shorter than minRun.
tf = logical(tf(:));
n = numel(tf);
if minRun <= 1 || n == 0 || ~any(tf)
    return
end

runs = local_runs_circular(tf);
if isempty(runs)
    return
end

out = false(n,1);
for k = 1:size(runs,1)
    if runs(k,3) >= minRun
        idx = local_indices_cyclic_arc(n, runs(k,1), runs(k,2));
        out(idx) = true;
    end
end
tf = out;
end

function idx = local_indices_cyclic_arc(n, i0, i1)
% Inclusive forward arc i0->i1 on circular indices.
i0 = local_wrap_index(i0, n);
i1 = local_wrap_index(i1, n);
if i0 <= i1
    idx = (i0:i1).';
else
    idx = [ (i0:n).'; (1:i1).' ];
end
end

function L = local_cyclic_len(n, i0, i1)
% Inclusive forward arc length i0->i1 on a circle.
i0 = local_wrap_index(i0, n);
i1 = local_wrap_index(i1, n);
if i0 <= i1
    L = i1 - i0 + 1;
else
    L = (n - i0 + 1) + i1;
end
end

function i = local_wrap_index(i, n)
i = mod(i-1, n) + 1;
end

function P = local_clean_polyline(P)
% Remove non-finite rows and consecutive duplicates
if isempty(P)
    return
end
P = P(all(isfinite(P),2),:);
if size(P,1) < 2
    return
end
d = diff(P,1,1);
keep = [true; any(abs(d) > 0, 2)];
P = P(keep,:);
end

function [sRun, eRun] = local_longest_true_run_fallback(tf)
n = numel(tf);
sRun = []; eRun = [];
if n == 0 || ~any(tf)
    return
end
if all(tf)
    sRun = 1; eRun = n; return
end

tf2 = [tf; tf];
bestLen = 0; bestS = []; bestE = [];

ii = 1;
while ii <= n
    if tf2(ii)
        jj = ii;
        while jj <= ii+n-1 && tf2(jj)
            jj = jj + 1;
        end
        len = jj - ii;
        if len > bestLen
            bestLen = len;
            bestS = mod(ii-1,n) + 1;
            bestE = mod(jj-2,n) + 1;
        end
        ii = jj;
    else
        ii = ii + 1;
    end
end

sRun = bestS;
eRun = bestE;
end
