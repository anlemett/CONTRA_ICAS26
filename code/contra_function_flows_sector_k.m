function [sector_ab, a_band, flows_j] = contra_function_flows_sector_k(k, main_sectors, adjacent_sectors)
% Build sector polygons per altitude band and the set of admissible flows (triplets)
% for sector k based on shared borders with other (main + adjacent) sectors.
%
% Inputs:
%   k               - index of the main sector of interest (scalar)
%   main_sectors    - 1xN cell; each cell contains a struct array of subsectors
%   adjacent_sectors- 1xM cell; each cell contains a struct array of adjacent subsectors
%
% Outputs:
%   sector_ab - cell(nab,1), sector polyshape for sector k at each altitude band
%   a_band    - vector of altitude bands (FL) considered for sector k
%   flows_j   - cell(nab,1), each contains struct array of flow definitions

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

% Build for each altitude band
for i = 1:nab

    % Polygons at this altitude band
    [sectors_pgon, adj_pgon] = contra_function_sector_adjacent_pgon_a_band(main_sectors, adjacent_sectors, a_band(i));

    if k > numel(sectors_pgon) || isempty(sectors_pgon{k}) || sectors_pgon{k}.NumRegions == 0
        sector_ab{i} = polyshape.empty(0,1);
        flows_j{i}   = [];
        continue
    end

    sector_ab{i} = sectors_pgon{k};

    % Determine which other polygons touch the boundary of sector k
    Vk = sectors_pgon{k}.Vertices;
    if isempty(Vk)
        flows_j{i} = [];
        continue
    end

    edge_data = zeros(size(Vk,1), nSec + nAdj);

    % Main sectors (excluding k)
    n_vec = 1:nSec;
    n_vec(n_vec == k) = [];

    for n = n_vec
        if n <= numel(sectors_pgon) && ~isempty(sectors_pgon{n}) && sectors_pgon{n}.NumRegions > 0
            in = inpolygon(Vk(:,1), Vk(:,2), sectors_pgon{n}.Vertices(:,1), sectors_pgon{n}.Vertices(:,2));
            edge_data(:, n) = in;
        end
    end

    % Adjacent sectors
    for n = 1:nAdj
        if n <= numel(adj_pgon) && ~isempty(adj_pgon{n}) && adj_pgon{n}.NumRegions > 0
            in = inpolygon(Vk(:,1), Vk(:,2), adj_pgon{n}.Vertices(:,1), adj_pgon{n}.Vertices(:,2));
            edge_data(:, n + nSec) = in;
        end
    end

    % All candidate edges that exist (touch boundary somewhere)
    all_edges = 1:(nSec + nAdj);
    all_edges = all_edges(any(edge_data, 1));

    if numel(all_edges) <= 1
        flows_j{i} = [];
        continue
    end

    % Triplets among edges
    triplets  = nchoosek(all_edges, 2);
    ntriplets = size(triplets, 1);

    flows = repmat(struct('triplet', [], 'S', [], 'D', [], 'T', [], 'B', [], 'Omincut', []), 2*ntriplets, 1);

    for n = 1:ntriplets
        e1 = triplets(n,1);
        e2 = triplets(n,2);

        % Forward flow
        flows(n).triplet = [k, e1, e2];
        flows(n).S = create_edge_fun(sectors_pgon{k}, logical(edge_data(:, e1)));
        flows(n).D = create_edge_fun(sectors_pgon{k}, logical(edge_data(:, e2)));

        [T_index, B_index] = top_bottom_index_fun(edge_data(:, e1), edge_data(:, e2));
        flows(n).T = create_edge_fun(sectors_pgon{k}, T_index);
        flows(n).B = create_edge_fun(sectors_pgon{k}, B_index);

        flows(n).Omincut = contra_function_Omincut(sectors_pgon{k}, flows(n).T, flows(n).B);

        % Inverse flow
        flows(n+ntriplets).triplet = [k, e2, e1];
        flows(n+ntriplets).S = flows(n).D;
        flows(n+ntriplets).D = flows(n).S;
        flows(n+ntriplets).T = flows(n).B;
        flows(n+ntriplets).B = flows(n).T;
        flows(n+ntriplets).Omincut = flows(n).Omincut;
    end

    flows_j{i} = flows;
end

end

%% ==================== Local helpers ====================

function edge = create_edge_fun(sector, mask)
%CREATE_EDGE_FUN  Build an ordered boundary polyline from a cyclic logical mask.
% Returns ONE contiguous arc (the longest true-run) in correct vertex order.

    mask = logical(mask(:));
    V = sector.Vertices;
    n = size(V,1);

    if n == 0 || ~any(mask)
        edge = zeros(0,2);
        return
    end

    % If only one vertex is selected, return it
    if nnz(mask) == 1
        edge = V(mask,:);
        return
    end

    % Find longest contiguous TRUE run on a circular list
    [sRun, eRun] = local_longest_true_run(mask);

    if isempty(sRun) || isempty(eRun)
        edge = V(mask,:); % safe fallback
        return
    end

    % Build ordered indices along boundary (inclusive), respecting wrap
    if sRun <= eRun
        idx = sRun:eRun;
    else
        idx = [sRun:n, 1:eRun];
    end

    edge = V(idx,:);
end


function [T_index, B_index] = top_bottom_index_fun(indexS, indexD)
%TOP_BOTTOM_INDEX_FUN  Robustly build "Top" and "Bottom" boundary index masks.
%
% Inputs
%   indexS, indexD : logical vectors (length n) indicating vertices belonging to
%                    the two "edge" regions (source/destination adjacency).
%
% Outputs
%   T_index, B_index : logical vectors (length n) selecting two complementary
%                      boundary arcs on the cyclic polygon vertex list.

    % --- normalize ---
    indexS = logical(indexS(:));
    indexD = logical(indexD(:));
    n = numel(indexS);

    T_index = false(n,1);
    B_index = false(n,1);

    if n < 2 || ~any(indexS) || ~any(indexD)
        return
    end

    % --- get circular runs (start/end indices) and pick the longest run ---
    [sS, eS] = local_longest_true_run(indexS);
    [sD, eD] = local_longest_true_run(indexD);

    % If something went wrong, exit safely
    if isempty(sS) || isempty(eS) || isempty(sD) || isempty(eD)
        return
    end

    % Top arc: from end(S-run) to start(D-run)
    % Bottom arc: from end(D-run) to start(S-run)
    ini_T = eS;   fin_T = sD;
    ini_B = eD;   fin_B = sS;

    % --- fill indices on a cyclic list, inclusive endpoints ---
    T_index = local_fill_cyclic_arc(n, ini_T, fin_T);
    B_index = local_fill_cyclic_arc(n, ini_B, fin_B);

end

% ============================ helpers ===================================

function [sRun, eRun] = local_longest_true_run(tf)
% Return start/end indices of the longest contiguous true-run in a circular vector.
    n = numel(tf);

    % Special case: all true
    if all(tf)
        sRun = 1; eRun = n;
        return
    end

    % Circular transitions
    d = diff([tf; tf(1)]);      % length n
    starts = find(d == 1) + 1;  % where run starts
    ends   = find(d == -1);     % where run ends

    % Wrap start index n+1 -> 1
    starts(starts == n+1) = 1;

    % If the vector starts inside a run, ensure a start at 1
    if tf(1) && (isempty(starts) || starts(1) ~= 1)
        starts = [1; starts(:)];
    else
        starts = starts(:);
    end
    ends = ends(:);

    % Defensive repair if mismatched counts (can happen in some degenerate masks)
    if numel(starts) ~= numel(ends)
        % fallback: choose the longest run in linear sense
        tf2 = [tf; tf];
        bestLen = 0; bestS = []; bestE = [];
        i = 1;
        while i <= n
            if tf2(i)
                j = i;
                while j <= i+n-1 && tf2(j)
                    j = j + 1;
                end
                len = j - i;
                if len > bestLen
                    bestLen = len;
                    bestS = mod(i-1,n) + 1;
                    bestE = mod(j-2,n) + 1;
                end
                i = j;
            else
                i = i + 1;
            end
        end
        sRun = bestS; eRun = bestE;
        return
    end

    % Compute cyclic lengths and select longest
    lens = mod(ends - starts, n) + 1;
    [~, k] = max(lens);
    sRun = starts(k);
    eRun = ends(k);
end

function mask = local_fill_cyclic_arc(n, i0, i1)
% Inclusive arc from i0 to i1 along forward direction in cyclic indexing.
    mask = false(n,1);

    if i0 <= i1
        mask(i0:i1) = true;
    else
        mask(i0:n) = true;
        mask(1:i1) = true;
    end
end
