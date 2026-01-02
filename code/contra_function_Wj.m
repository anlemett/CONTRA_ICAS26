function [Wj, Wij, total_ac] = contra_function_Wj(t_ini, t_fin, p_in, p_out, a_band, flows_j, AC_in)

% Collect all possible triplets robustly (avoid empty cells)
nonEmpty = ~cellfun(@isempty, flows_j);
if ~any(nonEmpty)
    nab = numel(a_band);
    Wj = zeros(1,0);
    Wij = zeros(nab,0);
    total_ac = 0;
    return
end

aux = [flows_j{nonEmpty}];
aux2 = vertcat(aux.triplet);
all_triplets = unique(aux2, 'rows');
[ntri, ~] = size(all_triplets);

nab = numel(a_band);
Wij = zeros(nab, ntri);

% Aircraft that fly sector k in the time interval [t_ini, t_fin]
ac_idx = (t_ini <= p_out(:,1)) & (t_fin >= p_in(:,1));

total_ac = sum(ac_idx);

if total_ac == 0
    Wj = zeros(1, ntri);
    return
end

p_in  = p_in(ac_idx,:);
p_out = p_out(ac_idx,:);
AC_in = AC_in(ac_idx);

for a = 1:total_ac
    [a_band_idx, tri] = contra_function_find_ac(p_in(a,:), p_out(a,:), a_band, flows_j);

    if ~isempty(tri) && ~isempty(a_band_idx)
        tri_idx = find(ismember(all_triplets, tri, 'rows'), 1, 'first');
        if ~isempty(tri_idx)
            Wij(a_band_idx, tri_idx) = Wij(a_band_idx, tri_idx) + 1;
        end
    end
end

total_ac = sum(Wij, 'all');

if total_ac > 0
    Wj  = sum(Wij, 1) / total_ac;
    Wij = Wij / total_ac;
else
    Wj = zeros(1, ntri);
end

end

function [a_band_idx, tri] = contra_function_find_ac(p_in, p_out, a_band, flows_j)

a_band_idx = [];
tri = [];

t_in   = p_in(1);
t_out  = p_out(1);
lon_in = p_in(2); lon_out = p_out(2);
lat_in = p_in(3); lat_out = p_out(3);
h_in   = p_in(4); h_out   = p_out(4);

h = (h_in + h_out) / 2;

a_band_idx = find((a_band >= (h-5)) & (a_band <= (h+5)), 1, 'first');
if isempty(a_band_idx)
    return
end

flow = flows_j{a_band_idx};
if isempty(flow)
    return
end

all_triplets = vertcat(flow.triplet);
ntri = size(all_triplets, 1);

% Avoid polyfit issues for vertical/near-vertical segments
dx = (lon_out - lon_in);
dy = (lat_out - lat_in);

% Parameterize line segment as P(s) = P0 + s*(P1-P0), s in [0,1]
P0 = [lon_in, lat_in];
P1 = [lon_out, lat_out];

for j = 1:(ntri/2)

    % Intersect with source edge
    [xS, yS] = polyxpoly([P0(1), P1(1)], [P0(2), P1(2)], flow(j).S(:,1), flow(j).S(:,2));
    if isempty(xS)
        continue
    end

    % Intersect with destination edge
    [xD, yD] = polyxpoly([P0(1), P1(1)], [P0(2), P1(2)], flow(j).D(:,1), flow(j).D(:,2));
    if isempty(xD)
        continue
    end

    % Use the closest intersection point to the entry point to decide direction
    dS = min(vecnorm([xS - lon_in, yS - lat_in], 2, 2));
    dD = min(vecnorm([xD - lon_in, yD - lat_in], 2, 2));

    if dS < dD
        tri_idx = j;
    else
        tri_idx = j + (ntri/2);
    end

    tri = flow(tri_idx).triplet;
    return
end

end
