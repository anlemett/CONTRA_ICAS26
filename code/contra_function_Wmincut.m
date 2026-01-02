function Wmincut = contra_function_Wmincut(sector_pgon, T, B, cost_polygons)
% Compute min-cut between boundaries T and B in a sector, given obstacle polygons.

if isempty(sector_pgon) || ~isa(sector_pgon, 'polyshape') || sector_pgon.NumRegions == 0
    Wmincut = NaN;
    return
end

if isempty(T) || isempty(B) || size(T,2) < 2 || size(B,2) < 2
    Wmincut = NaN;
    return
end

[lon_c, lat_c] = centroid(sector_pgon);

% Sector polygon in local azimuthal-equidistant coordinates
[sector_y, sector_x] = contra_function_spherical_to_eq_azimuth( ...
    sector_pgon.Vertices(:,2), sector_pgon.Vertices(:,1), lat_c, lon_c);
sector_pgon_xy = polyshape(sector_x, sector_y);

% Boundary polylines in local coordinates
[T_y, T_x] = contra_function_spherical_to_eq_azimuth(T(:,2), T(:,1), lat_c, lon_c);
[B_y, B_x] = contra_function_spherical_to_eq_azimuth(B(:,2), B(:,1), lat_c, lon_c);

% Robustness: remove NaN/Inf separators (polyshape vertices may contain NaNs)
[T_x, T_y] = local_drop_nan_xy(T_x, T_y);
[B_x, B_y] = local_drop_nan_xy(B_x, B_y);

nw = numel(cost_polygons);

% Select obstacle parts that lie within the sector
count = 0;
x_w = {};
y_w = {};

for i = 1:nw
    try
        polyout_ll = intersect(sector_pgon, cost_polygons(i));
    catch
        continue
    end

    if polyout_ll.NumRegions > 0
        w_vertices = cost_polygons(i).Vertices;
        w_vertices = w_vertices(all(isfinite(w_vertices),2),:); % drop NaN rows
        if size(w_vertices,1) < 3
            continue
        end

        [y, x] = contra_function_spherical_to_eq_azimuth(w_vertices(:,2), w_vertices(:,1), lat_c, lon_c);
        w_pgon_xy = polyshape(x, y);

        polyout_xy = intersect(sector_pgon_xy, w_pgon_xy);

        if polyout_xy.NumRegions > 0
            V = polyout_xy.Vertices;
            V = V(all(isfinite(V),2),:); % drop NaN separators (multi-ring)
            if size(V,1) < 2
                continue
            end
            count = count + 1;
            x_w{count} = V(:,1);
            y_w{count} = V(:,2);
        end
    end
end

NW = count;          % number of obstacle cells intersecting the sector (in XY)
NG = 2 + NW;         % nodes: 1=T, 2..NW+1=obstacles, NG=B

% Graph edges between all node pairs
nodes_comb = nchoosek(1:NG, 2);
sg = nodes_comb(:,1).';
tg = nodes_comb(:,2).';

weights = zeros(size(sg));

% Direct T-B distance
weights(tg == NG & sg == 1) = boundary_to_boundary_fun(B_x, B_y, T_x, T_y);

% Distances involving obstacle nodes
for i = 1:NW
    xi = x_w{i}; yi = y_w{i};
    [xi, yi] = local_drop_nan_xy(xi, yi);

    % T -> obstacle i (node i+1) only if no intersection
    if isempty(polyxpoly(T_x, T_y, xi, yi))
        weights((sg == 1) & (tg == (i+1))) = boundary_to_boundary_fun(T_x, T_y, xi, yi);
    else
        weights((sg == 1) & (tg == (i+1))) = 0;
    end

    % obstacle i -> B
    if isempty(polyxpoly(B_x, B_y, xi, yi))
        weights((sg == (i+1)) & (tg == NG)) = boundary_to_boundary_fun(B_x, B_y, xi, yi);
    else
        weights((sg == (i+1)) & (tg == NG)) = 0;
    end

    % obstacle-obstacle distances
    for j = (i+1):NW
        xj = x_w{j}; yj = y_w{j};
        [xj, yj] = local_drop_nan_xy(xj, yj);
        weights((sg == (i+1)) & (tg == (j+1))) = boundary_to_boundary_fun(xj, yj, xi, yi);
    end
end

% Robustness: convert invalid weights to Inf (unusable edges)
bad = ~isfinite(weights) | (weights < 0);
weights(bad) = Inf;

% Build graph and compute shortest path from T (1) to B (NG)
G = graph(sg, tg, weights);
[~, Wmincut] = shortestpath(G, 1, NG);

end

%% ========================== Local helpers ===============================

function [x,y] = local_drop_nan_xy(x,y)
% Remove NaN/Inf separators from polyline coordinate arrays.
if isempty(x) || isempty(y)
    return
end
x = x(:); y = y(:);
m = isfinite(x) & isfinite(y);
x = x(m);
y = y(m);
end

function d = boundary_to_boundary_fun(x1, y1, x2, y2)
% Distance between two polylines (boundaries) in 2D.

[x1,y1] = local_drop_nan_xy(x1,y1);
[x2,y2] = local_drop_nan_xy(x2,y2);

N1 = numel(x1);
N2 = numel(x2);

if N1 < 1 || N2 < 1
    d = NaN;
    return
end

d = inf;

if N1 == 1
    pt = [x1, y1];
    for j = 1:max(1, N2-1)
        if N2 == 1
            w1 = [x2(1), y2(1)];
            w2 = w1;
        else
            w1 = [x2(j),   y2(j)];
            w2 = [x2(j+1), y2(j+1)];
        end
        d_new = point_to_line(pt, w1, w2);
        if isfinite(d_new)
            d = min(d, d_new);
        end
    end

elseif N2 == 1
    pt = [x2, y2];
    for i = 1:max(1, N1-1)
        if N1 == 1
            v1 = [x1(1), y1(1)];
            v2 = v1;
        else
            v1 = [x1(i),   y1(i)];
            v2 = [x1(i+1), y1(i+1)];
        end
        d_new = point_to_line(pt, v1, v2);
        if isfinite(d_new)
            d = min(d, d_new);
        end
    end

else
    for i = 1:(N1-1)
        v1 = [x1(i),   y1(i)];
        v2 = [x1(i+1), y1(i+1)];
        for j = 1:(N2-1)
            w1 = [x2(j),   y2(j)];
            w2 = [x2(j+1), y2(j+1)];
            d_new = line_to_line(v1, v2, w1, w2);
            if isfinite(d_new)
                d = min(d, d_new);
            end
        end
    end
end

if isinf(d)
    d = NaN;
end

end

function d = line_to_line(v1, v2, w1, w2)
% Min distance between two line segments in 2D.
d1 = point_to_line(v1, w1, w2);
d2 = point_to_line(v2, w1, w2);
d3 = point_to_line(w1, v1, v2);
d4 = point_to_line(w2, v1, v2);
d = min([d1 d2 d3 d4], [], 'omitnan');
end

function d = point_to_line(pt, v1, v2)
% Distance from point pt to segment v1-v2 in 2D.

a = pt - v1;
b = v2 - v1;

len = b*b.';
if len == 0
    x = v1;
else
    param = (a*b.')/len;
    if param < 0
        x = v1;
    elseif param > 1
        x = v2;
    else
        x = v1 + param*b;
    end
end

d = sqrt((x-pt)*(x-pt).');

end
