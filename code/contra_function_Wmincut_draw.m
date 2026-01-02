function Wmincut = contra_function_Wmincut_draw(sector_pgon, T, B, cost_polygons)
%contra_function_Wmincut_draw  Draw version of mincut (map + graph).
% Computes the same shortest-path mincut value as contra_function_Wmincut,
% and draws:
%   1) Equal-azimuth map: sector, obstacles, Top (blue), Bottom (red),
%      and the mincut path as RED DASHED segments that connect
%      boundary-to-boundary closest points.
%   2) Graph view with highlighted shortest path.
%
% Inputs:
%   sector_pgon   : polyshape (lon/lat)
%   T, B          : [N x 2] edges in lon/lat (columns: lon, lat)
%   cost_polygons : polyshape array (lon/lat) of obstacles (may be empty)
%
% Output:
%   Wmincut       : shortest path length (same scalar as contra_function_Wmincut)

% ---------- Projection to equal-azimuth ----------
[lon_c, lat_c] = centroid(sector_pgon);

[sector_y, sector_x] = contra_function_spherical_to_eq_azimuth( ...
    sector_pgon.Vertices(:,2), sector_pgon.Vertices(:,1), lat_c, lon_c);
sector_pgon_xy = polyshape(sector_x, sector_y);

[T_y, T_x] = contra_function_spherical_to_eq_azimuth(T(:,2), T(:,1), lat_c, lon_c);
[B_y, B_x] = contra_function_spherical_to_eq_azimuth(B(:,2), B(:,1), lat_c, lon_c);

% ---------- Figure 1: Map ----------
f1 = figure('Name','Mincut in equal-azimuth projection', 'NumberTitle','off');
clf(f1);
axes('Parent', f1);
hold on; grid on;
axis equal;
title('Mincut in equal-azimuth projection');
xlabel('X (equal-azimuth)');
ylabel('Y (equal-azimuth)');

% Sector outline
plot(sector_pgon_xy, 'FaceAlpha', 0.0, 'LineWidth', 1.5);

% Obstacles inside sector (projected/intersected)
nw = numel(cost_polygons);
x_w = {}; y_w = {};
for i = 1:nw
    polyout_ll = intersect(sector_pgon, cost_polygons(i));
    if polyout_ll.NumRegions <= 0
        continue
    end

    w_vertices = cost_polygons(i).Vertices;
    [wy, wx] = contra_function_spherical_to_eq_azimuth( ...
        w_vertices(:,2), w_vertices(:,1), lat_c, lon_c);
    w_pgon_xy = polyshape(wx, wy);

    polyout_xy = intersect(sector_pgon_xy, w_pgon_xy);
    if polyout_xy.NumRegions <= 0
        continue
    end

    plot(polyout_xy, ...
        'EdgeColor', [199, 0, 57]/255, ...
        'FaceColor', [199, 0, 57]/255, ...
        'FaceAlpha', 0.25, ...
        'LineWidth', 1.0);

    x_w{end+1} = polyout_xy.Vertices(:,1);
    y_w{end+1} = polyout_xy.Vertices(:,2);
end

NW = numel(x_w);       % obstacles that actually intersect the sector (in XY)
NG = 2 + NW;           % nodes: 1=T, 2..NW+1=obstacles, NG=B

% ---------- Build complete graph (undirected) ----------
nodes_comb = nchoosek(1:NG, 2);
sg = nodes_comb(:,1)';  % source node indices (always < tg)
tg = nodes_comb(:,2)';  % target node indices

weights = zeros(size(sg));

% For plotting the chosen mincut as boundary-to-boundary segments, store endpoints:
p1_xy = cell(size(sg)); % [x y] on boundary of sg-node side
p2_xy = cell(size(sg)); % [x y] on boundary of tg-node side

% Compute weights and closest boundary points for each pair (u,v)
for e = 1:numel(sg)
    u = sg(e);
    v = tg(e);

    if u == 1 && v == NG
        % T - B
        [w, pB, pT] = boundary_to_boundary_fun(B_x, B_y, T_x, T_y);
        weights(e) = w;
        % store endpoints in consistent order with nodes (u=1 is T, v=NG is B)
        p1_xy{e} = pT; % on T
        p2_xy{e} = pB; % on B

    elseif u == 1 && (v >= 2 && v <= NW+1)
        % T - obstacle(v-1)
        xo = x_w{v-1}; yo = y_w{v-1};
        [w, pT, pO] = boundary_to_boundary_fun(T_x, T_y, xo, yo);
        weights(e) = w;
        p1_xy{e} = pT; % on T
        p2_xy{e} = pO; % on obstacle

    elseif (u >= 2 && u <= NW+1) && v == NG
        % obstacle(u-1) - B
        xo = x_w{u-1}; yo = y_w{u-1};
        [w, pB, pO] = boundary_to_boundary_fun(B_x, B_y, xo, yo);
        weights(e) = w;
        % store endpoints in consistent order with nodes (u is obstacle, v is B)
        p1_xy{e} = pO; % on obstacle
        p2_xy{e} = pB; % on B

    else
        % obstacle - obstacle
        xo1 = x_w{u-1}; yo1 = y_w{u-1};
        xo2 = x_w{v-1}; yo2 = y_w{v-1};
        [w, pO1, pO2] = boundary_to_boundary_fun(xo1, yo1, xo2, yo2);
        weights(e) = w;
        p1_xy{e} = pO1;
        p2_xy{e} = pO2;
    end
end

G = graph(sg, tg, weights);
[short_path, Wmincut] = shortestpath(G, 1, NG);

% ---------- Plot Top/Bottom (ALWAYS after obstacles, before mincut) ----------
% Top (blue)
plot(T_x, T_y, '-', 'Color', 'b', 'LineWidth', 4);
text(mean(T_x), mean(T_y), 'T', 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 14);

% Bottom (red)
plot(B_x, B_y, '-', 'Color', 'r', 'LineWidth', 4);
text(mean(B_x), mean(B_y), 'B', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 14);

% -------------- Plot mincut on the MAP as RED DASHED segments ---------------
if numel(short_path) >= 2
    for k = 1:(numel(short_path)-1)
        u = short_path(k);
        v = short_path(k+1);

        % locate undirected edge index in (sg,tg) (stored with sg<tg)
        uu = min(u,v);
        vv = max(u,v);
        e = find((sg == uu) & (tg == vv), 1, 'first');
        if isempty(e) || isempty(p1_xy{e}) || isempty(p2_xy{e})
            continue
        end

        % p1_xy{e} belongs to node sg(e) side, p2_xy{e} belongs to node tg(e) side.
        % If we traversed edge in reverse, swap endpoints so it still draws correctly.
        P1 = p1_xy{e};
        P2 = p2_xy{e};
        if u ~= sg(e)
            % traversing from tg -> sg
            tmp = P1; P1 = P2; P2 = tmp;
        end

        plot([P1(1), P2(1)], [P1(2), P2(2)], 'r--', 'LineWidth', 2.5);
    end
end

% ---------- Figure 2: Graph ----------
f2 = figure('Name','Graph for shortest path', 'NumberTitle','off');
clf(f2);
axes('Parent', f2);
title('Graph for shortest path');
H = plot(G, 'Layout', 'layered', 'EdgeLabel', G.Edges.Weight);
highlight(H, short_path, 'EdgeColor', 'r', 'LineWidth', 2);

% Nicer node labels (T, obstacles, B)
try
    nodeNames = strings(1, NG);
    nodeNames(1) = "T";
    for i = 1:NW
        nodeNames(1+i) = "O" + i;
    end
    nodeNames(NG) = "B";
    H.NodeLabel = cellstr(nodeNames);
catch
    % ignore if plotting backend does not support
end

end


% ======================================================================
% Geometry helpers (used by both mincut and draw)
% ======================================================================

function [d, p1, p2] = boundary_to_boundary_fun(x1, y1, x2, y2)
% Returns:
%   d  : minimum distance between boundaries
%   p1 : closest point [x y] on boundary 1
%   p2 : closest point [x y] on boundary 2

N1 = length(x1);
N2 = length(x2);

d = inf;
p1 = [NaN NaN];
p2 = [NaN NaN];

if N1 == 1
    pt = [x1, y1];
    p1 = pt;
    for j = 1:(N2-1)
        w1 = [x2(j), y2(j)];
        w2 = [x2(j+1), y2(j+1)];
        [d_new, x] = point_to_line(pt, w1, w2);
        if d_new < d
            d = d_new;
            p2 = x;
        end
    end

elseif N2 == 1
    pt = [x2, y2];
    p2 = pt;
    for i = 1:(N1-1)
        v1 = [x1(i), y1(i)];
        v2 = [x1(i+1), y1(i+1)];
        [d_new, x] = point_to_line(pt, v1, v2);
        if d_new < d
            d = d_new;
            p1 = x;
        end
    end

else
    for i = 1:(N1-1)
        v1 = [x1(i), y1(i)];
        v2 = [x1(i+1), y1(i+1)];
        for j = 1:(N2-1)
            w1 = [x2(j), y2(j)];
            w2 = [x2(j+1), y2(j+1)];
            [d_new, p1_new, p2_new] = line_to_line(v1, v2, w1, w2);
            if d_new < d
                d  = d_new;
                p1 = p1_new;
                p2 = p2_new;
            end
        end
    end
end
end

function [d, x1, x2] = line_to_line(v1, v2, w1, w2)
% Minimum distance between two 2D segments (v1-v2) and (w1-w2)

x1v = zeros(4,2); x1v(1,:) = v1; x1v(2,:) = v2;
x2v = zeros(4,2); x2v(3,:) = w1; x2v(4,:) = w2;

dv = zeros(1,4);

[dv(1), x2v(1,:)] = point_to_line(v1, w1, w2);
[dv(2), x2v(2,:)] = point_to_line(v2, w1, w2);
[dv(3), x1v(3,:)] = point_to_line(w1, v1, v2);
[dv(4), x1v(4,:)] = point_to_line(w2, v1, v2);

[d, idx] = min(dv);

x1 = x1v(idx,:);
x2 = x2v(idx,:);
end

function [d, x] = point_to_line(pt, v1, v2)
% Distance from point pt to segment v1-v2, and closest point x on segment

a = pt - v1;
b = v2 - v1;

len = b*b';
if len <= eps
    x = v1;
    d = sqrt((x-pt)*(x-pt)');
    return
end

param = (a*b')/len;

if param < 0
    x = v1;
elseif param > 1
    x = v2;
else
    x = v1 + param*b;
end

d = sqrt((x-pt)*(x-pt)');
end
