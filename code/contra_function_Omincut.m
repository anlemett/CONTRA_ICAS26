function Omincut = contra_function_Omincut(sector_pgon, T, B)
% Compute obstacle-free min-cut (distance) between the bottom boundary B and
% top boundary T of a sector polygon at a given altitude band.

if isempty(sector_pgon) || ~isa(sector_pgon, 'polyshape') || sector_pgon.NumRegions == 0
    Omincut = NaN;
    return
end

if isempty(T) || isempty(B) || size(T,2) < 2 || size(B,2) < 2
    Omincut = NaN;
    return
end

[lon_c, lat_c] = centroid(sector_pgon);

% Convert to local azimuthal-equidistant coordinates around the centroid
[T_y, T_x] = contra_function_spherical_to_eq_azimuth(T(:,2), T(:,1), lat_c, lon_c);
[B_y, B_x] = contra_function_spherical_to_eq_azimuth(B(:,2), B(:,1), lat_c, lon_c);

weights = boundary_to_boundary_fun(B_x, B_y, T_x, T_y);

if ~isfinite(weights) || weights <= 0
    Omincut = NaN;
    return
end

G = graph(1, 2, weights);
[~, short_path] = shortestpath(G, 1, 2);

Omincut = short_path;

end

function d = boundary_to_boundary_fun(x1, y1, x2, y2)
% Distance between two polylines (boundaries) in 2D.

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
        d = min(d, d_new);
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
        d = min(d, d_new);
    end

else
    for i = 1:(N1-1)
        v1 = [x1(i),   y1(i)];
        v2 = [x1(i+1), y1(i+1)];
        for j = 1:(N2-1)
            w1 = [x2(j),   y2(j)];
            w2 = [x2(j+1), y2(j+1)];
            d_new = line_to_line(v1, v2, w1, w2);
            d = min(d, d_new);
        end
    end
end

end

function d = line_to_line(v1, v2, w1, w2)
% Min distance between two line segments (v1-v2) and (w1-w2) in 2D.

d1 = point_to_line(v1, w1, w2);
d2 = point_to_line(v2, w1, w2);
d3 = point_to_line(w1, v1, v2);
d4 = point_to_line(w2, v1, v2);

d = min([d1 d2 d3 d4]);

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
