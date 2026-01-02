function [p_in, p_out, AC_in] = contra_function_p_in_out(AC, sector_ab_k, a_band)
% Computes entry (p_in) and exit (p_out) points of aircraft trajectories AC(a).WP
% through a sector defined by polygons sector_ab_k over altitude bands a_band.
%
% AC(a).WP is assumed to be [t, lon, lat, alt_ft_or_FL_units_consistent_with_a_band].
%
%  1) explicit boundary intersection (polyxpoly) + inside checks (inpolygon)
%  2) If a trajectory is fully inside the sector horizontally (no boundary crossing),
%     it is still kept and we set:
%        p_in  = first waypoint
%        p_out = last waypoint
%  3) Guard against empty unions and NaN separators in polyshape vertices.

Nac = numel(AC);

% ----------------- robust union of sector polygons -----------------
sector_big = polyshape.empty(0,1);
for ii = 1:numel(sector_ab_k)
    if isempty(sector_ab_k{ii}) || ~isa(sector_ab_k{ii},'polyshape') || sector_ab_k{ii}.NumRegions==0
        continue
    end
    if isempty(sector_big)
        sector_big = sector_ab_k{ii};
    else
        try
            sector_big = union(sector_big, sector_ab_k{ii});
        catch
            % ignore bad unions
        end
    end
end

if isempty(sector_big) || sector_big.NumRegions == 0
    % No valid sector geometry -> nothing can be inside
    p_in  = nan(0,4);
    p_out = nan(0,4);
    AC_in = AC([]);
    return
end

% Clean sector boundary vertices (may contain NaN separators)
Vsec = sector_big.Vertices;
Vsec = Vsec(all(isfinite(Vsec),2),:);
if size(Vsec,1) < 3
    p_in  = nan(0,4);
    p_out = nan(0,4);
    AC_in = AC([]);
    return
end
% Close boundary polyline for polyxpoly
if any(Vsec(1,:) ~= Vsec(end,:))
    Vsec(end+1,:) = Vsec(1,:);
end

in_sector  = false(Nac,1);

for a = 1:Nac

    if ~isfield(AC(a), 'WP') || isempty(AC(a).WP) || size(AC(a).WP,2) < 4
        continue
    end

    WP = AC(a).WP;

    % Horizontal polyline
    traj = [WP(:,2), WP(:,3)];
    traj = local_drop_consecutive_duplicates(traj);

    if size(traj,1) < 2
        continue
    end

    % Horizontal checks: any point inside OR boundary intersection
    inside_any = false;
    try
        inside_any = any(inpolygon(traj(:,1), traj(:,2), Vsec(:,1), Vsec(:,2)));
    catch
        inside_any = false;
    end

    hits_boundary = false;
    try
        [xi, yi] = polyxpoly(traj(:,1), traj(:,2), Vsec(:,1), Vsec(:,2));
        hits_boundary = ~isempty(xi);
    catch
        hits_boundary = false;
    end

    ok_horizontal = inside_any || hits_boundary;

    % Vertical check
    max_alt = max(WP(:,4));
    min_alt = min(WP(:,4));
    ok_vertical = ~((max_alt <= (a_band(1) - 5)) || (min_alt > (a_band(end) + 5)));

    if ok_horizontal && ok_vertical
        in_sector(a) = true;
    end
end

AC_in = AC(in_sector);

% - keep classic "cross boundary" cases
% - ALSO keep "fully inside" cases (start and end inside, no boundary crossing)
index_pass = contra_function_pass_through(AC_in, sector_big, Vsec, a_band);
AC_in = AC_in(index_pass);

Nac = numel(AC_in);

p_in  = nan(Nac, 4);
p_out = nan(Nac, 4);

for a = 1:Nac
    [p_in_a, p_out_a] = contra_function_sector_intersection(AC_in(a).WP, sector_ab_k, a_band, sector_big, Vsec);
    p_in(a,:)  = p_in_a;
    p_out(a,:) = p_out_a;
end

end

%% ========================================================================
%  PASS-THROUGH FILTER
% ========================================================================

function index = contra_function_pass_through(AC, sector_big, Vsec, a_band)
Nac = numel(AC);

index = false(Nac,1);

for a = 1:Nac

    if ~isfield(AC(a),'WP') || isempty(AC(a).WP) || size(AC(a).WP,2) < 4
        continue
    end

    WP = AC(a).WP;
    traj = [WP(:,2), WP(:,3)];
    traj = local_drop_consecutive_duplicates(traj);

    if size(traj,1) < 2
        continue
    end

    % Determine start/end inside (horizontal, using union polygon)
    in0 = false; in1 = false;
    try
        in0 = inpolygon(WP(1,2),   WP(1,3),   Vsec(:,1), Vsec(:,2));
        in1 = inpolygon(WP(end,2), WP(end,3), Vsec(:,1), Vsec(:,2));
    catch
        in0 = false; in1 = false;
    end

    % Boundary crossing?
    crosses = false;
    try
        [xi, yi] = polyxpoly(traj(:,1), traj(:,2), Vsec(:,1), Vsec(:,2));
        crosses = ~isempty(xi);
    catch
        crosses = false;
    end

    % Vertical window overlap check
    max_alt = max(WP(:,4));
    min_alt = min(WP(:,4));
    ok_vertical = ~((max_alt <= (a_band(1) - 5)) || (min_alt > (a_band(end) + 5)));

    if ~ok_vertical
        continue
    end

    % Keep if it crosses the boundary, OR if it stays fully inside
    if crosses || (in0 && in1)
        index(a) = true;
    end
end

end

%% ========================================================================
%  ENTRY / EXIT COMPUTATION
% ========================================================================

function [p_in, p_out] = contra_function_sector_intersection(WP, sector_ab, a_band, sector_big, Vsec)
% Compute p_in/p_out. If no boundary crossing but inside, use endpoints.

p_in  = nan(1,4);
p_out = nan(1,4);

if isempty(WP) || size(WP,2) < 4
    return
end

traj = [WP(:,2), WP(:,3)];
traj = local_drop_consecutive_duplicates(traj);
if size(traj,1) < 2
    return
end

% Check boundary intersections with union polygon
xi = []; yi = [];
try
    [xi, yi] = polyxpoly(traj(:,1), traj(:,2), Vsec(:,1), Vsec(:,2));
catch
    xi = []; yi = [];
end

% If no boundary crossing: if fully inside, define entry/exit as endpoints
if isempty(xi)
    in0 = false; in1 = false;
    try
        in0 = inpolygon(WP(1,2),   WP(1,3),   Vsec(:,1), Vsec(:,2));
        in1 = inpolygon(WP(end,2), WP(end,3), Vsec(:,1), Vsec(:,2));
    catch
        in0 = false; in1 = false;
    end

    if in0 && in1
        p_in  = WP(1,1:4);
        p_out = WP(end,1:4);
    end
    return
end

% Otherwise: compute first/last intersection along the trajectory order
% We will project intersections onto cumulative arclength of trajectory polyline.
dist_vec = [0; cumsum(hypot(diff(traj(:,1)), diff(traj(:,2))))];
if dist_vec(end) <= 0
    return
end

% For each intersection point, estimate its distance-from-start by nearest segment
s_int = nan(numel(xi),1);
for k = 1:numel(xi)
    s_int(k) = local_project_point_to_polyline_arclen([xi(k), yi(k)], traj, dist_vec);
end

% Sort intersections along the path
[s_sorted, ord] = sort(s_int, 'ascend', 'MissingPlacement','last');
xi = xi(ord); yi = yi(ord);

% Keep only finite projected intersections
good = isfinite(s_sorted);
xi = xi(good); yi = yi(good); s_sorted = s_sorted(good);

if numel(s_sorted) < 1
    return
end

% First and last intersection points
x_in  = xi(1);  y_in  = yi(1);  s_in  = s_sorted(1);
x_out = xi(end);y_out = yi(end);s_out = s_sorted(end);

% Interpolate time and altitude versus arclength
t_int  = interp1(dist_vec, WP(:,1), [s_in; s_out], 'linear', 'extrap');
h_int  = interp1(dist_vec, WP(:,4), [s_in; s_out], 'linear', 'extrap');

p_in  = [t_int(1),  x_in,  y_in,  h_int(1)];
p_out = [t_int(2), x_out, y_out, h_int(2)];

end

%% ========================================================================
%  GEOMETRY HELPERS
% ========================================================================

function s = local_project_point_to_polyline_arclen(P, XY, dist_vec)
% Project point P onto polyline XY and return arclength coordinate (0..end).
% dist_vec is cumulative arclength for XY vertices.

x = XY(:,1); y = XY(:,2);
s = NaN;

best_d2 = inf;
best_s  = NaN;

for i = 1:(size(XY,1)-1)
    A = [x(i), y(i)];
    B = [x(i+1), y(i+1)];
    AB = B - A;
    AP = P - A;

    denom = AB*AB.';
    if denom <= eps
        t = 0;
    else
        t = (AP*AB.') / denom;
        t = max(0, min(1, t));
    end

    Q = A + t*AB;
    d2 = sum((P - Q).^2);

    if d2 < best_d2
        best_d2 = d2;
        seglen = hypot(AB(1), AB(2));
        best_s = dist_vec(i) + t*seglen;
    end
end

s = best_s;
end

function xy = local_drop_consecutive_duplicates(xy)
if isempty(xy) || size(xy,1) < 2
    return
end
d = diff(xy, 1, 1);
keep = [true; any(abs(d) > 0, 2)];
xy = xy(keep,:);
end
