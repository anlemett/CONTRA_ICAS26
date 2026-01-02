function [p_in, p_out, AC_in] = contra_function_p_in_out(AC, sector_ab_k, a_band)
% Computes entry (p_in) and exit (p_out) points of aircraft trajectories AC(a).WP
% through a sector defined by polygons sector_ab_k over altitude bands a_band.
%
% AC(a).WP is assumed to be [t, lon, lat, alt_ft_or_FL_units_consistent_with_a_band].
%
% Notes / robustness:
% - Removes consecutive duplicate horizontal points before geometric intersection.
% - Guards against empty/degenerate unions and intersections.

Nac = numel(AC);

sector_big = union([sector_ab_k{:}]);
in_sector  = false(Nac,1);

for a = 1:Nac

    if ~isfield(AC(a), 'WP') || isempty(AC(a).WP) || size(AC(a).WP,2) < 4
        continue
    end

    WP = AC(a).WP;

    % Horizontal check
    lineseg = [WP(:,2), WP(:,3)];
    lineseg = local_drop_consecutive_duplicates(lineseg);

    ok_horizontal = false;
    if size(lineseg,1) >= 2 && ~isempty(sector_big) && sector_big.NumRegions > 0
        try
            in = intersect(sector_big, lineseg);
            ok_horizontal = ~isempty(in);
        catch
            ok_horizontal = false;
        end
    end

    % Vertical check
    max_alt = max(WP(:,4));
    min_alt = min(WP(:,4));
    ok_vertical = ~((max_alt <= (a_band(1) - 5)) || (min_alt > (a_band(end) + 5)));

    if ok_horizontal && ok_vertical
        in_sector(a) = true;
    end
end

AC_in = AC(in_sector);

index_pass = contra_function_pass_through(AC_in, sector_ab_k, a_band);
AC_in = AC_in(index_pass);

Nac = numel(AC_in);

p_in  = nan(Nac, 4);
p_out = nan(Nac, 4);

for a = 1:Nac
    [p_in_a, p_out_a] = contra_function_sector_intersection(AC_in(a).WP, sector_ab_k, a_band);
    p_in(a,:)  = p_in_a;
    p_out(a,:) = p_out_a;
end

end

function [sector_div, a_band_div] = contra_function_divide_sector(sector_ab, a_band)
nab = size(a_band,1);

sector_div = cell(size(sector_ab));
a_band_div = zeros(nab,2);

sector_div{1}   = sector_ab{1};
a_band_div(1,1) = 1;

i_aux  = 1;
i_aux2 = 1;

for i = 2:nab
    polyout = subtract(sector_ab{i_aux2}, sector_ab{i});
    if polyout.NumRegions > 0
        a_band_div(i_aux,2) = i;
        i_aux = i_aux + 1;
        sector_div{i_aux} = sector_ab{i};
        a_band_div(i_aux,1) = i;
        i_aux2 = i;
    end
end

a_band_div(cellfun(@isempty,sector_div),:) = [];
a_band_div(end,2) = nab;
sector_div(cellfun(@isempty,sector_div)) = [];

index = a_band_div;
a_band_div(:,1) = a_band(index(:,1));
a_band_div(:,2) = a_band(index(:,2));

end

function index = contra_function_pass_through(AC, sector_ab_k, a_band)
Nac = numel(AC);

index_first = false(Nac,1);
index_last  = false(Nac,1);

for a = 1:Nac

    if ~isfield(AC(a),'WP') || isempty(AC(a).WP) || size(AC(a).WP,2) < 4
        continue
    end

    % First waypoint
    WP0 = [AC(a).WP(1,2), AC(a).WP(1,3), AC(a).WP(1,4)];

    if (WP0(3) < a_band(1)) || (WP0(3) > a_band(end))
        index_first(a) = true;
    else
        sec_cell = sector_ab_k((a_band >= (WP0(3)-5)) & (a_band < (WP0(3)+5)));
        if isempty(sec_cell) || isempty(sec_cell{1}) || sec_cell{1}.NumRegions == 0
            index_first(a) = true;
        else
            sec_poly = sec_cell{1};
            in_idx = inpolygon(WP0(1), WP0(2), sec_poly.Vertices(:,1), sec_poly.Vertices(:,2));
            index_first(a) = ~in_idx;
        end
    end

    % Last waypoint
    WP1 = [AC(a).WP(end,2), AC(a).WP(end,3), AC(a).WP(end,4)];

    if (WP1(3) < a_band(1)) || (WP1(3) > a_band(end))
        index_last(a) = true;
    else
        sec_cell = sector_ab_k((a_band >= (WP1(3)-5)) & (a_band < (WP1(3)+5)));
        if isempty(sec_cell) || isempty(sec_cell{1}) || sec_cell{1}.NumRegions == 0
            index_last(a) = true;
        else
            sec_poly = sec_cell{1};
            in_idx = inpolygon(WP1(1), WP1(2), sec_poly.Vertices(:,1), sec_poly.Vertices(:,2));
            index_last(a) = ~in_idx;
        end
    end
end

index = index_first & index_last;

end

function [p_in, p_out] = contra_function_sector_intersection(WP, sector_ab, a_band)
p_in  = nan(1,4);
p_out = nan(1,4);

p_in_found  = false;
p_out_found = false;

segments = contra_function_divide_trajectory(WP);
n_seg = numel(segments);

k = 1;
while (k <= n_seg) && ~(p_in_found && p_out_found)

    if ~p_in_found
        [p_in_found, p_inaux, p_out_found, p_outaux] = contra_function_find_p_in( ...
            segments(k).WP, segments(k).phase, sector_ab, a_band);

        if p_in_found
            p_in = p_inaux;

            if p_out_found
                p_out = p_outaux;
            else
                [p_out_found, p_outaux] = contra_function_find_p_out( ...
                    segments(k).WP, segments(k).phase, sector_ab, a_band);
                if p_out_found
                    p_out = p_outaux;
                else
                    k = k + 1;
                end
            end
        else
            k = k + 1;
        end
    else
        [p_out_found, p_outaux] = contra_function_find_p_out( ...
            segments(k).WP, segments(k).phase, sector_ab, a_band);
        if p_out_found
            p_out = p_outaux;
        else
            k = k + 1;
        end
    end
end

end

function [p_in_found, p_in, p_out_found, p_out] = contra_function_find_p_in(WP, phase, sector_ab, a_band)
p_in_found  = false;
p_out_found = false;
p_in  = [];
p_out = [];

lineseg = [WP(:,2), WP(:,3)];
lineseg = local_drop_consecutive_duplicates(lineseg);

switch phase

    case 0  % level

        h = WP(1,4);
        ab_idx = ((a_band-5) <= h) & ((a_band+5) > h);

        if any(ab_idx) && size(lineseg,1) >= 2
            sec_poly = sector_ab{ab_idx};

            in = [];
            try
                in = intersect(sec_poly, lineseg);
            catch
                in = [];
            end

            if ~isempty(in)
                p_in_found = true;

                in_end = inpolygon(WP(end,2), WP(end,3), sec_poly.Vertices(:,1), sec_poly.Vertices(:,2));

                lon_int = in([1,end],1);
                lat_int = in([1,end],2);

                dist_vec = sqrt((WP(:,2)-WP(1,2)).^2 + (WP(:,3)-WP(1,3)).^2);
                dist_int = sqrt((lon_int-WP(1,2)).^2 + (lat_int-WP(1,3)).^2);
                t_int = interp1(dist_vec, WP(:,1), dist_int);

                [time, idxT] = sort(t_int);
                longitude    = lon_int(idxT);
                latitude     = lat_int(idxT);

                p_in = [time(1), longitude(1), latitude(1), h];

                if ~in_end
                    p_out_found = true;
                    p_out = [time(2), longitude(2), latitude(2), h];
                end
            end
        end

    case 1  % climb

        h_ini = WP(1,4);
        h_end = WP(end,4);

        if (h_end >= a_band(1)) && (h_ini <= a_band(end)) && size(lineseg,1) >= 2

            [sector_div, a_band_div] = contra_function_divide_sector(sector_ab, a_band);
            n_div = numel(sector_div);

            for kk = 1:n_div
                if h_ini < a_band_div(kk,2)
                    sec_poly = sector_div{kk};

                    in = [];
                    try
                        in = intersect(sec_poly, lineseg);
                    catch
                        in = [];
                    end

                    if ~isempty(in)

                        lon_int = in([1,end],1);
                        lat_int = in([1,end],2);

                        dist_vec = sqrt((WP(:,2)-WP(1,2)).^2 + (WP(:,3)-WP(1,3)).^2);
                        dist_int = sqrt((lon_int-WP(1,2)).^2 + (lat_int-WP(1,3)).^2);
                        t_int = interp1(dist_vec, WP(:,1), dist_int);
                        h_int = interp1(dist_vec, WP(:,4), dist_int);

                        [time, idxT] = sort(t_int);
                        altitude     = h_int(idxT);
                        longitude    = lon_int(idxT);
                        latitude     = lat_int(idxT);

                        in_ini = inpolygon(WP(1,2), WP(1,3), sec_poly.Vertices(:,1), sec_poly.Vertices(:,2));

                        if in_ini
                            if WP(1,4) < (a_band_div(kk,1)-5)
                                h_bottom = a_band_div(kk,1)-5;
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_bottom);
                                time_int = int_values(1); lon_i = int_values(2); lat_i = int_values(3);
                                if inpolygon(lon_i, lat_i, sec_poly.Vertices(:,1), sec_poly.Vertices(:,2))
                                    p_in_found = true;
                                    p_in = [time_int, lon_i, lat_i, h_bottom];
                                    break
                                end
                            end
                        else
                            if (altitude(1) >= (a_band_div(kk,1)-5)) && (altitude(1) < (a_band_div(kk,2)+5))
                                p_in_found = true;
                                p_in = [time(1), longitude(1), latitude(1), altitude(1)];
                                break
                            elseif altitude(1) < (a_band_div(kk,1)-5)
                                h_bottom = a_band_div(kk,1)-5;
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_bottom);
                                time_int = int_values(1); lon_i = int_values(2); lat_i = int_values(3);
                                if inpolygon(lon_i, lat_i, sec_poly.Vertices(:,1), sec_poly.Vertices(:,2))
                                    p_in_found = true;
                                    p_in = [time_int, lon_i, lat_i, h_bottom];
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end

    case -1 % descent

        h_ini = WP(1,4);
        h_end = WP(end,4);

        if ((h_end <= a_band(end)) || (h_ini >= a_band(1))) && size(lineseg,1) >= 2

            [sector_div, a_band_div] = contra_function_divide_sector(sector_ab, a_band);
            n_div = numel(sector_div);

            for kk = n_div:-1:1
                if h_ini > a_band_div(kk,1)
                    sec_poly = sector_div{kk};

                    in = [];
                    try
                        in = intersect(sec_poly, lineseg);
                    catch
                        in = [];
                    end

                    if ~isempty(in)

                        lon_int = in([1,end],1);
                        lat_int = in([1,end],2);

                        dist_vec = sqrt((WP(:,2)-WP(1,2)).^2 + (WP(:,3)-WP(1,3)).^2);
                        dist_int = sqrt((lon_int-WP(1,2)).^2 + (lat_int-WP(1,3)).^2);
                        t_int = interp1(dist_vec, WP(:,1), dist_int);
                        h_int = interp1(dist_vec, WP(:,4), dist_int);

                        [time, idxT] = sort(t_int);
                        altitude     = h_int(idxT);
                        longitude    = lon_int(idxT);
                        latitude     = lat_int(idxT);

                        in_ini = inpolygon(WP(1,2), WP(1,3), sec_poly.Vertices(:,1), sec_poly.Vertices(:,2));

                        if in_ini
                            if WP(1,4) > (a_band_div(kk,2)+5)
                                h_top = a_band_div(kk,2)+5;
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_top);
                                time_int = int_values(1); lon_i = int_values(2); lat_i = int_values(3);
                                if inpolygon(lon_i, lat_i, sec_poly.Vertices(:,1), sec_poly.Vertices(:,2))
                                    p_in_found = true;
                                    p_in = [time_int, lon_i, lat_i, h_top];
                                    break
                                end
                            end
                        else
                            if (altitude(1) >= (a_band_div(kk,1)-5)) && (altitude(1) < (a_band_div(kk,2)+5))
                                p_in_found = true;
                                p_in = [time(1), longitude(1), latitude(1), altitude(1)];
                                break
                            elseif altitude(1) > (a_band_div(kk,2)+5)
                                h_top = a_band_div(kk,1)-5; % kept as in legacy code
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_top);
                                time_int = int_values(1); lon_i = int_values(2); lat_i = int_values(3);
                                if inpolygon(lon_i, lat_i, sec_poly.Vertices(:,1), sec_poly.Vertices(:,2))
                                    p_in_found = true;
                                    p_in = [time_int, lon_i, lat_i, h_top];
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
end

end

function [p_out_found, p_out] = contra_function_find_p_out(WP, phase, sector_ab, a_band)
p_out_found = false;
p_out = [];

lineseg = [WP(:,2), WP(:,3)];
lineseg = local_drop_consecutive_duplicates(lineseg);

switch phase

    case 0  % level

        h = WP(1,4);
        ab_idx = ((a_band-5) <= h) & ((a_band+5) > h);

        if any(ab_idx) && size(lineseg,1) >= 2
            sec_poly = sector_ab{ab_idx};

            in = [];
            try
                in = intersect(sec_poly, lineseg);
            catch
                in = [];
            end

            if ~isempty(in)
                in_end = inpolygon(WP(end,2), WP(end,3), sec_poly.Vertices(:,1), sec_poly.Vertices(:,2));
                if ~in_end
                    p_out_found = true;

                    lon_int = in([1,end],1);
                    lat_int = in([1,end],2);

                    dist_vec = sqrt((WP(:,2)-WP(1,2)).^2 + (WP(:,3)-WP(1,3)).^2);
                    dist_int = sqrt((lon_int-WP(1,2)).^2 + (lat_int-WP(1,3)).^2);
                    t_int = interp1(dist_vec, WP(:,1), dist_int);

                    [time, idxT] = sort(t_int);
                    longitude    = lon_int(idxT);
                    latitude     = lat_int(idxT);

                    p_out = [time(2), longitude(2), latitude(2), h];
                end
            end
        end

    case 1  % climb

        h_ini = WP(1,4);
        h_end = WP(end,4);

        if (h_end >= a_band(1)) && (h_ini <= a_band(end)) && size(lineseg,1) >= 2

            [sector_div, a_band_div] = contra_function_divide_sector(sector_ab, a_band);
            n_div = numel(sector_div);

            for kk = n_div:-1:1
                if h_end > a_band_div(kk,1)
                    sec_poly = sector_div{kk};

                    in = [];
                    try
                        in = intersect(sec_poly, lineseg);
                    catch
                        in = [];
                    end

                    if ~isempty(in)

                        lon_int = in([1,end],1);
                        lat_int = in([1,end],2);

                        dist_vec = sqrt((WP(:,2)-WP(1,2)).^2 + (WP(:,3)-WP(1,3)).^2);
                        dist_int = sqrt((lon_int-WP(1,2)).^2 + (lat_int-WP(1,3)).^2);
                        t_int = interp1(dist_vec, WP(:,1), dist_int);
                        h_int = interp1(dist_vec, WP(:,4), dist_int);

                        [time, idxT] = sort(t_int);
                        altitude     = h_int(idxT);
                        longitude    = lon_int(idxT);
                        latitude     = lat_int(idxT);

                        in_end = inpolygon(WP(end,2), WP(end,3), sec_poly.Vertices(:,1), sec_poly.Vertices(:,2));

                        if in_end
                            if WP(end,4) > (a_band_div(kk,2)+5)
                                h_top = a_band_div(kk,2)+5;
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_top);
                                time_int = int_values(1); lon_i = int_values(2); lat_i = int_values(3);
                                if inpolygon(lon_i, lat_i, sec_poly.Vertices(:,1), sec_poly.Vertices(:,2))
                                    p_out_found = true;
                                    p_out = [time_int, lon_i, lat_i, h_top];
                                    break
                                end
                            end
                        else
                            if (altitude(2) >= (a_band_div(kk,1)-5)) && (altitude(2) < (a_band_div(kk,2)+5))
                                p_out_found = true;
                                p_out = [time(2), longitude(2), latitude(2), altitude(2)];
                                break
                            elseif altitude(2) > (a_band_div(kk,2)+5)
                                h_top = a_band_div(kk,1)-5;
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_top);
                                time_int = int_values(1); lon_i = int_values(2); lat_i = int_values(3);
                                if inpolygon(lon_i, lat_i, sec_poly.Vertices(:,1), sec_poly.Vertices(:,2))
                                    p_out_found = true;
                                    p_out = [time_int, lon_i, lat_i, h_top];
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end

    case -1  % descent

        h_ini = WP(1,4);
        h_end = WP(end,4);

        if ((h_end <= a_band(end)) || (h_ini >= a_band(1))) && size(lineseg,1) >= 2

            [sector_div, a_band_div] = contra_function_divide_sector(sector_ab, a_band);
            n_div = numel(sector_div);

            for kk = 1:n_div
                if h_end < a_band_div(kk,2)
                    sec_poly = sector_div{kk};

                    in = [];
                    try
                        in = intersect(sec_poly, lineseg);
                    catch
                        in = [];
                    end

                    if ~isempty(in)

                        lon_int = in([1,end],1);
                        lat_int = in([1,end],2);

                        dist_vec = sqrt((WP(:,2)-WP(1,2)).^2 + (WP(:,3)-WP(1,3)).^2);
                        dist_int = sqrt((lon_int-WP(1,2)).^2 + (lat_int-WP(1,3)).^2);
                        t_int = interp1(dist_vec, WP(:,1), dist_int);
                        h_int = interp1(dist_vec, WP(:,4), dist_int);

                        [time, idxT] = sort(t_int);
                        altitude     = h_int(idxT);
                        longitude    = lon_int(idxT);
                        latitude     = lat_int(idxT);

                        in_end = inpolygon(WP(end,2), WP(end,3), sec_poly.Vertices(:,1), sec_poly.Vertices(:,2));

                        if in_end
                            if WP(end,4) < (a_band_div(kk,1)-5)
                                h_bottom = a_band_div(kk,1)-5;
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_bottom);
                                time_int = int_values(1); lon_i = int_values(2); lat_i = int_values(3);
                                if inpolygon(lon_i, lat_i, sec_poly.Vertices(:,1), sec_poly.Vertices(:,2))
                                    p_out_found = true;
                                    p_out = [time_int, lon_i, lat_i, h_bottom];
                                    break
                                end
                            end
                        else
                            if (altitude(2) >= (a_band_div(kk,1)-5)) && (altitude(2) < (a_band_div(kk,2)+5))
                                p_out_found = true;
                                p_out = [time(2), longitude(2), latitude(2), altitude(2)];
                                break
                            elseif altitude(2) < (a_band_div(kk,1)-5)
                                h_bottom = a_band_div(kk,1)-5;
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_bottom);
                                time_int = int_values(1); lon_i = int_values(2); lat_i = int_values(3);
                                if inpolygon(lon_i, lat_i, sec_poly.Vertices(:,1), sec_poly.Vertices(:,2))
                                    p_out_found = true;
                                    p_out = [time_int, lon_i, lat_i, h_bottom];
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
end

end

function segments = contra_function_divide_trajectory(WP)
[nwp,~] = size(WP);

alt_dif = diff(WP(:,4));
if isempty(alt_dif)
    segments = struct('phase', 0, 'WP', WP);
    return
end

% Robust segmentation: treat zeros as previous non-zero sign if possible
sgn = sign(alt_dif);
for i = 2:numel(sgn)
    if sgn(i) == 0
        sgn(i) = sgn(i-1);
    end
end
if sgn(1) == 0
    nz = find(sgn ~= 0, 1);
    if ~isempty(nz)
        sgn(1) = sgn(nz);
    else
        sgn(:) = 0;
    end
end

breaks = find(diff(sgn) ~= 0);
nseg = numel(breaks) + 1;

segments(nseg) = struct('phase', [], 'WP', []);

i_vec = 1;
for i = 1:nseg
    type = sgn(max(1, min(i_vec, numel(sgn))));
    segments(i).phase = type;

    i_ini = i_vec;
    while i_vec <= numel(sgn) && sgn(i_vec) == type
        i_vec = i_vec + 1;
        if i_vec == nwp
            break
        end
    end

    segments(i).WP = WP(i_ini:i_vec, :);
end

end

function xy = local_drop_consecutive_duplicates(xy)
if isempty(xy) || size(xy,1) < 2
    return
end
d = diff(xy, 1, 1);
keep = [true; any(abs(d) > 0, 2)];
xy = xy(keep,:);
end
