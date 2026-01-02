function [sectors_pgon, adjacent_pgon] = contra_function_sector_adjacent_pgon_a_band(sector_data, adjacent_sectors_data, altitude_band)
% Build 2D polyshape polygons for main sectors and adjacent sectors at a given altitude band (FL).

% Main sectors
N = numel(sector_data);
sectors_pgon = cell(size(sector_data));

for i = 1:N
    sector_i = sector_data{i};
    if isempty(sector_i)
        sectors_pgon{i} = polyshape.empty(0,1);
        continue
    end

    % Filter subsectors by altitude band
    [sector_i, Ni] = local_filter_by_altitude(sector_i, altitude_band);

    if Ni <= 0
        sectors_pgon{i} = polyshape.empty(0,1);
        continue
    end

    % Union all subsector polygons
    pgon = polyshape.empty(0,1);
    for j = 1:Ni
        [lon, lat] = local_get_lonlat(sector_i(j));
        if numel(lon) < 3 || numel(lat) < 3
            continue
        end

        % Ensure clockwise orientation
        if ~ispolycw(lon, lat)
            [lon, lat] = poly2cw(lon, lat);
        end

        try
            %pj = polyshape(lon, lat, 'Simplify', true, 'KeepCollinearPoints', false);
            pj = polyshape(lon, lat, 'Simplify', false, 'KeepCollinearPoints', true);
        catch
            continue
        end

        if pj.NumRegions == 0
            continue
        end

        if isempty(pgon)
            pgon = pj;
        else
            pgon = union(pgon, pj);
        end
    end

    sectors_pgon{i} = pgon;
end

% Adjacent sectors
M = numel(adjacent_sectors_data);
adjacent_pgon = cell(size(adjacent_sectors_data));

for i = 1:M
    sector_i = adjacent_sectors_data{i};
    if isempty(sector_i)
        adjacent_pgon{i} = polyshape.empty(0,1);
        continue
    end

    [sector_i, Ni] = local_filter_by_altitude(sector_i, altitude_band);

    if Ni <= 0
        adjacent_pgon{i} = polyshape.empty(0,1);
        continue
    end

    pgon = polyshape.empty(0,1);
    for j = 1:Ni
        [lon, lat] = local_get_lonlat(sector_i(j));
        if numel(lon) < 3 || numel(lat) < 3
            continue
        end

        if ~ispolycw(lon, lat)
            [lon, lat] = poly2cw(lon, lat);
        end

        try
            %pj = polyshape(lon, lat, 'Simplify', true, 'KeepCollinearPoints', false);
            pj = polyshape(lon, lat, 'Simplify', false, 'KeepCollinearPoints', true);
        catch
            continue
        end

        if pj.NumRegions == 0
            continue
        end

        if isempty(pgon)
            pgon = pj;
        else
            pgon = union(pgon, pj);
        end
    end

    adjacent_pgon{i} = pgon;
end

end

%% ==================== Local helpers ====================

function [sector_f, Ni] = local_filter_by_altitude(sector_i, altitude_band)
% Return subsectors of sector_i that include altitude_band (FL).
% Supports:
%   - sector_i as struct array with .properties fields
%   - sector_i as scalar struct (single subsector)

if ~isstruct(sector_i)
    sector_f = [];
    Ni = 0;
    return
end

% If sector_i is scalar, treat as 1 subsector
if numel(sector_i) == 1
    lower_FL = sector_i.properties.LOWER_LIMIT_VALUE;
    upper_FL = sector_i.properties.UPPER_LIMIT_VALUE;
    keep = (lower_FL < altitude_band) && (upper_FL > altitude_band);
    if keep
        sector_f = sector_i;
        Ni = 1;
    else
        sector_f = [];
        Ni = 0;
    end
    return
end

aux = [sector_i.properties];
upper_FL = [aux.UPPER_LIMIT_VALUE]';
lower_FL = [aux.LOWER_LIMIT_VALUE]';

A = altitude_band * ones(size(lower_FL));
sec_index = (lower_FL < A) & (upper_FL > A);

sector_f = sector_i(sec_index);
Ni = numel(sector_f);

end

function [lon, lat] = local_get_lonlat(sec)
% Extract lon/lat vectors

coords = sec.geometry.coordinates;

% Expected shapes:coords is (1 x N x 2) or (N x 1 x 2) or numeric (N x 2)
if iscell(coords)
    % If coordinates are stored as cell, take first
    coords = coords{1};
end

if isnumeric(coords) && ndims(coords) == 3
    lon = squeeze(coords(:,:,1));
    lat = squeeze(coords(:,:,2));
    lon = lon(:);
    lat = lat(:);
elseif isnumeric(coords) && size(coords,2) == 2
    lon = coords(:,1);
    lat = coords(:,2);
else
    lon = [];
    lat = [];
end

end
