function adjacent_sectors_data = contra_function_create_adjacent_sectors()
% contra_function_create_adjacent_sectors
% Creates 4 dummy adjacent sectors around ESMM31 using ESMM31 border segments
% between extreme points, derived from code_input/ESMM31.geojson.
%
% Adjacent sectors:
%   1) DUMMY_SW: ESMM31 border from West -> South + 2 outer points
%   2) DUMMY_NW: ESMM31 border from West -> North + 2 outer points
%   3) DUMMY_NE: ESMM31 border from North -> East + 2 outer points
%   4) DUMMY_SE: ESMM31 border from East -> South + 2 outer points
%
% Output format:
% adjacent_sectors_data is a 1x4 cell.
% Each cell contains a 1x1 struct with fields:
%   properties.DESIGNATOR, LOWER_LIMIT_VALUE, UPPER_LIMIT_VALUE
%   geometry.coordinates (1 x P x 2) with [lon, lat]

%% ==================== Load ESMM31 boundary from GeoJSON ====================

thisFileDir = fileparts(mfilename('fullpath'));
sector_file = fullfile(thisFileDir, 'code_input', 'ESMM31.geojson');

if exist(sector_file, 'file') ~= 2
    error('Sector file not found: %s', sector_file);
end

raw = fileread(sector_file);
G = jsondecode(raw);

if ~isfield(G, 'features') || isempty(G.features)
    error('Invalid GeoJSON: missing features.');
end

feat = G.features(1);

if ~isfield(feat, 'geometry') || ~isfield(feat.geometry, 'coordinates')
    error('Invalid GeoJSON: missing geometry.coordinates.');
end

coords = feat.geometry.coordinates;

if iscell(coords)
    ring = coords{1};                 % exterior ring [P x 2]
elseif isnumeric(coords) && ndims(coords) == 3
    ring = squeeze(coords(1,:,:));    % [P x 2]
elseif isnumeric(coords) && ndims(coords) == 2
    ring = coords;                    % [P x 2]
else
    error('Unsupported geometry.coordinates format.');
end

lon = ring(:,1);
lat = ring(:,2);

% Ensure closed polygon
if lon(1) ~= lon(end) || lat(1) ~= lat(end)
    lon(end+1) = lon(1);
    lat(end+1) = lat(1);
end

n = numel(lon);

%% ==================== Extreme points ====================

[~, iS] = min(lat);
[~, iN] = max(lat);
[~, iW] = min(lon);
[~, iE] = max(lon);

idx_SW = local_wrap_path(iW, iS, n);  % West -> South
idx_NW = local_wrap_path(iW, iN, n);  % West -> North
idx_NE = local_wrap_path(iN, iE, n);  % North -> East
idx_SE = local_wrap_path(iE, iS, n);  % East -> South

%% ==================== Outer points (extend outward) ====================

cLat = mean(lat);
cLon = mean(lon);

padLat = 0.60;
padLon = 0.80;

outer_S = [cLon, min(lat) - padLat];
outer_N = [cLon, max(lat) + padLat];
outer_W = [min(lon) - padLon,  cLat];
outer_E = [max(lon) + padLon,  cLat];

%% ============================ Build 4 adjacent sectors =============================

adjacent_sectors_data = cell(1,4);

adjacent_sectors_data{1} = local_make_sector( ...
    'DUMMY_SW', 0, 660, [lon(idx_SW); outer_S(1); outer_W(1); lon(idx_SW(1))], ...
                        [lat(idx_SW); outer_S(2); outer_W(2); lat(idx_SW(1))]);

adjacent_sectors_data{2} = local_make_sector( ...
    'DUMMY_NW', 0, 660, [lon(idx_NW); outer_N(1); outer_W(1); lon(idx_NW(1))], ...
                        [lat(idx_NW); outer_N(2); outer_W(2); lat(idx_NW(1))]);

adjacent_sectors_data{3} = local_make_sector( ...
    'DUMMY_NE', 0, 660, [lon(idx_NE); outer_E(1); outer_N(1); lon(idx_NE(1))], ...
                        [lat(idx_NE); outer_E(2); outer_N(2); lat(idx_NE(1))]);

adjacent_sectors_data{4} = local_make_sector( ...
    'DUMMY_SE', 0, 660, [lon(idx_SE); outer_S(1); outer_E(1); lon(idx_SE(1))], ...
                        [lat(idx_SE); outer_S(2); outer_E(2); lat(idx_SE(1))]);

end

%% ==================== Local helpers ====================

function idx = local_wrap_path(iStart, iEnd, n)
% Returns indices along boundary going forward (with wrap) from iStart to iEnd.
if iStart <= iEnd
    idx = (iStart:iEnd).';
else
    idx = [ (iStart:n) 1:iEnd ].';
end
end

function sec = local_make_sector(des, fl_lo, fl_hi, lon_vec, lat_vec)
% Ensure closed ring
if lon_vec(1) ~= lon_vec(end) || lat_vec(1) ~= lat_vec(end)
    lon_vec(end+1) = lon_vec(1);
    lat_vec(end+1) = lat_vec(1);
end

sec = struct('properties', [], 'geometry', []);
sec.properties.DESIGNATOR = des;
sec.properties.LOWER_LIMIT_VALUE = fl_lo;
sec.properties.UPPER_LIMIT_VALUE = fl_hi;

sec.geometry.coordinates = zeros(1, numel(lat_vec), 2);
sec.geometry.coordinates(:,:,1) = lon_vec(:).';
sec.geometry.coordinates(:,:,2) = lat_vec(:).';
end
