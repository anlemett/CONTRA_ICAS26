function contra_plot_dummy_adj_sectors()
clc;

%% ==================== Path (relative to this script) ====================
thisFileDir = fileparts(mfilename('fullpath'));
sector_file = fullfile(thisFileDir, 'code_input', 'ESMM31.geojson');

if exist(sector_file, 'file') ~= 2
    error('Sector file not found: %s', sector_file);
end

fprintf('Using sector file:\n  %s\n', sector_file);

%% ======================== Read ESMM31 GeoJSON ===========================
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

if lon(1) ~= lon(end) || lat(1) ~= lat(end)
    lon(end+1) = lon(1);
    lat(end+1) = lat(1);
end

%% ================ Dummy adjacent sectors (border segments) ==============
[~, iS] = min(lat);
[~, iN] = max(lat);
[~, iW] = min(lon);
[~, iE] = max(lon);

n = numel(lon);

idx_SW = local_wrap_path(iW, iS, n);  % West -> South
idx_NW = local_wrap_path(iW, iN, n);  % West -> North
idx_NE = local_wrap_path(iN, iE, n);  % North -> East
idx_SE = local_wrap_path(iE, iS, n);  % East -> South

cLat = mean(lat);
cLon = mean(lon);

padLat = 0.60;
padLon = 0.80;

outer_S = [cLon, min(lat) - padLat];
outer_N = [cLon, max(lat) + padLat];
outer_W = [min(lon) - padLon,  cLat];
outer_E = [max(lon) + padLon,  cLat];

adj(1).name = "DUMMY_SW";
adj(1).lon  = [lon(idx_SW); outer_S(1); outer_W(1); lon(idx_SW(1))];
adj(1).lat  = [lat(idx_SW); outer_S(2); outer_W(2); lat(idx_SW(1))];

adj(2).name = "DUMMY_NW";
adj(2).lon  = [lon(idx_NW); outer_N(1); outer_W(1); lon(idx_NW(1))];
adj(2).lat  = [lat(idx_NW); outer_N(2); outer_W(2); lat(idx_NW(1))];

adj(3).name = "DUMMY_NE";
adj(3).lon  = [lon(idx_NE); outer_E(1); outer_N(1); lon(idx_NE(1))];
adj(3).lat  = [lat(idx_NE); outer_E(2); outer_N(2); lat(idx_NE(1))];

adj(4).name = "DUMMY_SE";
adj(4).lon  = [lon(idx_SE); outer_S(1); outer_E(1); lon(idx_SE(1))];
adj(4).lat  = [lat(idx_SE); outer_S(2); outer_E(2); lat(idx_SE(1))];

%% ============================ Plot ======================================
figure('Color','w'); hold on; grid on;

plot(lon, lat, 'k-', 'LineWidth', 2);
text(mean(lon,'omitnan'), mean(lat,'omitnan'), "ESMM31", ...
    'FontWeight','bold', 'HorizontalAlignment','center');

for i = 1:4
    pg = polyshape(adj(i).lon, adj(i).lat, 'Simplify', true);
    plot(pg, 'FaceAlpha', 0.20, 'LineWidth', 1.5);
    text(mean(adj(i).lon,'omitnan'), mean(adj(i).lat,'omitnan'), adj(i).name, ...
        'HorizontalAlignment','center');
end

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('ESMM31 (from ESMM31.geojson) and Dummy Adjacent Sectors');
axis equal;

end

%% ========================= Local function ===============================

function idx = local_wrap_path(iStart, iEnd, n)
% Returns indices along the boundary going forward (with wrap) from iStart to iEnd.
if iStart <= iEnd
    idx = (iStart:iEnd).';
else
    idx = [ (iStart:n) 1:iEnd ].';
end
end
