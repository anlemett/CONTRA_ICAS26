clear; clc;

csv_file = fullfile('.', 'code_input', 'grid_era5_smoothed_day28_ESMM31.csv');

opts = detectImportOptions(csv_file, "Delimiter", ",");
opts.Delimiter = ",";
opts.VariableNamesLine = 1;
opts.VariableNamingRule = "preserve";

disp(opts.VariableNames);

% Exact match using char vectors
iTs = find(strcmpi(opts.VariableNames, 'timestamp'), 1);
if isempty(iTs)
    error("No variable named 'timestamp' found in opts.VariableNames.");
end

% Force timestamp to be read as text BEFORE readtable
opts.VariableTypes{iTs} = 'char';

T = readtable(csv_file, opts);
fprintf("%s\n", T.Properties.VariableNames{:});

disp(T(1,:))
tsKey  = string(T.timestamp); 

% Fix numeric columns explicitly (works for cell/string/categorical/numeric)
toNum = @(x) str2double(strrep(string(x), ",", "."));  % also handles decimal comma

T.latitude  = toNum(T.latitude);
T.longitude = toNum(T.longitude);
T.cost      = toNum(T.cost);
T.altitude  = toNum(T.altitude);


% Remove rows with invalid values
valid = ~ismissing(tsKey) & ~isnan(T.altitude) & ...
        ~isnan(T.latitude) & ~isnan(T.longitude) & ~isnan(T.cost);

% Select first valid timestamp and altitude (just an example)
tsVals  = unique(tsKey(valid), "stable");
altVals = unique(T.altitude(valid), "stable");

ts0  = tsVals(1);
alt0 = altVals(1);

% Filter
Tf = T(tsKey == ts0 & T.altitude == alt0, :);

allLatitudes = sort(unique(Tf.latitude));
allLongitudes = sort(unique(Tf.longitude));


tic

data_dir = fullfile('.', 'code_input');

 
json_name = 'ESMM31_rect_obstacles.geojson';
full_file_name = fullfile(data_dir, json_name);
        
fid=fopen(full_file_name,'w');
        

% rectangles around points where cost > thr

thr = 0.700;
mask = Tf.cost > thr;

% Grid spacing (degrees). Works for regular lat/lon grids.
dlat = median(diff(allLatitudes));
dlon = median(diff(allLongitudes));

% Build per-point rectangles in lon/lat around each selected grid point
sel = find(mask);
nSel = numel(sel);

features = struct('type', {}, 'geometry', {}, 'properties', {});
kFeat = 0;

for k = 1:nSel
    i = sel(k);

    lat = Tf.latitude(i);
    lon = Tf.longitude(i);
    cst = Tf.cost(i);

    % Rectangle corners (lon,lat). Close ring.
    coords = [ ...
        lon - dlon/2, lat - dlat/2; ...
        lon + dlon/2, lat - dlat/2; ...
        lon + dlon/2, lat + dlat/2; ...
        lon - dlon/2, lat + dlat/2; ...
        lon - dlon/2, lat - dlat/2  ...
    ];

    kFeat = kFeat + 1;
    geom = struct("type","Polygon", ...
                  "coordinates", reshape(coords, 1, size(coords,1), 2));

    features(kFeat) = struct( ...
        "type","Feature", ...
        "geometry", geom, ...
        "properties", struct("cost", cst, "thr", thr) ...
    );
end

% FeatureCollection (one rectangle per point)
name = "cost_cells";
geo = struct("type","FeatureCollection","name",name,"features",features);

% Write GeoJSON
data_dir = fullfile('.', 'code_input');
if ~exist(data_dir,'dir'), mkdir(data_dir); end
full_file_name = fullfile(data_dir, 'cost_cells.geojson');
fid = fopen(full_file_name,'w');
fprintf(fid, "%s", jsonencode(geo));
fclose(fid);

% Visualization: rectangles + underlying points
figure; hold on;
rectColor = [1 0 0];
for k = 1:kFeat
    coord = features(k).geometry.coordinates;  % 1 x N x 2
    x = squeeze(coord(1,:,1));
    y = squeeze(coord(1,:,2));
    pg = polyshape(x, y);
    plot(pg, ...
        'FaceColor', rectColor, ...
        'FaceAlpha', 0.5, ...     % transparency (0–1)
        'EdgeColor', rectColor, ...
        'LineWidth', 1);
end
plot(Tf.longitude(mask), Tf.latitude(mask), 'k.', 'MarkerSize', 12);
xlim([min(allLongitudes)-dlon, max(allLongitudes)+dlon]);
ylim([min(allLatitudes)-dlat,  max(allLatitudes)+dlat]);
xlabel('Longitude'); ylabel('Latitude');
title(sprintf('Rectangles around points with cost > %.4g', thr));
