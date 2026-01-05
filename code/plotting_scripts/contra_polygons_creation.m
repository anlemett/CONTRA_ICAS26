clear; clc;

%% ------------------------ Paths / settings ------------------------------
csv_file = fullfile('.', 'code_input', 'grid_era5_smoothed_day28_ESMM31.csv');

threshold  = 0.7;

conn       = 4;      % 4- or 8-connectivity
round_prec = 6;      % rounding precision for lat/lon indexing
dpi        = 200;    % PNG resolution

% Root folders
fig_root_base  = fullfile('.', 'figures');
json_root_base = fullfile('.', 'jsons');

% Add threshold-level subfolder
thr_token   = regexprep(sprintf('thr_%0.3f', threshold), '\.', 'p');
fig_root    = fullfile(fig_root_base,  thr_token);
json_root   = fullfile(json_root_base, thr_token);

%% --------------------------- Helpers ------------------------------------
toNum = @(x) str2double(strrep(string(x), ",", "."));

% Safe filename token: keep [A-Za-z0-9_-], replace everything else with "_"
sanitizeToken = @(s) regexprep(string(s), '[^A-Za-z0-9_\-]', '_');

% Altitude token: avoid '.' and '-' in folder/file names
altToken = @(a) regexprep(regexprep(sprintf('alt_%g', a), '\.', 'p'), '-', 'm');

%% ------------------------ Sector boundary  ------------------------------
raw_data = [ ...
"575020N 0141939E 565727N 0155231E" newline ...
"565200N 0153959E 564101N 0151432E" newline ...
"563741N 0150704E 563006N 0145321E" newline ...
"561900N 0143212E 555250N 0134214E" newline ...
"554600N 0132050E 553831N 0130536E" newline ...
"553233N 0125337E 553101N 0125032E" newline ...
"553238N 0124831E 560433N 0120806E" newline ...
"563124N 0115721E 563601N 0120410E" newline ...
"565522N 0123308E 570504N 0125135E" newline ...
"573059N 0134209E 575020N 0141939E" ];

% Convert DMS strings to decimal degrees
dms_to_dd = @(s) local_dms_to_dd(s);

lines = splitlines(strtrim(raw_data));
tokens = {};
for i = 1:numel(lines)
    if strlength(strtrim(lines(i))) == 0, continue; end
    tokens = [tokens; cellstr(split(strtrim(lines(i))))];
end
tokens = tokens(~cellfun(@isempty, tokens));

if mod(numel(tokens), 2) ~= 0
    error("Sector coordinate tokens count is not even; cannot form lat/lon pairs.");
end

sector_lat = zeros(numel(tokens)/2, 1);
sector_lon = zeros(numel(tokens)/2, 1);

k = 0;
for i = 1:2:numel(tokens)
    k = k + 1;
    lat_str = string(tokens{i});
    lon_str = string(tokens{i+1});
    sector_lat(k) = dms_to_dd(lat_str);
    sector_lon(k) = dms_to_dd(lon_str);
end

% Ensure closed ring for plotting
if sector_lat(1) ~= sector_lat(end) || sector_lon(1) ~= sector_lon(end)
    sector_lat(end+1) = sector_lat(1);
    sector_lon(end+1) = sector_lon(1);
end

%% ---------------- Create root output folders if needed ------------------
if ~exist(fig_root_base, 'dir');  mkdir(fig_root_base);  end
if ~exist(json_root_base,'dir');  mkdir(json_root_base); end
if ~exist(fig_root, 'dir');       mkdir(fig_root);       end
if ~exist(json_root,'dir');       mkdir(json_root);      end

%% ----------------------- Read table robustly ----------------------------
if ~exist(csv_file, 'file')
    error("CSV file not found: %s", csv_file);
end

opts = detectImportOptions(csv_file, "Delimiter", ",");
opts.Delimiter = ",";
opts.VariableNamesLine = 1;
opts.VariableNamingRule = "preserve";

% timestamp as text
iTs = find(strcmpi(opts.VariableNames, 'timestamp'), 1);
if isempty(iTs)
    error("No variable named 'timestamp' found in opts.VariableNames.");
end
opts.VariableTypes{iTs} = 'char';

T = readtable(csv_file, opts);

requiredCols = ["latitude","longitude","cost","altitude","timestamp"];
missingCols = setdiff(requiredCols, string(T.Properties.VariableNames));
if ~isempty(missingCols)
    error("Missing required columns in CSV: %s", strjoin(missingCols, ", "));
end

% Convert numeric columns
T.latitude  = toNum(T.latitude);
T.longitude = toNum(T.longitude);
T.cost      = toNum(T.cost);
T.altitude  = toNum(T.altitude);

% Round lat/lon once for stable grid indexing
T.latitude  = round(T.latitude,  round_prec);
T.longitude = round(T.longitude, round_prec);

% Timestamp key
tsKey = string(T.timestamp);

%% -------------- Determine unique timestamp/altitude values --------------
valid = ~ismissing(tsKey) & ~isnan(T.altitude) & ...
        ~isnan(T.latitude) & ~isnan(T.longitude) & ~isnan(T.cost);

tsVals  = unique(tsKey(valid), "stable");
altVals = unique(T.altitude(valid), "stable");

fprintf('Unique timestamps (valid): %d\n', numel(tsVals));
fprintf('Unique altitudes  (valid): %d\n', numel(altVals));

if isempty(tsVals) || isempty(altVals)
    error("No valid timestamp/altitude values after filtering.");
end

%% -------------- Loop over all (timestamp, altitude) pairs ---------------
% Suppress the specific polyshape "repaired" warning (still repairs, just silences)
warnState = warning('off', 'MATLAB:polyshape:repairedBySimplify');
cleanupWarn = onCleanup(@() warning(warnState));

totalPairs = numel(tsVals) * numel(altVals);
pairCount  = 0;
madeCount  = 0;

% accumulate polygon counts per altitude across timestamps
polygons_per_altitude = zeros(numel(altVals), 1);

for a = 1:numel(altVals)
    alt0 = altVals(a);
    alt_dirname = altToken(alt0);

    % Figures: threshold folder -> altitude folder
    fig_alt_dir = fullfile(fig_root, alt_dirname);
    if ~exist(fig_alt_dir, 'dir'); mkdir(fig_alt_dir); end

    % JSONs: threshold folder -> altitude folder
    json_alt_dir = fullfile(json_root, alt_dirname);
    if ~exist(json_alt_dir, 'dir'); mkdir(json_alt_dir); end

    for t = 1:numel(tsVals)
        ts0 = tsVals(t);
        pairCount = pairCount + 1;

        % Filter
        idx = valid & (T.altitude == alt0) & (tsKey == ts0);
        Tf = T(idx, :);
        if isempty(Tf)
            continue
        end

        % Build grid
        allLatitudes  = sort(unique(Tf.latitude));
        allLongitudes = sort(unique(Tf.longitude));
        if numel(allLatitudes) < 1 || numel(allLongitudes) < 1
            continue
        end

        % Grid spacing for plot bounds
        if numel(allLatitudes) >= 2,  dlat = median(diff(allLatitudes)); else, dlat = 0.01; end
        if numel(allLongitudes) >= 2, dlon = median(diff(allLongitudes)); else, dlon = 0.01; end

        [tfLat, latIdx] = ismember(Tf.latitude,  allLatitudes);
        [tfLon, lonIdx] = ismember(Tf.longitude, allLongitudes);
        ok = tfLat & tfLon;

        Tf_ok = Tf(ok,:);
        if isempty(Tf_ok)
            continue
        end

        occ = (Tf_ok.cost > threshold);

        cost_matrix = false(numel(allLatitudes), numel(allLongitudes));
        if any(occ)
            linIdx = sub2ind(size(cost_matrix), latIdx(ok & occ), lonIdx(ok & occ));
            cost_matrix(linIdx) = true;
        end

        % Connected components -> polygons (repair with polyshape; export repaired boundary)
        CC = bwconncomp(cost_matrix, conn);

        features = struct('type', {}, 'geometry', {}, 'properties', {});
        polygon_number = 0;

        for component = 1:CC.NumObjects
            idxList = CC.PixelIdxList{component};
            [r, c] = ind2sub(size(cost_matrix), idxList);

            if numel(r) < 3
                continue
            end

            try
                bnd = boundary(c, r, 1);
            catch
                continue
            end
            if isempty(bnd) || numel(bnd) < 3
                continue
            end

            rb = r(bnd);
            cb = c(bnd);

            raw_coords = [allLongitudes(cb), allLatitudes(rb)]; % [lon, lat]

            % Remove consecutive duplicates
            if size(raw_coords,1) >= 2
                keep = [true; any(diff(raw_coords,1,1) ~= 0, 2)];
                raw_coords = raw_coords(keep,:);
            end
            if size(raw_coords,1) < 3
                continue
            end

            try
                pg = polyshape(raw_coords(:,1), raw_coords(:,2), ...
                               'Simplify', true, 'KeepCollinearPoints', false);
            catch
                continue
            end

            if pg.NumRegions ~= 1
                continue
            end

            [vx, vy] = boundary(pg);
            good = ~(isnan(vx) | isnan(vy));
            vx = vx(good); vy = vy(good);

            if numel(vx) < 3
                continue
            end

            coords = [vx(:), vy(:)];

            % Close ring
            if ~isequal(coords(1,:), coords(end,:))
                coords(end+1,:) = coords(1,:);
            end
            if size(coords,1) < 4
                continue
            end

            polygon_number = polygon_number + 1;

            geometry_struct = struct( ...
                "type", "Polygon", ...
                "coordinates", { { coords } } ...
            );

            features(polygon_number) = struct( ...
                "type", "Feature", ...
                "geometry", geometry_struct, ...
                "properties", struct() ...
            );
        end

        % add polygons for this timestamp to altitude total
        polygons_per_altitude(a) = polygons_per_altitude(a) + polygon_number;

        % ----------------------- File names ------------------------------
        ts_token = sanitizeToken(ts0);

        % base includes threshold + altitude + timestamp
        base = sprintf('%s__%s__%s', thr_token, alt_dirname, "ts_" + ts_token);

        json_name = base + ".geojson";
        json_path = fullfile(json_alt_dir, json_name);

        png_name  = base + ".png";
        png_path  = fullfile(fig_alt_dir, png_name);

        % --------------------- Write GeoJSON -----------------------------
        fc_name = string(base);
        if polygon_number == 0
            json_structure = struct("type","FeatureCollection","name",fc_name,"features",[]);
        else
            json_structure = struct("type","FeatureCollection","name",fc_name,"features",features);
        end

        % encodedJSON = jsonencode(json_structure);
        % fid = fopen(json_path, 'w');
        % if fid < 0
        %     error("Could not open output file for writing: %s", json_path);
        % end
        % fprintf(fid, '%s', encodedJSON);
        % fclose(fid);

        % ----------------------- Create plot -----------------------------
        fig = figure('Visible','off'); hold on;

        for k = 1:polygon_number
            ring = features(k).geometry.coordinates{1}; % N x 2
            pgon = polyshape(ring(:,1), ring(:,2), 'Simplify', true, 'KeepCollinearPoints', false);
            plot(pgon);
        end

        plot(Tf_ok.longitude(occ), Tf_ok.latitude(occ), 'k.', 'MarkerSize', 10);

        % --------- Plot ATC sector boundary on top (lon,lat) -------------
        plot(sector_lon, sector_lat, 'LineWidth', 2);

        % -- Axis limits must include BOTH cost grid and sector boundary --
        x_min = min([allLongitudes(:); sector_lon(:)]);
        x_max = max([allLongitudes(:); sector_lon(:)]);
        y_min = min([allLatitudes(:);  sector_lat(:)]);
        y_max = max([allLatitudes(:);  sector_lat(:)]);

        xlim([x_min - dlon, x_max + dlon]);
        ylim([y_min - dlat, y_max + dlat]);

        % Geographic-looking proportions for lon/lat degrees
        % At latitude ~57N, 1 deg lon is shorter than 1 deg lat by cos(lat)
        meanLat = mean([Tf_ok.latitude; sector_lat], "omitnan");
        if ~isnan(meanLat)
            daspect([1/cosd(meanLat), 1, 1]);  % compress x relative to y
        end

        xlabel('Longitude');
        ylabel('Latitude');
        title(sprintf('thr=%0.3f | alt=%g | ts=%s', threshold, alt0, ts0), 'Interpreter', 'none');

        exportgraphics(fig, png_path, 'Resolution', dpi);
        close(fig);

        madeCount = madeCount + 1;

        if mod(pairCount, 50) == 0
            fprintf('Processed %d / %d pairs; created %d outputs.\n', pairCount, totalPairs, madeCount);
        end
    end

    fprintf('Altitude %g: total polygons across all timestamps = %d\n', alt0, polygons_per_altitude(a));
end

fprintf('Done. Created %d JSON+PNG pairs.\n', madeCount);

%% ------------------------ Local function --------------------------------
function dd = local_dms_to_dd(dms_string)
    s = char(dms_string);
    dirc = upper(s(end));
    dms_value = s(1:end-1);

    n = length(dms_value);
    if n == 6
        D = str2double(dms_value(1:2));
        M = str2double(dms_value(3:4));
        S = str2double(dms_value(5:6));
    elseif n == 7
        D = str2double(dms_value(1:3));
        M = str2double(dms_value(4:5));
        S = str2double(dms_value(6:7));
    else
        error("Invalid DMS string length: %s", s);
    end

    dd = D + M/60 + S/3600;
    if dirc == 'S' || dirc == 'W'
        dd = -dd;
    end
end
