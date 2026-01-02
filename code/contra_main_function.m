function [ASCR, sector_name] = contra_main_function()
% contra_main_function
% Computes ASCR for a single sector (ESMM31) and a single time window using:
%   - Trajectories from a CSV file (converted into AC(a).WP)
%   - Contrail cost obstacles derived from an ERA5 cost grid CSV (thresholded)
% Sector geometry is read from code_input/ESMM31.geojson.

%% ==================== USER SETTINGS ====================

sector_name = "ESMM31";
ASCR = NaN;  % safe default; will be overwritten if computation succeeds

thisFileDir = fileparts(mfilename('fullpath'));

sector_file = fullfile(thisFileDir, 'code_input', 'ESMM31.geojson');

traj_file = fullfile(thisFileDir, 'code_input', 'original_all_flights_in_ESMM31_day28.csv');
cost_file = fullfile(thisFileDir, 'code_input', 'grid_era5_smoothed_day28_ESMM31.csv');

TARGET_FL        = 370;
ALT_HALF_BAND_FT = 500;
MIN_FRACTION_IN_BAND = 0.80;

threshold  = 0.700;
conn       = 4;
round_prec = 6;

tz = "UTC";

%% ==================== Time interval (single window) ====================

t_ini_str = '2024-01-28 18:00:00';
t_fin_str = '2024-01-28 19:00:00';

t_vec_ini = datenum(t_ini_str, 'yyyy-mm-dd HH:MM:SS');
t_vec_fin = datenum(t_fin_str, 'yyyy-mm-dd HH:MM:SS');

t0_dt = datetime(t_ini_str, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', tz);
t1_dt = datetime(t_fin_str, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', tz);

%% ==================== Derived settings ====================

target_ft   = TARGET_FL * 100;
alt_low_ft  = target_ft - ALT_HALF_BAND_FT;
alt_high_ft = target_ft + ALT_HALF_BAND_FT;

COST_ALT_FT = target_ft;

toNum = @(x) str2double(strrep(string(x), ",", "."));

%% ========================= Read ESMM31 sector ===========================

if exist(sector_file, 'file') ~= 2
    error("Sector file not found: %s", sector_file);
end

G = jsondecode(fileread(sector_file));
if ~isfield(G, "features") || isempty(G.features)
    error("Invalid GeoJSON: missing features.");
end

feat0 = G.features(1);

if ~isfield(feat0, "geometry") || ~isfield(feat0.geometry, "coordinates")
    error("Invalid GeoJSON: missing geometry.coordinates.");
end

coords = feat0.geometry.coordinates;

if iscell(coords)
    ring = coords{1};
elseif isnumeric(coords) && ndims(coords) == 3
    ring = squeeze(coords(1,:,:));
elseif isnumeric(coords) && ndims(coords) == 2
    ring = coords;
else
    error("Unsupported geometry.coordinates format.");
end

lon = ring(:,1);
lat = ring(:,2);

if lon(1) ~= lon(end) || lat(1) ~= lat(end)
    lon(end+1) = lon(1);
    lat(end+1) = lat(1);
end

feat = struct();
feat.type = "Feature";

feat.properties = struct();
if isfield(feat0, "properties") && isfield(feat0.properties, "DESIGNATOR")
    feat.properties.DESIGNATOR = feat0.properties.DESIGNATOR;
else
    feat.properties.DESIGNATOR = "ESMM31";
end
if isfield(feat0, "properties") && isfield(feat0.properties, "LOWER_LIMIT_VALUE")
    feat.properties.LOWER_LIMIT_VALUE = feat0.properties.LOWER_LIMIT_VALUE;
else
    feat.properties.LOWER_LIMIT_VALUE = 0;
end
if isfield(feat0, "properties") && isfield(feat0.properties, "UPPER_LIMIT_VALUE")
    feat.properties.UPPER_LIMIT_VALUE = feat0.properties.UPPER_LIMIT_VALUE;
else
    feat.properties.UPPER_LIMIT_VALUE = 660;
end

feat.geometry = struct();
feat.geometry.type = "Polygon";
feat.geometry.coordinates = zeros(1, numel(lat), 2);
feat.geometry.coordinates(:,:,1) = lon(:).';
feat.geometry.coordinates(:,:,2) = lat(:).';

main_sectors = cell(1,1);
main_sectors{1} = feat;

sector_names = {'ESMM31'};
sector_time  = true;
sector_data  = main_sectors;

%% ================== Read TRAJECTORIES and build AC(a).WP ================

if exist(traj_file,'file') ~= 2
    error("Trajectory file not found: %s", traj_file);
end

fid = fopen(traj_file, 'r');
if fid < 0
    error("Could not open trajectory file: %s", traj_file);
end

header_line = fgetl(fid);
if ~ischar(header_line) || isempty(strtrim(header_line))
    fclose(fid);
    error("Could not read header line from trajectory file: %s", traj_file);
end

header_tokens = regexp(strtrim(header_line), '\s+', 'split');
ncols = numel(header_tokens);
fmt = repmat('%s', 1, ncols);

C = textscan(fid, fmt, ...
    'Delimiter', {' ','\t'}, ...
    'MultipleDelimsAsOne', true, ...
    'CollectOutput', true);
fclose(fid);

if isempty(C) || isempty(C{1})
    error("No trajectory data rows were read from: %s", traj_file);
end

Ttr = cell2table(C{1}, 'VariableNames', matlab.lang.makeUniqueStrings(header_tokens));

if ~all(ismember(["id","timestamp","latitude","longitude","altitude"], string(Ttr.Properties.VariableNames)))
    error("Trajectory file must contain columns: id, timestamp, latitude, longitude, altitude");
end

flight_id_all = string(Ttr.id);
epoch_ts_all  = toNum(Ttr.timestamp);
lat_all       = toNum(Ttr.latitude);
lon_all       = toNum(Ttr.longitude);
alt_ft_all    = toNum(Ttr.altitude);

good = ~ismissing(flight_id_all) & ~isnan(epoch_ts_all) & ...
       ~isnan(lat_all) & ~isnan(lon_all) & ~isnan(alt_ft_all);

flight_id_all = flight_id_all(good);
epoch_ts_all  = epoch_ts_all(good);
lat_all       = lat_all(good);
lon_all       = lon_all(good);
alt_ft_all    = alt_ft_all(good);

in_band_alt_all = (alt_ft_all >= alt_low_ft) & (alt_ft_all <= alt_high_ft);

uIDs = unique(flight_id_all, 'stable');
keep_flight = false(size(uIDs));
for i = 1:numel(uIDs)
    idx = (flight_id_all == uIDs(i));
    n = nnz(idx);
    if n == 0
        continue
    end
    keep_flight(i) = (nnz(in_band_alt_all(idx)) / n) >= MIN_FRACTION_IN_BAND;
end
keptIDs = uIDs(keep_flight);

if isempty(keptIDs)
    error("No trajectories satisfy the altitude band criterion.");
end

% Keep ALL points of kept flights (full trajectories)
keep_rows_tr = ismember(flight_id_all, keptIDs);

flight_id = flight_id_all(keep_rows_tr);
epoch_ts  = epoch_ts_all(keep_rows_tr);
lat_tr    = lat_all(keep_rows_tr);
lon_tr    = lon_all(keep_rows_tr);
alt_ft_tr = alt_ft_all(keep_rows_tr);

ts_tr = datetime(epoch_ts, 'ConvertFrom', 'posixtime', 'TimeZone', tz);

% Convert to AC(a).WP = [t, lon, lat, alt]
% Only keep points that fall within the requested time window (robustness: flights may have no points inside)
in_window = (ts_tr >= t0_dt) & (ts_tr < t1_dt);

if ~any(in_window)
    error("No trajectory points fall inside the requested time window %s to %s (UTC).", t_ini_str, t_fin_str);
end

flight_id_w = flight_id(in_window);
epoch_ts_w  = epoch_ts(in_window);
lon_w       = lon_tr(in_window);
lat_w       = lat_tr(in_window);
alt_w       = alt_ft_tr(in_window);

uWin = unique(flight_id_w, 'stable');
AC = struct('id', cell(numel(uWin),1), 'WP', cell(numel(uWin),1));

for i = 1:numel(uWin)
    id0 = uWin(i);
    idx = (flight_id_w == id0);
    if nnz(idx) < 2
        continue
    end
    t0 = epoch_ts_w(idx);
    [t0s, ord] = sort(t0);
    AC(i).id = id0;
    AC(i).WP = [t0s, lon_w(idx)];
    % Fix construction (expand columns cleanly)
    lon_i = lon_w(idx); lon_i = lon_i(ord);
    lat_i = lat_w(idx); lat_i = lat_i(ord);
    alt_i = alt_w(idx); alt_i = alt_i(ord);
    AC(i).WP = [t0s, lon_i, lat_i, alt_i];
end

% Drop empties (flights with <2 points in window)
emptyAC = arrayfun(@(s) isempty(s.WP) || size(s.WP,1) < 2, AC);
AC(emptyAC) = [];

if isempty(AC)
    error("After filtering to the time window, no flights had >= 2 points to form a trajectory.");
end

%% ========================= Read COST GRID ===============================

if exist(cost_file, 'file') ~= 2
    error("Cost CSV file not found: %s", cost_file);
end

opts = detectImportOptions(cost_file, "Delimiter", ",");
opts.Delimiter = ",";
opts.VariableNamesLine = 1;
opts.VariableNamingRule = "preserve";

iTs = find(strcmpi(opts.VariableNames, 'timestamp'), 1);
if isempty(iTs)
    error("No variable named 'timestamp' found in cost file.");
end
opts.VariableTypes{iTs} = 'char';

Tc = readtable(cost_file, opts);

requiredCols = ["latitude","longitude","cost","altitude","timestamp"];
missingCols = setdiff(requiredCols, string(Tc.Properties.VariableNames));
if ~isempty(missingCols)
    error("Missing required columns in cost CSV: %s", strjoin(missingCols, ", "));
end

Tc.latitude  = round(toNum(Tc.latitude),  round_prec);
Tc.longitude = round(toNum(Tc.longitude), round_prec);
Tc.cost      = toNum(Tc.cost);
Tc.altitude  = toNum(Tc.altitude);
tsKey_cost   = string(Tc.timestamp);

validc = ~ismissing(tsKey_cost) & ~isnan(Tc.latitude) & ~isnan(Tc.longitude) & ...
         ~isnan(Tc.cost) & ~isnan(Tc.altitude);
Tc = Tc(validc,:);
tsKey_cost = tsKey_cost(validc);

cost_alt_mask_all = (Tc.altitude == COST_ALT_FT);

%% ============= Create cost obstacles for the selected hour ==============

key0 = local_cost_key_utc(t0_dt);
in_hour_cost = (tsKey_cost == key0) & cost_alt_mask_all;

polygons = polyshape.empty(0,1);

if any(in_hour_cost)
    Th = Tc(in_hour_cost, :);

    allLatitudes  = sort(unique(Th.latitude));
    allLongitudes = sort(unique(Th.longitude));

    if ~isempty(allLatitudes) && ~isempty(allLongitudes)
        [tfLat, latIdx] = ismember(Th.latitude,  allLatitudes);
        [tfLon, lonIdx] = ismember(Th.longitude, allLongitudes);
        ok = tfLat & tfLon;

        Th_ok = Th(ok,:);
        latIdx_ok = latIdx(ok);
        lonIdx_ok = lonIdx(ok);

        occ = (Th_ok.cost > threshold);

        cost_matrix = false(numel(allLatitudes), numel(allLongitudes));
        if any(occ)
            linIdx = sub2ind(size(cost_matrix), latIdx_ok(occ), lonIdx_ok(occ));
            cost_matrix(linIdx) = true;
        end

        CC = bwconncomp(cost_matrix, conn);

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

            raw_coords = [allLongitudes(cb), allLatitudes(rb)];

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

            polygons(end+1,1) = pg;
        end
    end
end

cost_obstacles = {};
if ~isempty(polygons)
    cost_obstacles = cell(numel(polygons), 1);
    for i = 1:numel(polygons)
        cost_obstacles{i} = struct('pgon', polygons(i));
    end
end

%% ========================= ASCR computation =============================

adjacent_sectors = contra_function_create_adjacent_sectors();

ASCR = cell(1, 1);
ASCR{1} = nan(1, 1);

[sector_ab, a_band, flows_j] = contra_function_flows_sector_k( ...
    1, main_sectors, adjacent_sectors);

[p_in, p_out, AC_in] = contra_function_p_in_out( ...
    AC, sector_ab, a_band);

Wij_prev = contra_function_Wij_ini(flows_j);

[~, Wij, total_ac] = contra_function_Wj( ...
    t_vec_ini, t_vec_fin, p_in, p_out, a_band, flows_j, AC_in);

if total_ac == 0
    Wij = Wij_prev;
end

ASCR{1}(1, 1) = contra_function_ASCR_k( ...
    sector_ab, flows_j, cost_obstacles, Wij, a_band, adjacent_sectors);

end

%% ========================= Local functions ==============================

function key = local_cost_key_utc(tUTC)
try
    tUTC.TimeZone = "UTC";
catch
end
key = string(datestr(tUTC, 'yyyy-mm-dd HH:MM:SS')) + "+00:00";
end
