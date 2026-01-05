clear; clc;

%% ============================ SETTINGS ==================================
% Files
traj_file = fullfile('.', 'code_input', 'original_all_flights_in_ESMM31_day28_extended.csv');

% Target flight level and band (feet)
TARGET_FL        = 390;
ALT_HALF_BAND_FT = 500;      % ±500 ft

% Time / output
tz  = "UTC";
dpi = 200;

fig_root_base = fullfile('.', 'figures');
fig_root_traj = fullfile(fig_root_base, 'Trajectories');
fig_root      = fullfile(fig_root_traj, sprintf('FL%d', TARGET_FL));

%% ======================== Derived settings ==============================
target_ft   = TARGET_FL * 100;
alt_low_ft  = target_ft - ALT_HALF_BAND_FT;
alt_high_ft = target_ft + ALT_HALF_BAND_FT;

%% ======================= Create output folders ==========================
if ~exist(fig_root_base, 'dir'); mkdir(fig_root_base); end
if ~exist(fig_root_traj, 'dir'); mkdir(fig_root_traj); end
if ~exist(fig_root, 'dir');      mkdir(fig_root);      end

%% ============================== Helpers =================================
toNum = @(x) str2double(strrep(string(x), ",", "."));

%% ========================== Sector boundary =============================
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
    sector_lat(k) = local_dms_to_dd(string(tokens{i}));
    sector_lon(k) = local_dms_to_dd(string(tokens{i+1}));
end

if sector_lat(1) ~= sector_lat(end) || sector_lon(1) ~= sector_lon(end)
    sector_lat(end+1) = sector_lat(1);
    sector_lon(end+1) = sector_lon(1);
end

%% ========== Read TRAJECTORIES (whitespace-delimited with header) ========
if ~exist(traj_file,'file')
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

C = textscan(fid, fmt, 'Delimiter', {' ','\t'}, 'MultipleDelimsAsOne', true, 'CollectOutput', true);
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

good = ~ismissing(flight_id_all) & ~isnan(epoch_ts_all) & ~isnan(lat_all) & ~isnan(lon_all) & ~isnan(alt_ft_all);

flight_id_all = flight_id_all(good);
epoch_ts_all  = epoch_ts_all(good);
lat_all       = lat_all(good);
lon_all       = lon_all(good);
alt_ft_all    = alt_ft_all(good);

fprintf("Traj: rows read (valid):    %d\n", numel(lon_all));
fprintf("Traj: unique flights:       %d\n", numel(unique(flight_id_all)));

% Datetime for all points
ts_all = datetime(epoch_ts_all, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
day0   = dateshift(min(ts_all), 'start', 'day');

in_band_alt_all = (alt_ft_all >= alt_low_ft) & (alt_ft_all <= alt_high_ft);
in_sector_all   = inpolygon(lon_all, lat_all, sector_lon, sector_lat);

uIDs = unique(flight_id_all, 'stable');
keep_flight = false(size(uIDs));

for i = 1:numel(uIDs)
    idx = (flight_id_all == uIDs(i));
    if nnz(idx) < 2
        continue
    end

    ii = find(idx);
    [~, ord] = sort(epoch_ts_all(ii));
    ii = ii(ord);

    tfi = ts_all(ii);
    inBandPts = in_band_alt_all(ii);
    inSecPts  = in_sector_all(ii);

    segInBand = inBandPts(1:end-1) & inBandPts(2:end);
    segInSec  = inSecPts(1:end-1)  & inSecPts(2:end);

    hasSolidBlueSomeHour = false;
    for h = 0:23
        t0 = day0 + hours(h);
        t1 = day0 + hours(h+1);
        inHrPts = (tfi >= t0) & (tfi < t1);
        segInHr = inHrPts(1:end-1) & inHrPts(2:end);

        if any(segInHr & segInBand & segInSec)
            hasSolidBlueSomeHour = true;
            break
        end
    end

    keep_flight(i) = hasSolidBlueSomeHour;
end

keptIDs = uIDs(keep_flight);

fprintf("Traj: altitude band:        FL%d ± %d ft => [%g..%g] ft\n", TARGET_FL, ALT_HALF_BAND_FT, alt_low_ft, alt_high_ft);
fprintf("Traj: kept flights:         %d\n", numel(keptIDs));
fprintf("Traj: dropped flights:      %d\n", numel(uIDs) - numel(keptIDs));

if isempty(keptIDs)
    error("No trajectories produce a SOLID BLUE inside-sector in-hour segment. Widen ALT_HALF_BAND_FT or adjust TARGET_FL.");
end

% Keep ALL points of kept flights (full trajectories)
keep_rows_tr = ismember(flight_id_all, keptIDs);

flight_id = flight_id_all(keep_rows_tr);
epoch_ts  = epoch_ts_all(keep_rows_tr);
lat_tr    = lat_all(keep_rows_tr);
lon_tr    = lon_all(keep_rows_tr);
alt_ft_tr = alt_ft_all(keep_rows_tr);
ts_tr     = ts_all(keep_rows_tr);

% Inside-sector mask for kept flights (point-wise)
in_sector_tr = inpolygon(lon_tr, lat_tr, sector_lon, sector_lat);

%% ===================== Global plot limits / aspect ======================
x_all = sector_lon(:);
y_all = sector_lat(:);

xlim_all = [min(x_all), max(x_all)];
ylim_all = [min(y_all), max(y_all)];

meanLatAll = mean(y_all, "omitnan");
if isnan(meanLatAll), meanLatAll = 57; end
asp = [1/cosd(meanLatAll), 1, 1];

%% ========================== Hourly plots ================================
% Only INSIDE-SECTOR segments are drawn
% Color encodes time:
%   - Blue  => within this hour
%   - Grey  => outside this hour
% Dash encodes altitude band:
%   - Solid => inside altitude band
%   - Dashed=> outside altitude band

for h = 0:23
    t0 = day0 + hours(h);
    t1 = day0 + hours(h+1);

    % count/plot flights that have at least one INSIDE-SECTOR
    % segment in this hour AND within the altitude band
    candidates = unique(flight_id((ts_tr >= t0) & (ts_tr < t1)), 'stable');
    keep = false(size(candidates));

    for k = 1:numel(candidates)
        id = candidates(k);

        idx = (flight_id == id);
        if nnz(idx) < 2
            continue
        end

        ii = find(idx);
        [~, ord] = sort(epoch_ts(ii));
        ii = ii(ord);

        pts = (ts_tr(ii) >= t0) & (ts_tr(ii) < t1) & in_sector_tr(ii) & (alt_ft_tr(ii) >= alt_low_ft) & (alt_ft_tr(ii) <= alt_high_ft);
        keep(k) = any(pts(1:end-1) & pts(2:end));
    end

    flights_in_hour = candidates(keep);

    fig = figure('Visible','off');
    ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');

    for i = 1:numel(flights_in_hour)
        id0 = flights_in_hour(i);

        idxAll = (flight_id == id0);
        if nnz(idxAll) < 2
            continue
        end

        iiAll = find(idxAll);
        [~, ordAll] = sort(epoch_ts(iiAll));
        iiAll = iiAll(ordAll);

        % Point-wise masks for this flight
        inHrPts   = (ts_tr(iiAll) >= t0) & (ts_tr(iiAll) < t1);
        inBandPts = (alt_ft_tr(iiAll) >= alt_low_ft) & (alt_ft_tr(iiAll) <= alt_high_ft);
        inSecPts  = in_sector_tr(iiAll);

        % Segment-wise masks
        seg_inHr   = inHrPts(1:end-1)   & inHrPts(2:end);
        seg_inBand = inBandPts(1:end-1) & inBandPts(2:end);
        seg_inSec  = inSecPts(1:end-1)  & inSecPts(2:end);

        % Only inside-sector segments
        seg = seg_inSec;

        local_plot_runs(ax, lon_tr, lat_tr, iiAll, seg &  seg_inHr &  seg_inBand, '-',  [0 0.4470 0.7410], 2.6); % blue solid
        local_plot_runs(ax, lon_tr, lat_tr, iiAll, seg &  seg_inHr & ~seg_inBand, '--', [0 0.4470 0.7410], 2.2); % blue dashed
        local_plot_runs(ax, lon_tr, lat_tr, iiAll, seg & ~seg_inHr &  seg_inBand, '-',  [0.55 0.55 0.55], 1.8);  % grey solid
        local_plot_runs(ax, lon_tr, lat_tr, iiAll, seg & ~seg_inHr & ~seg_inBand, '--', [0.55 0.55 0.55], 1.8);  % grey dashed
    end

    % Sector boundary
    plot(ax, sector_lon, sector_lat, 'LineWidth', 2);

    xlabel(ax, 'Longitude');
    ylabel(ax, 'Latitude');
    axis(ax, 'equal');
    xlim(ax, xlim_all);
    ylim(ax, ylim_all);
    daspect(ax, asp);

    title(ax, sprintf('In-sector only | FL%d ± %d ft | UTC %s–%s | flights=%d', ...
        TARGET_FL, ALT_HALF_BAND_FT, datestr(t0,'HH:MM'), datestr(t1,'HH:MM'), numel(flights_in_hour)), ...
        'Interpreter','none');

    fname = sprintf('traj_inSector_FL%d_%s_%02d00_%02d00.png', ...
        TARGET_FL, datestr(day0,'yyyymmdd'), h, mod(h+1,24));
    exportgraphics(fig, fullfile(fig_root, fname), 'Resolution', dpi);
    close(fig);
end

fprintf("Done. Saved 24 hourly plots to: %s\n", fig_root);

%% ========================== Local functions =============================

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

function local_plot_runs(ax, lon, lat, iiAll, segMask, lineStyle, col, lw)
    % Plot contiguous runs of segments.
    % segMask length = numel(iiAll)-1, describing segments between consecutive points.
    if isempty(segMask) || ~any(segMask)
        return
    end
    runs = local_true_runs(segMask);
    for r = 1:size(runs,1)
        s = runs(r,1);
        e = runs(r,2);
        pts = iiAll(s:(e+1));
        if numel(pts) >= 2
            plot(ax, lon(pts), lat(pts), lineStyle, 'Color', col, 'LineWidth', lw);
        end
    end
end

function runs = local_true_runs(tf)
    tf = tf(:);
    if isempty(tf) || ~any(tf)
        runs = zeros(0,2);
        return
    end
    d = diff([false; tf; false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    runs = [starts ends];
end
