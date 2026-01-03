clear; clc;

%% ---------------------------- SETTINGS ----------------------------------
traj_file = fullfile('.', 'code_input', 'original_all_flights_in_ESMM31_day28.csv');

TARGET_FL = 390;

% Band around TARGET_FL used for:
% 1) deciding if a flight is "mostly at TARGET_FL"
% 2) selecting which points to plot
FL_HALF_BAND = 5;              % ±5 FL => ±500 ft (FL375..FL385 for TARGET_FL=380)

MIN_FRACTION_IN_BAND = 0.80;   % keep flight if >=80% of its points are inside the band

% Output folders
fig_root_base = fullfile('.', 'figures');
fig_root_traj = fullfile(fig_root_base, 'Trajectories');
fig_root      = fullfile(fig_root_traj, sprintf('FL%d', TARGET_FL));

% Plot export settings
dpi = 200;
tz  = "UTC";

%% ------------------------- Derived settings -----------------------------
FL_target_ft     = TARGET_FL * 100;
FL_halfwidth_ft  = FL_HALF_BAND * 100;

alt_cut_low_ft   = FL_target_ft - FL_halfwidth_ft;
alt_cut_high_ft  = FL_target_ft + FL_halfwidth_ft;

%% -------------------- Create output folders if needed -------------------
if ~exist(fig_root_base, 'dir'); mkdir(fig_root_base); end
if ~exist(fig_root_traj, 'dir'); mkdir(fig_root_traj); end
if ~exist(fig_root, 'dir');      mkdir(fig_root);      end

%% -------- Read trajectories (whitespace-delimited with header) ----------
if ~exist(traj_file,'file')
    error("File not found: %s", traj_file);
end

fid = fopen(traj_file, 'r');
if fid < 0
    error("Could not open file: %s", traj_file);
end

header_line = fgetl(fid);
if ~ischar(header_line) || isempty(strtrim(header_line))
    fclose(fid);
    error("Could not read header line from file: %s", traj_file);
end

header_tokens = regexp(strtrim(header_line), '\s+', 'split');
ncols = numel(header_tokens);

fmt = repmat('%s', 1, ncols);
C = textscan(fid, fmt, ...
    'Delimiter', {' ','\t'}, ...
    'MultipleDelimsAsOne', true, ...
    'CollectOutput', true);
fclose(fid);

data = C{1};
if isempty(data)
    error("No data rows were read from: %s", traj_file);
end

T = cell2table(data, 'VariableNames', matlab.lang.makeUniqueStrings(header_tokens));

%% ------------------------- Extract / convert ----------------------------
toNum = @(x) str2double(strrep(string(x), ",", "."));

flight_id = string(T.id);
epoch_ts  = toNum(T.timestamp);
lat       = toNum(T.latitude);
lon       = toNum(T.longitude);
alt_ft    = toNum(T.altitude);

good = ~ismissing(flight_id) & ~isnan(epoch_ts) & ~isnan(lat) & ~isnan(lon) & ~isnan(alt_ft);

flight_id = flight_id(good);
epoch_ts  = epoch_ts(good);
lat       = lat(good);
lon       = lon(good);
alt_ft    = alt_ft(good);

fprintf("Total rows read: %d\n", numel(lon));
fprintf("Unique flights:  %d\n", numel(unique(flight_id)));

%% ---------- Decide which flights to keep (mostly TARGET_FL) -------------
uIDs = unique(flight_id, 'stable');
in_band = alt_ft >= alt_cut_low_ft & alt_ft <= alt_cut_high_ft;

keep_flight = false(size(uIDs));
for i = 1:numel(uIDs)
    idx = (flight_id == uIDs(i));
    n = nnz(idx);
    if n == 0, continue; end
    keep_flight(i) = (nnz(in_band(idx)) / n) >= MIN_FRACTION_IN_BAND;
end

keptIDs = uIDs(keep_flight);

fprintf("Target band:     FL%d ± %d (ft: [%g..%g])\n", TARGET_FL, FL_HALF_BAND, alt_cut_low_ft, alt_cut_high_ft);
fprintf("Flights kept:    %d\n", numel(keptIDs));
fprintf("Flights dropped: %d\n", numel(uIDs) - numel(keptIDs));

if isempty(keptIDs)
    error("No flights satisfy the criterion. Lower MIN_FRACTION_IN_BAND or widen FL_HALF_BAND.");
end

%% -------------------- Keep all points of kept flights -------------------
keep_rows = ismember(flight_id, keptIDs);

flight_id = flight_id(keep_rows);
epoch_ts  = epoch_ts(keep_rows);
lat       = lat(keep_rows);
lon       = lon(keep_rows);
alt_ft    = alt_ft(keep_rows);

%% -------------------- Plot only points inside the band ------------------
keep_alt = alt_ft >= alt_cut_low_ft & alt_ft <= alt_cut_high_ft;

flight_id = flight_id(keep_alt);
epoch_ts  = epoch_ts(keep_alt);
lat       = lat(keep_alt);
lon       = lon(keep_alt);
alt_ft    = alt_ft(keep_alt);

fprintf("Rows after altitude band filter: %d\n", numel(lon));

if isempty(epoch_ts)
    error("No data left after filtering.");
end

%% -------------- Build UTC datetimes and determine the day ---------------
ts = datetime(epoch_ts, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
day0 = dateshift(min(ts), 'start', 'day');

% Consistent axes limits across all 24 plots
xlim_all = [min(lon), max(lon)];
ylim_all = [min(lat), max(lat)];

%% -------------- Create 24 hourly plots (saved, not shown) ---------------
for h = 0:23
    t0 = day0 + hours(h);
    t1 = day0 + hours(h+1);

    in_hour = (ts >= t0) & (ts < t1);
    if ~any(in_hour)
        flights_in_hour = string.empty(0,1);
    else
        flights_in_hour = unique(flight_id(in_hour), 'stable');
    end

    fig = figure('Visible','off');
    ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');

    for i = 1:numel(flights_in_hour)
        idx = in_hour & (flight_id == flights_in_hour(i));
        if nnz(idx) < 2
            continue
        end

        ii = find(idx);
        [~, ord] = sort(epoch_ts(ii));
        ii = ii(ord);

        plot(ax, lon(ii), lat(ii), 'LineWidth', 1);
    end

    xlabel(ax, 'Longitude');
    ylabel(ax, 'Latitude');
    axis(ax, 'equal');
    xlim(ax, xlim_all);
    ylim(ax, ylim_all);

    title(ax, sprintf('Trajectories %s (UTC) %s–%s', ...
        sprintf('FL%d±%d', TARGET_FL, FL_HALF_BAND), datestr(t0,'HH:MM'), datestr(t1,'HH:MM')), ...
        'Interpreter','none');

    fname = sprintf('trajectories_FL%d_%s_%02d00_%02d00.png', ...
        TARGET_FL, datestr(day0,'yyyymmdd'), h, mod(h+1,24));
    fpath = fullfile(fig_root, fname);

    exportgraphics(fig, fpath, 'Resolution', dpi);
    close(fig);
end

fprintf("Done. Saved 24 hourly trajectory plots to: %s\n", fig_root);
