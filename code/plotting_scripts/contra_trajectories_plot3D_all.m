clear; clc;

%% -------------------------- Paths / settings ----------------------------
traj_file = fullfile('.', 'code_input', 'original_all_flights_in_ESMM31_day28.csv');

% UNIX epoch seconds -> interpret as UTC
tz = "UTC";

z_exaggeration = 10;   % keep the same overall vertical scale
view_elev = 22;
view_azim = -55;

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
    'CollectOutput', true, ...
    'ReturnOnError', false);
fclose(fid);

data = C{1};
if isempty(data)
    error("No data rows were read from: %s", traj_file);
end

T = cell2table(data, 'VariableNames', matlab.lang.makeUniqueStrings(header_tokens));
vars = string(T.Properties.VariableNames);

req = ["id","timestamp","latitude","longitude","altitude"];
missing = req(~ismember(req, vars));
if ~isempty(missing)
    error("Missing required columns: %s\nFound variables:\n%s", ...
        strjoin(missing, ", "), strjoin(vars, ", "));
end

%% ------------------------ Extract / convert -----------------------------
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

% Flight level (FL) from feet
alt_fl = alt_ft / 100;

% scale for vizualization:
fl_scale = (100 * 0.0003048) * z_exaggeration;
z_plot = alt_fl * fl_scale;

fprintf("Altitude range (ft):  %.0f .. %.0f\n", min(alt_ft), max(alt_ft));
fprintf("Altitude range (FL):  %.1f .. %.1f\n", min(alt_fl), max(alt_fl));
fprintf("Z range (scaled):     %.3f .. %.3f\n", min(z_plot), max(z_plot));

%% ---------------------- 3D plot (ALL trajectories) ----------------------
figure;
ax = axes; hold(ax, 'on');
grid(ax, 'on');

uIDs = unique(flight_id, 'stable');

for i = 1:numel(uIDs)
    idx = find(flight_id == uIDs(i));
    if numel(idx) < 2
        continue
    end

    [~, ord] = sort(epoch_ts(idx));
    idx = idx(ord);

    plot3(ax, lon(idx), lat(idx), z_plot(idx), 'LineWidth', 1);
end

xlabel(ax, 'Longitude');
ylabel(ax, 'Latitude');
zlabel(ax, sprintf('Flight level (FL) (scaled ×%.5g)', fl_scale));
title(ax, 'All trajectories (3D)');

% Camera/view
view(ax, view_azim, view_elev);

axis(ax, 'vis3d');
rotate3d(ax, 'on');

fl_ticks = [370 380 390 400 410];
set(ax, 'ZTick', fl_ticks * fl_scale);
set(ax, 'ZTickLabel', string(fl_ticks));

fig = gcf;
savefig(fig, 'all_trajectories_3D.fig');
