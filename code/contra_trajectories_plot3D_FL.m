clear; clc;

%% -------------------------- Paths / settings ----------------------------
traj_file = fullfile('.', 'code_input', 'original_all_flights_in_ESMM31_day28.csv');

z_exaggeration = 10;
view_elev = 22;
view_azim = -55;

%% --------------- Keep flights that are (mostly) on one FL ---------------
FL_target_ft        = 37000;
FL_halfwidth_ft     = 500;
min_fraction_in_FL  = 0.80;    % keep flight if >=80% of points in the band

%% ---------------------- Simple altitude cuts ----------------------------
alt_cut_low_ft  = FL_target_ft - FL_halfwidth_ft; % drop trajectory parts below this altitude
alt_cut_high_ft = FL_target_ft + FL_halfwidth_ft; % drop trajectory parts above this altitude

% Z-axis tick labels in FL
FL_ticks = 360:5:380;

%% ---------------------- Read trajectories -------------------------------
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

good = ~ismissing(flight_id) & ~isnan(epoch_ts) & ...
       ~isnan(lat) & ~isnan(lon) & ~isnan(alt_ft);

flight_id = flight_id(good);
epoch_ts  = epoch_ts(good);
lat       = lat(good);
lon       = lon(good);
alt_ft    = alt_ft(good);

fprintf("Total rows read: %d\n", numel(lon));
fprintf("Unique flights:  %d\n", numel(unique(flight_id)));

%% --------------- Decide which flights to keep (mostly FL) ---------------
uIDs = unique(flight_id, 'stable');
is_in_FL = alt_ft >= (FL_target_ft - FL_halfwidth_ft) & ...
           alt_ft <= (FL_target_ft + FL_halfwidth_ft);

keep_flight = false(size(uIDs));
for i = 1:numel(uIDs)
    idx = (flight_id == uIDs(i));
    n = nnz(idx);
    if n == 0, continue; end
    keep_flight(i) = (nnz(is_in_FL(idx)) / n) >= min_fraction_in_FL;
end

keptIDs = uIDs(keep_flight);

fprintf("Flights kept (mostly FL370): %d\n", numel(keptIDs));
fprintf("Flights dropped:            %d\n", numel(uIDs) - numel(keptIDs));

if isempty(keptIDs)
    error("No flights satisfy the FL criterion. Lower min_fraction_in_FL (e.g., 0.5).");
end

%% -------------------- Keep ALL points of kept flights -------------------
keep_rows = ismember(flight_id, keptIDs);

flight_id = flight_id(keep_rows);
epoch_ts  = epoch_ts(keep_rows);
lat       = lat(keep_rows);
lon       = lon(keep_rows);
alt_ft    = alt_ft(keep_rows);

%% ------------------------ Apply altitude window -------------------------
keep_alt = (alt_ft >= alt_cut_low_ft) & (alt_ft <= alt_cut_high_ft);

flight_id = flight_id(keep_alt);
epoch_ts  = epoch_ts(keep_alt);
lat       = lat(keep_alt);
lon       = lon(keep_alt);
alt_ft    = alt_ft(keep_alt);

fprintf("Rows after altitude window [%g .. %g] ft: %d\n", ...
    alt_cut_low_ft, alt_cut_high_ft, numel(lon));

alt_km = alt_ft * 0.0003048;
z_plot = alt_km * z_exaggeration;

%% ---------------------------- 3D plot -----------------------------------
figure;
ax = axes; hold(ax, 'on');
grid(ax, 'on');

uIDs2 = unique(flight_id, 'stable');

for i = 1:numel(uIDs2)
    idx = find(flight_id == uIDs2(i));
    if numel(idx) < 2
        continue
    end

    [~, ord] = sort(epoch_ts(idx));
    idx = idx(ord);

    plot3(ax, lon(idx), lat(idx), z_plot(idx), 'LineWidth', 1);
end

xlabel(ax, 'Longitude');
ylabel(ax, 'Latitude');
zlabel(ax, 'Flight level (FL)');

z_ticks = (FL_ticks * 100) * 0.0003048 * z_exaggeration;
set(ax, 'ZTick', z_ticks);
set(ax, 'ZTickLabel', compose('%d', FL_ticks));

title(ax, sprintf('Flights around FL370 (FL365–FL375)'), 'Interpreter','none');

view(ax, view_azim, view_elev);
axis(ax, 'vis3d');
rotate3d(ax, 'on');

savefig(gcf, 'trajectories_FL370_window_365_375.fig');
