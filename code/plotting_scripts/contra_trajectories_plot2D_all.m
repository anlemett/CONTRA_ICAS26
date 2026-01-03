clear; clc;

%% ============================ SETTINGS ==================================
traj_file = fullfile('.', 'code_input', 'original_all_flights_in_ESMM31_day28_extended.csv');

% Output folders
fig_root_base = fullfile('.', 'figures');
fig_root_traj = fullfile(fig_root_base, 'Trajectories_All_2D');
fig_root      = fig_root_traj;

% Plot export settings
dpi = 200;

%% ======================= Create output folders ==========================
if ~exist(fig_root_base, 'dir'); mkdir(fig_root_base); end
if ~exist(fig_root_traj, 'dir'); mkdir(fig_root_traj); end
if ~exist(fig_root, 'dir');      mkdir(fig_root);      end

%% ============================ Helpers ===================================
toNum = @(x) str2double(strrep(string(x), ",", "."));

%% ========================= Sector boundary ==============================
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

% Close ring
if sector_lat(1) ~= sector_lat(end) || sector_lon(1) ~= sector_lon(end)
    sector_lat(end+1) = sector_lat(1);
    sector_lon(end+1) = sector_lon(1);
end

%% ======== Read TRAJECTORIES (whitespace-delimited with header) ==========
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

T = cell2table(C{1}, 'VariableNames', matlab.lang.makeUniqueStrings(header_tokens));

requiredCols = ["id","timestamp","latitude","longitude"];
missingCols = setdiff(requiredCols, string(T.Properties.VariableNames));
if ~isempty(missingCols)
    error("Missing required columns in trajectory file: %s", strjoin(missingCols, ", "));
end

flight_id = string(T.id);
epoch_ts  = toNum(T.timestamp);
lat       = toNum(T.latitude);
lon       = toNum(T.longitude);

good = ~ismissing(flight_id) & ~isnan(epoch_ts) & ~isnan(lat) & ~isnan(lon);
flight_id = flight_id(good);
epoch_ts  = epoch_ts(good);
lat       = lat(good);
lon       = lon(good);

fprintf("Rows read (valid): %d\n", numel(lon));
uIDs = unique(flight_id, 'stable');
fprintf("Unique flights:    %d\n", numel(uIDs));

%% ================= Plot ALL trajectories + sector boundary ==============
fig = figure('Visible','off');
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');

for i = 1:numel(uIDs)
    idx = (flight_id == uIDs(i));
    if nnz(idx) < 2
        continue
    end

    ii = find(idx);
    [~, ord] = sort(epoch_ts(ii));
    ii = ii(ord);

    %plot(ax, lon(ii), lat(ii), 'LineWidth', 1);
end

% Sector boundary on top
plot(ax, sector_lon, sector_lat, 'LineWidth', 2);

adjacent_sectors = function_create_adjacent_sectors()
plot_dummy_adjacent_sectors(ax, adjacent_sectors)

xlabel(ax, 'Longitude');
ylabel(ax, 'Latitude');
axis(ax, 'equal');

% Limits include both sector and trajectories
%x_all = [sector_lon(:); lon(:)];
%y_all = [sector_lat(:); lat(:)];
%xlim(ax, [min(x_all), max(x_all)]);
%ylim(ax, [min(y_all), max(y_all)]);

% Limits based on adjacent sectors
adjacent_sectors = function_create_adjacent_sectors;

x_all = [];
y_all = [];

for i = 1:numel(adjacent_sectors)
    secs = adjacent_sectors{i};
    for j = 1:numel(secs)
        coords = squeeze(secs(j).geometry.coordinates); % N x 2
        x_all = [x_all; coords(:,1)];
        y_all = [y_all; coords(:,2)];
    end
end

xlim(ax, [min(x_all), max(x_all)]);
ylim(ax, [min(y_all), max(y_all)]);

% Geographic-looking proportions (lon scaled by cos(lat))
meanLatAll = mean(y_all, "omitnan");
if ~isnan(meanLatAll)
    daspect(ax, [1/cosd(meanLatAll), 1, 1]);
end

%title(ax, "All trajectories (2D) with ESMM31 sector boundary", 'Interpreter','none');
title(ax, "Dummy adjacent sectors to ESMM31 sector", 'Interpreter','none');

%out_path = fullfile(fig_root, 'all_trajectories_2D.png');
out_path = fullfile(fig_root, 'dummy_adj_sectors.png');
exportgraphics(fig, out_path, 'Resolution', dpi);
close(fig);

fprintf("Done. Saved plot to: %s\n", out_path);

%% ========================= Local function ===============================
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


function plot_dummy_adjacent_sectors(ax, adjacent_sectors_data)
% Plot dummy adjacent sectors created by function_create_adjacent_sectors
% Simple lon/lat line plotting only.

if isempty(adjacent_sectors_data)
    return
end

hold(ax, 'on');

for i = 1:numel(adjacent_sectors_data)
    secs = adjacent_sectors_data{i};
    if isempty(secs), continue; end

    for j = 1:numel(secs)
        if ~isfield(secs(j),'geometry') || ~isfield(secs(j).geometry,'coordinates')
            continue
        end

        % coordinates is 1 x N x 2
        coords = secs(j).geometry.coordinates;
        coords = squeeze(coords);    % N x 2

        lon = coords(:,1);
        lat = coords(:,2);

        plot(ax, lon, lat, ...
            'k--', ...        % dashed black outline
            'LineWidth', 1);
    end
end
end



function adjacent_sectors_data = function_create_adjacent_sectors
% Create variable with adjacent sector data

% adjacent_sectors_data is a 1xN cell (N is the total number of
% adjacent sectors). Each cell contains the information of the M subsectors
% that form each sector, wich is a Mx1 structure with fields:
% properties: (1x1 structure with fields DESIGNATOR, UPPER_LIMIT_VALUE, LOWER_LIMIT_VALUE)
% geometry (1x1 structure with field coordinates)

N = 4;
adjacent_sectors_data = cell(1,N); % Initialize output variable


adjacent_sectors_data{1}(1,1) = struct('properties',[], 'geometry', []);

adjacent_sectors_data{1}(1).properties.DESIGNATOR = 'dummy1';
adjacent_sectors_data{1}(1).properties.LOWER_LIMIT_VALUE = 0; 
adjacent_sectors_data{1}(1).properties.UPPER_LIMIT_VALUE = 660;

% South-West
dummy1_lon = [12.842222 12.808611 12.135000 11.955833 11.355833 12.842222 12.842222];
dummy1_lat = [55.516944 55.543889 56.075833 56.523333 56.523333 54.916944 55.516944];

adjacent_sectors_data{1}(1).geometry.coordinates = zeros(1,length(dummy1_lat),2);
adjacent_sectors_data{1}(1).geometry.coordinates(:,:,1) = dummy1_lon';
adjacent_sectors_data{1}(1).geometry.coordinates(:,:,2) = dummy1_lat';


adjacent_sectors_data{2}(1,1) = struct('properties',[], 'geometry', []);

adjacent_sectors_data{2}(1).properties.DESIGNATOR = 'dummy2';
adjacent_sectors_data{2}(1).properties.LOWER_LIMIT_VALUE = 0; 
adjacent_sectors_data{2}(1).properties.UPPER_LIMIT_VALUE = 660;

% South-East
dummy2_lon = [15.875278 15.666389 15.242222 15.117778 14.889167 14.536667 13.703889 13.347222 13.093333 12.893611 12.842222 12.842222 16.475278 15.875278];
dummy2_lat = [56.957500 56.866667 56.683611 56.628056 56.501667 56.316667 55.880556 55.766667 55.642500 55.542500 55.516944 54.916944 56.957500 56.957500];

adjacent_sectors_data{2}(1).geometry.coordinates = zeros(1,length(dummy2_lat),2);
adjacent_sectors_data{2}(1).geometry.coordinates(:,:,1) = dummy2_lon';
adjacent_sectors_data{2}(1).geometry.coordinates(:,:,2) = dummy2_lat';


adjacent_sectors_data{3}(1,1) = struct('properties',[], 'geometry', []);

adjacent_sectors_data{3}(1).properties.DESIGNATOR = 'dummy3';
adjacent_sectors_data{3}(1).properties.LOWER_LIMIT_VALUE = 0; 
adjacent_sectors_data{3}(1).properties.UPPER_LIMIT_VALUE = 660;

% North-East
dummy3_lon = [14.327500 15.875278 16.475278 14.327500 14.327500];
dummy3_lat = [57.838889 56.957500 56.957500 58.438889 57.838889];

adjacent_sectors_data{3}(1).geometry.coordinates = zeros(1,length(dummy3_lat),2);
adjacent_sectors_data{3}(1).geometry.coordinates(:,:,1) = dummy3_lon';
adjacent_sectors_data{3}(1).geometry.coordinates(:,:,2) = dummy3_lat';


adjacent_sectors_data{4}(1,1) = struct('properties',[], 'geometry', []);

adjacent_sectors_data{4}(1).properties.DESIGNATOR = 'dummy4';
adjacent_sectors_data{4}(1).properties.LOWER_LIMIT_VALUE = 0; 
adjacent_sectors_data{4}(1).properties.UPPER_LIMIT_VALUE = 660;

% North-West
dummy4_lon = [11.955833 12.069444 12.552222 12.859722 13.702500 14.327500 14.327500 11.355833 11.955833];
dummy4_lat = [56.523333 56.600278 56.922778 57.084722 57.516389 57.838889 58.438889 56.523333 56.523333];

adjacent_sectors_data{4}(1).geometry.coordinates = zeros(1,length(dummy4_lat),2);
adjacent_sectors_data{4}(1).geometry.coordinates(:,:,1) = dummy4_lon';
adjacent_sectors_data{4}(1).geometry.coordinates(:,:,2) = dummy4_lat';

end


