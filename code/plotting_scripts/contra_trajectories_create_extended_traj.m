clear; clc;

%% ============================ SETTINGS ==================================
traj_file = fullfile('.', 'code_input', 'original_all_flights_in_ESMM31_day28.csv');

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

%% ================= Extend trajectories to cross sector border ===========
% For each flight, we keep the portion inside the sector (longest
% contiguous inside run) and extend it in both directions so it crosses the
% sector boundary and continues for a few additional points outside.
extra_points_outside = 5;   % number of points to keep after crossing border (each direction)
max_extend_steps     = 200; % safety cap per direction

in_sector_all = inpolygon(lon, lat, sector_lon, sector_lat);

% Build an extended trajectory table (same variable names as input T)
varNames = T.Properties.VariableNames;

% Build an extended trajectory table (same variable names as input T)
Tout = T([],:); % empty table with same variables

for i = 1:numel(uIDs)
    idx = (flight_id == uIDs(i));
    if nnz(idx) < 2
        continue
    end

    ii = find(idx);

    % Determine contiguous "inside sector" runs and keep the longest
    inside_i = in_sector_all(ii);

    if ~any(inside_i)
        continue
    end

    % Find runs of true values in inside_i
    d = diff([false; inside_i; false]);
    run_starts = find(d == 1);
    run_ends   = find(d == -1) - 1;

    run_lengths = run_ends - run_starts + 1;
    [~, kmax] = max(run_lengths);
    rs = run_starts(kmax);
    re_ = run_ends(kmax);

    ii_seg = ii(rs:re_);

    % Segment arrays
    lon_seg = lon(ii_seg);
    lat_seg = lat(ii_seg);
    t_seg   = epoch_ts(ii_seg);

    if numel(lon_seg) < 2
        continue
    end

    % Estimate time step (median)
    dt = median(diff(t_seg));
    if ~isfinite(dt) || dt <= 0
        dt = 1;
    end

    % Backward heading (use first two points)
    v_b = [lon_seg(2)-lon_seg(1), lat_seg(2)-lat_seg(1)];
    if norm(v_b) == 0
        v_b = [1, 0];
    end
    v_b = v_b / norm(v_b);

    % Forward heading (use last two points)
    v_f = [lon_seg(end)-lon_seg(end-1), lat_seg(end)-lat_seg(end-1)];
    if norm(v_f) == 0
        v_f = [1, 0];
    end
    v_f = v_f / norm(v_f);

    % Extend backward
    lon_b = [];
    lat_b = [];
    n_out = 0;
    p = [lon_seg(1), lat_seg(1)];
    for s = 1:max_extend_steps
        p = p - v_b;  % step (in same units as lon/lat)
        lon_b(end+1,1) = p(1);
        lat_b(end+1,1) = p(2);
        inside = inpolygon(p(1), p(2), sector_lon, sector_lat);
        if ~inside
            n_out = n_out + 1;
            if n_out >= extra_points_outside
                break
            end
        end
    end
    lon_b = flipud(lon_b);
    lat_b = flipud(lat_b);
    t_b = (t_seg(1) - dt*(numel(lon_b):-1:1))';

    % Extend forward
    lon_f = [];
    lat_f = [];
    n_out = 0;
    p = [lon_seg(end), lat_seg(end)];
    for s = 1:max_extend_steps
        p = p + v_f;  % step
        lon_f(end+1,1) = p(1);
        lat_f(end+1,1) = p(2);
        inside = inpolygon(p(1), p(2), sector_lon, sector_lat);
        if ~inside
            n_out = n_out + 1;
            if n_out >= extra_points_outside
                break
            end
        end
    end
    t_f = (t_seg(end) + dt*(1:numel(lon_f)))';

    lon_ext = [lon_b; lon_seg; lon_f];
    lat_ext = [lat_b; lat_seg; lat_f];
    t_ext   = [t_b;   t_seg;   t_f];

    % Create an output block that preserves all original columns
    n_b = numel(lon_b);
    n_s = numel(lon_seg);
    n_f = numel(lon_f);

    idx_block = [repmat(ii_seg(1), n_b, 1); ii_seg(:); repmat(ii_seg(end), n_f, 1)];
    block = T(idx_block, :);

    % Overwrite mandatory fields; other fields (e.g., altitude) are preserved by replication
    if iscell(block.id)
        block.id = repmat({char(uIDs(i))}, height(block), 1);
    else
        block.id = repmat(uIDs(i), height(block), 1);
    end

    block.timestamp = arrayfun(@(x)sprintf('%.0f', x), t_ext, 'UniformOutput', false);
    block.latitude  = arrayfun(@(x)sprintf('%.8f', x), lat_ext, 'UniformOutput', false);
    block.longitude = arrayfun(@(x)sprintf('%.8f', x), lon_ext, 'UniformOutput', false);

    Tout = [Tout; block];
end

% Save extended trajectories CSV next to the input file
[in_dir, in_base, ~] = fileparts(traj_file);
out_csv = fullfile(in_dir, in_base + "_extended.csv");
writetable(Tout, out_csv, 'Delimiter',' ', 'WriteVariableNames', true);
fprintf("Saved extended trajectories CSV to: %s\n", out_csv);

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

function seg = local_longest_true_run(mask)
% Return indices (within mask) of the longest contiguous run of true.
    if isempty(mask) || ~any(mask)
        seg = [];
        return
    end

    d = diff([false; mask(:); false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    [~, k] = max(ends - starts + 1);
    seg = starts(k):ends(k);
end

function [lon_ext, lat_ext, t_ext] = local_extend_polyline_across_border(lon_seg, lat_seg, t_seg, sector_lon, sector_lat, extra_outside, max_steps)
% Extend a polyline segment in both directions so it crosses the polygon
% border and continues for extra_outside points beyond.

    lon_seg = lon_seg(:);
    lat_seg = lat_seg(:);
    t_seg   = t_seg(:);

    if numel(lon_seg) < 2
        lon_ext = lon_seg; lat_ext = lat_seg; t_ext = t_seg;
        return
    end

    % Time step for synthetic points
    dt = median(diff(t_seg), 'omitnan');
    if isempty(dt) || ~isfinite(dt) || dt <= 0
        dt = 1;
    end

    % Typical spatial step near ends
    d_back = [lon_seg(2)-lon_seg(1), lat_seg(2)-lat_seg(1)];
    d_forw = [lon_seg(end)-lon_seg(end-1), lat_seg(end)-lat_seg(end-1)];
    s_back = norm(d_back);
    s_forw = norm(d_forw);
    if s_back == 0
        s_back = max(norm([lon_seg(min(3,end))-lon_seg(1), lat_seg(min(3,end))-lat_seg(1)]), eps);
    end
    if s_forw == 0
        s_forw = max(norm([lon_seg(end)-lon_seg(max(end-2,1)), lat_seg(end)-lat_seg(max(end-2,1))]), eps);
    end

    u_back = -d_back ./ max(s_back, eps);
    u_forw =  d_forw ./ max(s_forw, eps);

    % Backward extension
    p = [lon_seg(1), lat_seg(1)];
    n_out = 0;
    lon_b = []; lat_b = [];
    for k = 1:max_steps
        p = p + u_back .* s_back;
        lon_b(end+1,1) = p(1);
        lat_b(end+1,1) = p(2);
        inside = inpolygon(p(1), p(2), sector_lon, sector_lat);
        if ~inside
            n_out = n_out + 1;
            if n_out >= extra_outside
                break
            end
        end
    end
    lon_b = flipud(lon_b);
    lat_b = flipud(lat_b);
    t_b   = (t_seg(1) - dt*(numel(lon_b):-1:1))';

    % Forward extension
    p = [lon_seg(end), lat_seg(end)];
    n_out = 0;
    lon_f = []; lat_f = [];
    for k = 1:max_steps
        p = p + u_forw .* s_forw;
        lon_f(end+1,1) = p(1);
        lat_f(end+1,1) = p(2);
        inside = inpolygon(p(1), p(2), sector_lon, sector_lat);
        if ~inside
            n_out = n_out + 1;
            if n_out >= extra_outside
                break
            end
        end
    end
    t_f = (t_seg(end) + dt*(1:numel(lon_f)))';

    lon_ext = [lon_b; lon_seg; lon_f];
    lat_ext = [lat_b; lat_seg; lat_f];
    t_ext   = [t_b;   t_seg;   t_f];
end
