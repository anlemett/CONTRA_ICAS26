% INPUT:
%   LOWER_FL, UPPER_FL : flight levels (inclusive range), e.g. 300 and 340
%   t0_dt, t1_dt       : datetimes
%
% OUTPUT:
%   cost_obstacles is cell(T,1), where T is the number of hourly timestamps
%   between t0_dt and t1_dt (inclusive, hour-aligned).
%   Each cost_obstacles{t} is a cell(1, num_poly). Each entry is a struct:
%       .pgon     (polyshape in lon/lat degrees)
%       .altitude (ft)
function cost_obstacles = function_get_contrail_polygons(LOWER_FL, UPPER_FL, t0_dt, t1_dt)

lower_ft = LOWER_FL * 100;
upper_ft = UPPER_FL * 100;

thisFileDir = fileparts(mfilename('fullpath'));
cost_file = fullfile(thisFileDir, 'code_input', 'grid_era5_smoothed_day28_ESMM31.csv');

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

toNum = @(x) str2double(strrep(string(x), ",", "."));

Tc.latitude  = toNum(Tc.latitude);
Tc.longitude = toNum(Tc.longitude);
Tc.cost      = toNum(Tc.cost);
Tc.altitude  = toNum(Tc.altitude);
tsKey_cost   = string(Tc.timestamp);

validc = ~ismissing(tsKey_cost) & ~isnan(Tc.latitude) & ~isnan(Tc.longitude) & ...
         ~isnan(Tc.cost) & ~isnan(Tc.altitude);
Tc = Tc(validc,:);
tsKey_cost = tsKey_cost(validc);

% Inclusive altitude band (ft)
cost_alt_mask_all = (Tc.altitude >= lower_ft) & (Tc.altitude <= upper_ft);

%% ============= Create cost obstacles for each hour in [t0_dt, t1_dt] ==============
threshold  = 0.700;
conn       = 4;

keys_win = local_cost_keys_utc_window(t0_dt, t1_dt);

%cost in the beginnig of the hour defines the obstacles for this hour
T = numel(keys_win)-1; % last time is not inkcluded

cost_obstacles = cell(T, 1);

for it = 1:T
    key_t = keys_win(it);

    in_time_cost = (tsKey_cost == key_t) & cost_alt_mask_all;

    cost_obstacles{it} = cell(1, 0);
    if ~any(in_time_cost)
        continue
    end

    Th_all = Tc(in_time_cost, :);
    alts = sort(unique(Th_all.altitude));

    tmp = {};
    for ia = 1:numel(alts)
        alt_ft = alts(ia);
        Th = Th_all(Th_all.altitude == alt_ft, :);
        if isempty(Th)
            continue
        end

        allLatitudes  = sort(unique(Th.latitude));
        allLongitudes = sort(unique(Th.longitude));
        if isempty(allLatitudes) || isempty(allLongitudes)
            continue
        end

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

            tmp{end+1} = struct('pgon', pg, 'altitude', alt_ft);
        end
    end

    if ~isempty(tmp)
        cost_obstacles{it} = reshape(tmp, 1, []);
    end
end

end

%% ========================= Local functions ==============================

function key = local_cost_key_utc(tUTC)
try
    tUTC.TimeZone = "UTC";
catch
end
key = string(datestr(tUTC, 'yyyy-mm-dd HH:MM:SS')) + "+00:00";
end

function keys = local_cost_keys_utc_window(t0_dt, t1_dt)
try
    t0_dt.TimeZone = "UTC";
catch
end
try
    t1_dt.TimeZone = "UTC";
catch
end

t0h = dateshift(t0_dt, 'start', 'hour');
t1h = dateshift(t1_dt, 'start', 'hour');

if t1h < t0h
    keys = string.empty(0,1);
    return
end

tvec = (t0h:hours(1):t1h).';
keys = string(datestr(tvec, 'yyyy-mm-dd HH:MM:SS')) + "+00:00";
end
