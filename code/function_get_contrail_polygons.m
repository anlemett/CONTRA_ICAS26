% Read weather data and store it as polygons

function cost_obstacles = function_get_contrail_polygons(TARGET_FL, t0_dt, t1_dt)

target_ft   = TARGET_FL * 100;

COST_ALT_FT = target_ft;

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

cost_alt_mask_all = (Tc.altitude == COST_ALT_FT);

%% ============= Create cost obstacles for the selected hour ==============
threshold  = 0.700;
conn       = 4;

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

end

%% ========================= Local functions ==============================

function key = local_cost_key_utc(tUTC)
try
    tUTC.TimeZone = "UTC";
catch
end
key = string(datestr(tUTC, 'yyyy-mm-dd HH:MM:SS')) + "+00:00";
end
