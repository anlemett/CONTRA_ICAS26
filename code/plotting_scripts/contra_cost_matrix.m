clear; clc;

csv_file = fullfile('.', 'code_input', 'grid_era5_smoothed_day28_ESMM31.csv');

opts = detectImportOptions(csv_file, "Delimiter", ",");
opts.Delimiter = ",";
opts.VariableNamesLine = 1;
opts.VariableNamingRule = "preserve";

disp(opts.VariableNames);

% Exact match using char vectors
iTs = find(strcmpi(opts.VariableNames, 'timestamp'), 1);
if isempty(iTs)
    error("No variable named 'timestamp' found in opts.VariableNames.");
end

% Force timestamp to be read as text BEFORE readtable
opts.VariableTypes{iTs} = 'char';

T = readtable(csv_file, opts);
fprintf("%s\n", T.Properties.VariableNames{:});

disp(T(1,:))
tsKey  = string(T.timestamp); 

% --- Fix numeric columns explicitly (works for cell/string/categorical/numeric) ---
toNum = @(x) str2double(strrep(string(x), ",", "."));  % also handles decimal comma

T.latitude  = toNum(T.latitude);
T.longitude = toNum(T.longitude);
T.cost      = toNum(T.cost);
T.altitude  = toNum(T.altitude);


% --- Remove rows with invalid values ---
valid = ~ismissing(tsKey) & ~isnan(T.altitude) & ...
        ~isnan(T.latitude) & ~isnan(T.longitude) & ~isnan(T.cost);

% --- Select first valid timestamp and altitude (just for example) ---
tsVals  = unique(tsKey(valid), "stable");
altVals = unique(T.altitude(valid), "stable");

ts0  = tsVals(1);
alt0 = altVals(1);

% --- Filter ---
Tf = T(tsKey == ts0 & T.altitude == alt0, :);

% --- Build cost matrix ---
allLatitudes = sort(unique(Tf.latitude));
allLongitudes = sort(unique(Tf.longitude));

[~, latIdx] = ismember(Tf.latitude, allLatitudes);
[~, lonIdx] = ismember(Tf.longitude, allLongitudes);

% --- Build numeric cost matrix (not binary) ---
cost_matrix = nan(numel(allLatitudes), numel(allLongitudes));

linIdx = sub2ind(size(cost_matrix), latIdx, lonIdx);
cost_matrix(linIdx) = Tf.cost;

% --- Visualize the cost matrix as a heatmap ---
figure;
heatmap(allLongitudes, allLatitudes, cost_matrix, 'ColorMap', parula);
xlabel('Longitude');
ylabel('Latitude');
title('Cost Matrix Heatmap');