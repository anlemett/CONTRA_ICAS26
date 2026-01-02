clear; clc;

%% -------------------- Paths --------------------
csv_file = fullfile('.', 'code_input', 'grid_era5_smoothed_day28_ESMM31.csv');
out_dir  = fullfile('.', 'figures');

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

out_png = fullfile(out_dir, 'cost_histogram_all_altitudes_all_timestamps.png');

%% -------------------- Read CSV --------------------
if ~exist(csv_file, 'file')
    error("CSV file not found: %s", csv_file);
end

opts = detectImportOptions(csv_file, "Delimiter", ",");
opts.Delimiter = ",";
opts.VariableNamesLine = 1;
opts.VariableNamingRule = "preserve";

T = readtable(csv_file, opts);

% Required column
if ~ismember("cost", string(T.Properties.VariableNames))
    error("Column 'cost' not found in CSV file.");
end

%% -------------------- Convert cost values --------------------
% Handle numeric / string / decimal comma consistently
toNum = @(x) str2double(strrep(string(x), ",", "."));
cost_vals = toNum(T.cost);

% Remove invalid values
cost_vals = cost_vals(~isnan(cost_vals) & ~isinf(cost_vals));

if isempty(cost_vals)
    error("No valid cost values found.");
end

%% -------------------- Histogram --------------------
fig = figure('Visible','on');

histogram(cost_vals, ...
    'BinMethod', 'fd', ...      % Freedman–Diaconis rule
    'Normalization', 'count');

xlabel('Cost value');
ylabel('Count');
title('Histogram of cost values (all altitudes, all timestamps)');

grid on;

%% -------------------- Save figure --------------------
exportgraphics(fig, out_png, 'Resolution', 200);

fprintf('Histogram saved to:\n%s\n', out_png);
