clear; clc;

%% -------------------- Paths / settings --------------------
parquet_file = fullfile('.', 'code_input', 'grid_era5_smoothed_day28_ESMM31.parquet');

%% -------------------- Read parquet --------------------
if ~exist(parquet_file, 'file')
    error("Parquet file not found: %s", parquet_file);
end

P = parquetread(parquet_file);

%% -------------------- Validate timestamp column --------------------
if ~ismember("timestamp", string(P.Properties.VariableNames))
    error("Column 'timestamp' not found in the Parquet file.");
end

%% -------------------- Extract unique timestamps --------------------
ts = string(P.timestamp);
ts = ts(~ismissing(ts));

tsVals = unique(ts, "stable");

fprintf("Number of unique timestamps: %d\n", numel(tsVals));
disp(tsVals);

%% -------------------- Check time spacing (if parseable) --------------------
try
    tsDT = datetime(tsVals);   % let MATLAB auto-detect format
    dt = diff(tsDT);

    fprintf("Unique time differences between consecutive timestamps:\n");
    disp(unique(dt));

    fprintf("First timestamp: %s\n", string(tsDT(1)));
    fprintf("Last  timestamp: %s\n", string(tsDT(end)));

catch ME
    warning(ME.identifier, '%s', ME.message);
    fprintf("Timestamps were printed as raw strings above.\n");
end
