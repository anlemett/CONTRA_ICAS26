function AC = function_read_trajectories()

thisFileDir = fileparts(mfilename('fullpath'));
traj_file = fullfile(thisFileDir, 'code_input', 'original_all_flights_in_ESMM31_day28_extended.csv');
assert(exist(traj_file,'file')==2, "Trajectory file not found: %s", traj_file);

% Read header + data as strings (space/tab delimited)
fid = fopen(traj_file,'r');
assert(fid>0, "Could not open trajectory file: %s", traj_file);

hdr = strtrim(fgetl(fid));
assert(ischar(hdr) && ~isempty(hdr), "Could not read header line from: %s", traj_file);

names = regexp(hdr, '\s+', 'split');
fmt   = repmat('%s', 1, numel(names));
C = textscan(fid, fmt, ...
    'Delimiter', {' ','\t'}, ...
    'MultipleDelimsAsOne', true, ...
    'CollectOutput', true);
fclose(fid);

assert(~isempty(C) && ~isempty(C{1}), "No trajectory data rows were read from: %s", traj_file);
T = cell2table(C{1}, 'VariableNames', matlab.lang.makeUniqueStrings(names));

req = ["id","timestamp","latitude","longitude","altitude"];
assert(all(ismember(req, string(T.Properties.VariableNames))), ...
    "Trajectory file must contain columns: id, timestamp, latitude, longitude, altitude");

% Parse + clean invalid/missing values
toNum = @(x) str2double(strrep(string(x), ",", "."));

id  = string(T.id);
ts  = toNum(T.timestamp);     % epoch seconds
lat = toNum(T.latitude);
lon = toNum(T.longitude);
alt = toNum(T.altitude);      % feet

good = ~ismissing(id) & ~isnan(ts) & ~isnan(lat) & ~isnan(lon) & ~isnan(alt);
id=id(good); ts=ts(good); lat=lat(good); lon=lon(good); alt=alt(good);

% Group by flight id and build AC
[uIDs, ~, g] = unique(id, 'stable');
AC = struct('id', cell(numel(uIDs),1), 'WP', cell(numel(uIDs),1));

for k = 1:numel(uIDs)
    idx = (g == k);

    t  = ts(idx);
    lo = lon(idx);
    la = lat(idx);
    al = alt(idx) / 100;   % convert ft -> FL

    [t, ord] = sort(t);
    AC(k).id = uIDs(k);
    AC(k).WP = [t, lo(ord), la(ord), al(ord)];
end

end
