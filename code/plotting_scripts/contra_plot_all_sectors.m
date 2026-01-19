%% Plot all sectors from GeoJSON files (script, not a function)
clear; clc;

% Folder containing this script
thisFileDir = fileparts(mfilename('fullpath'));

% Input folder with .geojson sector files
inputDir = fullfile(thisFileDir, 'code_input');

% Sector files to plot (added ESMM31.geojson)
sector_files = {
    'ESMM31.geojson'
    'ESMM41.geojson'
    'ESMM53.geojson'
    'ESMM54.geojson'
    'ESMM61.geojson'
    'ESMM91.geojson'
    'ESMMW1.geojson'
    %'DENMARK_DUMMY.geojson'
    'KOBENHAVN_CTA_A.geojson'
};

% Figure setup
figure('Color','w'); hold on; grid on;
xlabel('Longitude'); ylabel('Latitude');

% Plot each sector
for f = 1:numel(sector_files)

    fp = fullfile(inputDir, sector_files{f});
    if ~isfile(fp)
        warning('File not found: %s (skipping)', fp);
        continue;
    end

    G = jsondecode(fileread(fp));

    if ~isfield(G,'features') || isempty(G.features)
        warning('No features in: %s (skipping)', sector_files{f});
        continue;
    end

    feat = G.features(1);

    if ~isfield(feat,'geometry') || ~isfield(feat.geometry,'coordinates')
        warning('No geometry/coordinates in: %s (skipping)', sector_files{f});
        continue;
    end

    coords = feat.geometry.coordinates;

    % Expecting GeoJSON Polygon ring:
    % coords is typically 1xN x2 (lon/lat) or Nx2 depending on decoder.
    % Handle common cases robustly.
    try
        lon = coords(:,:,1);
        lat = coords(:,:,2);

        lon = lon(:);
        lat = lat(:);
    catch
        % Fallback if coords is Nx2
        lon = coords(:,1);
        lat = coords(:,2);
    end

    % Ensure ring is closed for plotting
    if ~isempty(lon) && (lon(1) ~= lon(end) || lat(1) ~= lat(end))
        lon(end+1) = lon(1);
        lat(end+1) = lat(1);
    end

    % Plot outline
    plot(lon, lat, '-', 'LineWidth', 1.5);

    % Optional label (uses properties.id/name if present; otherwise filename)
    [~, label, ~] = fileparts(sector_files{f});

    if isfield(feat,'properties') && isstruct(feat.properties)
        if isfield(feat.properties,'id') && ~isempty(feat.properties.id)
            label = string(feat.properties.id);
        elseif isfield(feat.properties,'name') && ~isempty(feat.properties.name)
            label = string(feat.properties.name);
        end
    end

    % Place label at centroid (mean of vertices; sufficient for display)
    if ~isempty(lon) && ~isempty(lat)
        text(mean(lon,'omitnan'), mean(lat,'omitnan'), char(label), ...
            'FontSize', 9, 'Interpreter', 'none', ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
end

xl = xlim;          % get current x-axis limits [xmin xmax]
xlim([xl(1)-1, xl(2)]);

axis equal;
hold off;
