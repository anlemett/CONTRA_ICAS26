function adjacent_sectors_data = function_create_adjacent_sectors

thisFileDir = fileparts(mfilename('fullpath'));
inputDir    = fullfile(thisFileDir, 'code_input');

sector_files = {
    'ESMM41.geojson'
    'ESMM53.geojson'
    'ESMM54.geojson'
    'ESMM61.geojson'
    'ESMM91.geojson'
    'ESMMW1.geojson'
    'DENMARK_DUMMY.geojson'
};

N = 7;
adjacent_sectors_data = cell(1,N); % Initialize output variable

for f = 1:N

    G = jsondecode(fileread(fullfile(inputDir, sector_files{f})));

    feat = G.features(1);

    adjacent_sectors_data{f}(1,1) = struct('properties',[], 'geometry', []);

    % properties
    adjacent_sectors_data{f}(1).properties = feat.properties;

    % geometry (single ring, already closed)
    coords = feat.geometry.coordinates;

    lon = coords(:,:,1);
    lat = coords(:,:,2);

    adjacent_sectors_data{f}(1).geometry.coordinates = zeros(1,length(lat),2);
    adjacent_sectors_data{f}(1).geometry.coordinates(:,:,1) = lon';
    adjacent_sectors_data{f}(1).geometry.coordinates(:,:,2) = lat';

end

end
