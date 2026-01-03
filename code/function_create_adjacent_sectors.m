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

