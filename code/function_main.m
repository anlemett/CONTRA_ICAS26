function [ASCR, sector_names, sector_time, sector_data] = function_main()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contra
% MATLAB version: MATLAB R2025a
% 
% Inputs:
%   1) Contrails cost matrix
% 
%   2) Flight data: Flight trajectories from OpenSky
%
%   3) ESMM31 sector coordinates
% 
% Output:
%   Available sector capacity ratio ASCR{k}(t,1) for the main sector, time t

% Time 

tz = "UTC";

minut_vec = 0:60:60;

t0_dt = datetime(2024,1,28,18,0,0,'TimeZone',tz);

t_vec_ini_dt   = t0_dt + minutes(minut_vec);
t_vec_fin_dt   = t_vec_ini_dt + hours(1);

t_vec_ini = posixtime(t_vec_ini_dt); % convert to Unix timestamps
t_vec_fin = posixtime(t_vec_fin_dt);

t1_dt = t_vec_fin_dt(end);

global TARGET_FL

%% ========================= Read ESMM31 sector ===========================

thisFileDir = fileparts(mfilename('fullpath'));

sector_file = fullfile(thisFileDir, 'code_input', 'ESMM31.geojson');

if exist(sector_file, 'file') ~= 2
    error("Sector file not found: %s", sector_file);
end

G = jsondecode(fileread(sector_file));
if ~isfield(G, "features") || isempty(G.features)
    error("Invalid GeoJSON: missing features.");
end

feat = G.features(1);
feat.properties.LOWER_LIMIT_VALUE = TARGET_FL - 5;
feat.properties.UPPER_LIMIT_VALUE = TARGET_FL + 5;
main_sectors = cell(1,1);
main_sectors{1} = feat;

sector_names = {main_sectors{1}.properties.DESIGNATOR};
num_of_hours = 2;
sector_time  = 1:num_of_hours;
sector_data = repmat(main_sectors(1), num_of_hours, 1);

adjacent_sectors = function_create_adjacent_sectors;


%% ========================= Read COST GRID ===============================

weather_polygons = function_get_contrail_polygons(TARGET_FL, TARGET_FL, t0_dt, t1_dt);

[nT,nM] = size(weather_polygons); % nT: times - nM: 1?


%% ================== Read TRAJECTORIES and build AC(a).WP ================

AC = function_read_trajectories();

%% ========================= ASCR computation =============================

    nk = numel(sector_names);
    ASCR = cell(nk, 1);
    
    %for k = 1:nk % For each sector
    for k = 1 % For one sector
    
        sectors_t = 0; % initialize variable sectors_t
    
        ASCR{k} = nan(nT, nM);

        for t = 1:nT % For each time
        %for t = 1 % Only for the first time (18:00)
    
                if ~isequal(sector_time(t), sectors_t) 
                    sectors_t = sector_time(t);
                    sector_data_t = sector_data(sectors_t);
                    [sector_ab, a_band, flows_j] = function_flows_sector_k(k, sector_data_t, adjacent_sectors);
                    [p_in, p_out, AC_in] = function_p_in_out(AC, sector_ab, a_band);
                end
        
                if t == 1 % initialies Wij_m
                    Wij_m1 = function_Wij_ini(flows_j);
                end
        
                [~, Wij, total_ac] = function_Wj(t_vec_ini(t), t_vec_fin(t), p_in, p_out, a_band, flows_j, AC_in);
        
                if total_ac == 0 % If there are no aircraft in the sector, use the previous weights
                    Wij = Wij_m1;
                end
        
                for m = 1 % For one member
                    %disp([k,t,m])
                    ASCR_k = function_ASCR_k(sector_ab, flows_j, weather_polygons{t,m}, Wij, a_band);
                    ASCR{k}(t,m) = ASCR_k;
                end
        
                Wij_m1 = Wij;
        
        end

    end

end
