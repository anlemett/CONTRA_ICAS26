function [p_in, p_out, AC_in] = function_p_in_out(AC, sector_ab_k, a_band)

% This function computes the entry (p_in) and exit (p_out) points of the
% aircraft trajectories described by the flight plans (AC) in the
% sector (sector_ab_k)

% Flight plan data
Nac = numel(AC);

% Select aircraft that (can) intersect with the sector
sector_big = union([sector_ab_k{:}]);
in_sector  = false(Nac,1);

for a = 1:Nac
    % Does the aircraft intersect the sector horizontaly?
    lineseg = [AC(a).WP(:,2),AC(a).WP(:,3)];
    [in,~] = intersect(sector_big,lineseg);
    ok_horizontal = ~isempty(in);
    
    % Does the aircraft fly within the sector altitude range?
    max_alt = max(AC(a).WP(:,4)); % max altitude
    min_alt = min(AC(a).WP(:,4)); % min altitude
    ok_vertical = ~((max_alt<=(a_band(1)-5))||(min_alt>(a_band(end)+5)));
    if ok_horizontal&&ok_vertical
        in_sector(a) = true;
    end
end

AC_in = AC(in_sector);


% Detect those trajectories that begin or end within the sector: Not
% considered -> We should have 0 of these! (the aircraft do not materialize
% or dissapear) If these is any, it should be a problem with the flight
% plan data.

index_pass = function_pass_through(AC_in, sector_ab_k, a_band);

AC_in = AC_in(index_pass); % Discard those trajectories that start or end within the sector
Nac = numel(AC_in);

% Trajectories intersections with sector k
p_in  = nan(Nac, 4); % Entrance point
p_out = p_in;          % Exit point

% figure, hold on
% plot(sector_ab_k{1})
% plot(sector_ab_k{40})
for a  = 1:Nac

    [p_in_a, p_out_a] = function_sector_intersection(AC_in(a).WP, sector_ab_k, a_band);
    
    p_in(a,:)  = p_in_a;
    p_out(a,:) = p_out_a;
    
    
%     plot3(p_in_a(2), p_in_a(3), p_in_a(4), 'r*')
%     plot3(p_out_a(2), p_out_a(3),p_out_a(4), 'bo')
%     plot3(AC_in(a).WP(:,2),AC_in(a).WP(:,3),AC_in(a).WP(:,4))

end

end

function [sector_div, a_band_div] = function_divide_sector(sector_ab, a_band)
% This function divides the 3D sector (sector_ab) into groups of 2D polygons with the
% same shape (sector_div) and the corresponding altitude division
% (a_band_div) index

[nab,~] = size(a_band); % number of altitude bands

sector_div = cell(size(sector_ab));
a_band_div = zeros(nab,2);

sector_div{1}   = sector_ab{1};
a_band_div(1,1) = 1;

i_aux  = 1;
i_aux2 = 1;
for i = 2:nab
   
   polyout = subtract(sector_ab{i_aux2},sector_ab{i});
   
   if polyout.NumRegions>0 % There is a shape difference
       a_band_div(i_aux,2) = i;
       i_aux = i_aux + 1;
       sector_div{i_aux} = sector_ab{i};
       a_band_div(i_aux,1) = i;
       
       i_aux2 = i;
   end
   
end

a_band_div(cellfun(@isempty,sector_div),:) = [];
a_band_div(end,2) = nab;
sector_div(cellfun(@isempty,sector_div)) = []; 

index = a_band_div;

a_band_div(:,1) = a_band(index(:,1));
a_band_div(:,2) = a_band(index(:,2));

end

function index = function_pass_through(AC, sector_ab_k, a_band)
% Function that detects those trajectories that begin or end within the
% sector:
% index: 1 - the aircraft passes through the sector
%        0 - The airctraft trajectory starts or ends within the sector

Nac = numel(AC);
index_first = false(Nac,1);
index_last  = false(Nac,1);
for a = 1:Nac
    
    % Is the first waypoint inside the sector
    WP = [AC(a).WP(1,2), AC(a).WP(1,3), AC(a).WP(1,4)];
    
    if (WP(3)<a_band(1))||((WP(3)>a_band(end))) % If the waypoint altitude is outside the altitude bands limits
        index_first(a) = true;
    else
        % Find WP altitude band
        sector_poly = sector_ab_k{(a_band>=(WP(3)-5))&(a_band<(WP(3)+5))};
        in_index = inpolygon(WP(1),WP(2),sector_poly.Vertices(:,1),sector_poly.Vertices(:,2));
        if ~in_index % If the point is not inside the sector
            index_first(a) = true;
        end
    end
    % Is the last waypoint inside the sector
    WP  = [AC(a).WP(end,2), AC(a).WP(end,3), AC(a).WP(end,4)];
    
    if (WP(3)<a_band(1))||((WP(3)>a_band(end))) % If the waypoint altitude is outside the altitude bands limits
        index_last(a) = true;
    else
        % Find WP altitude band
        sector_poly = sector_ab_k{(a_band>=(WP(3)-5))&(a_band<(WP(3)+5))};
        in_index = inpolygon(WP(1),WP(2),sector_poly.Vertices(:,1),sector_poly.Vertices(:,2));
        if ~in_index % If the point is not inside the sector
            index_last(a) = true;
        end
    end
    
end

index = index_first&index_last;
end

function [p_in, p_out] = function_sector_intersection(WP, sector_ab, a_band)

% Function that calculates the entry (p_in) and exit (p_out) points of a
% tajectory WP(t, lon, lat, h) in sector k, given by polygons (sector_ab_k)
% at different altitude bands (a_band)
%
% If there is no intersection, p_in and p_out are nan.

% [nab,~] = size(a_band); % number of altitude bands

p_in = nan(1,4); p_out = nan(1,4); 

p_in_found = false; p_out_found = false; % Have we found p_in and p_out?
    
% Sector divisions in different geometries
% [sector_div, a_band_div] = function_divide_sector(sector_ab, a_band); 
% n_div = numel(sector_div); % number of divisions

% Horizontal trajectory
% lineseg = [WP(:,2),WP(:,3)]; 

% Divide trajectory in segments: level(0), climb(1), descend(-1)
segments = function_divide_trajectory(WP);

n_seg = numel(segments);

k = 1;

while (k<=n_seg)&&(~(p_in_found&&p_out_found))   
    if ~p_in_found
        
        [p_in_found, p_inaux, p_out_found, p_outaux] = function_find_p_in(segments(k).WP, segments(k).phase, sector_ab, a_band);
        if p_in_found
            p_in = p_inaux;
            
            if p_out_found
                p_out = p_outaux;
            else
                [p_out_found, p_outaux] = function_find_p_out(segments(k).WP, segments(k).phase, sector_ab, a_band);
                if p_out_found
                    p_out = p_outaux;
                else
                    k = k+1;
                end
            end
        else
            k = k+1;
        end
    else
        [p_out_found, p_outaux] = function_find_p_out(segments(k).WP, segments(k).phase, sector_ab, a_band);
        if p_out_found
            p_out = p_outaux; 
        else
            k = k+1;
        end
    end 
end

end

function [p_in_found, p_in, p_out_found, p_out] = function_find_p_in(WP, phase ,sector_ab, a_band)
% Is the entry point in segment?
% Compute entry point

p_in_found  = false;
p_out_found = false;
p_in  = [];
p_out = [];

lineseg = [WP(:,2),WP(:,3)]; % horizontal trajectory

switch phase
    case 0  % level
        h = WP(1,4); % Trajectory altitude
        ab_idx = ((a_band-5)<=h)&((a_band+5)>h);
        
        if any(ab_idx) % if the segment is in the altitude range (otherwise is_p_in = flase)
            sec_poly = sector_ab{ab_idx};
            [in,~] = intersect(sec_poly,lineseg);
            if ~isempty(in)
                p_in_found = true;
                % Two options: a) only p_in in segment b) both p_in and
                % p_out in segment
                in_end = inpolygon(WP(end,2),WP(end,3),sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                
                % Interpolate time
                
                lon_int = in([1,end],1);
                lat_int = in([1,end],2);
                
%                 F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,1));
%                 t_int = F(lon_int, lat_int);

                dist_vec = sqrt((WP(:,2)-WP(1,2)).^2+(WP(:,3)-WP(1,3)).^2);
                dist_int = sqrt((lon_int-WP(1,2)).^2+(lat_int-WP(1,3)).^2);
                t_int = interp1(dist_vec, WP(:,1), dist_int);
                
                % Identify entry and exit point
                [time, index_t] = sort(t_int);
                longitude       = lon_int(index_t);
                latitude        = lat_int(index_t);
                    
                if in_end % a) Only p_in in segment
                    p_in = [time(1), longitude(1), latitude(1), h];
                else % b) Both p_in and p_out in segment
                    p_in = [time(1), longitude(1), latitude(1), h];
                    p_out_found = true;
                    p_out = [time(2), longitude(2), latitude(2), h];
                end

            end
        end
        
    case 1  % climb
        
        h_ini = WP(1,4); % Initial trajectory altitude
        h_end = WP(end,4); % Final trajectory altitude
        
        if (h_end>=a_band(1))&&(h_ini<=a_band(end))
            
            % Sector division into different geometries
            [sector_div, a_band_div] = function_divide_sector(sector_ab, a_band);
            n_div = numel(sector_div); % number of divisions
            
            for k = 1:n_div
                if h_ini<a_band_div(k,2) 
                    sec_poly = sector_div{k};
                    [in,~] = intersect(sec_poly,lineseg);
                    if ~isempty(in)
                        
                        % Interpolate altitude and time
                        
                        lon_int = in([1,end],1);
                        lat_int = in([1,end],2);
                        
%                         F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,4));
%                         h_int = F(lon_int, lat_int);
                        
                        %                         F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,1));
                        %                         t_int = F(lon_int, lat_int);
                        dist_vec = sqrt((WP(:,2)-WP(1,2)).^2+(WP(:,3)-WP(1,3)).^2);
                        dist_int = sqrt((lon_int-WP(1,2)).^2+(lat_int-WP(1,3)).^2);
                        t_int = interp1(dist_vec, WP(:,1), dist_int);
                        h_int = interp1(dist_vec, WP(:,4), dist_int);
                        
                        % Identify entry and exit point
                        [time, index_t] = sort(t_int);
                        altitude        = h_int(index_t);
                        longitude       = lon_int(index_t);
                        latitude        = lat_int(index_t);
                        
                        % Two options: a) p_in in the edge b) p_in in the
                        % inside of the polygon
                        in_ini = inpolygon(WP(1,2),WP(1,3),sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                        
                        if in_ini
                            if (WP(1,4)<(a_band_div(k,1)-5)) % p_in in the inside of the polygon
                            h_bottom = a_band_div(k,1)-5; % bottom altitude of the sector division
                            int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_bottom);
                            time_int = int_values(1); lon_int = int_values(2); lat_int = int_values(3);
                            in_aux = inpolygon(lon_int,lat_int,sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));  
                            if in_aux
                                p_in_found = true;
                                p_in = [time_int, lon_int, lat_int, h_bottom];
                        
                                break % break for loop
                            end
                            end
                        else % p_in in the edge
                            
                            if (altitude(1)>=(a_band_div(k,1)-5))&&(altitude(1)<(a_band_div(k,2)+5))
                                p_in_found = true;
                                p_in = [time(1), longitude(1), latitude(1), altitude(1)];
                                
                                break
                                
                            elseif altitude(1)<(a_band_div(k,1)-5)
                                h_bottom = a_band_div(k,1)-5; % bottom altitude of the sector division
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_bottom);
                                time_int = int_values(1); lon_int = int_values(2); lat_int = int_values(3);
                                in_aux = inpolygon(lon_int,lat_int,sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                                if in_aux
                                    p_in_found = true;
                                    p_in = [time_int, lon_int, lat_int, h_bottom];
                                    
                                    break % break for loop
                                end              
                            end                           
                        end                
                    end
                end
            end         
        end
        
    case -1 % Descent
        h_ini = WP(1,4); % Initial trajectory altitude
        h_end = WP(end,4); % Final trajectory altitude
        
        if (h_end<=a_band(end))||(h_ini>=a_band(1))
            
            % Sector division into different geometries
            [sector_div, a_band_div] = function_divide_sector(sector_ab, a_band);
            n_div = numel(sector_div); % number of divisions
            
            for k = n_div:-1:1
                if h_ini>a_band_div(k,1) 
                    sec_poly = sector_div{k};
                    [in,~] = intersect(sec_poly,lineseg);
                    if ~isempty(in)
                        % Interpolate altitude and time
                        
                        lon_int = in([1,end],1);
                        lat_int = in([1,end],2);
                        
%                         F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,4));
%                         h_int = F(lon_int, lat_int);
%                         
%                         F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,1));
%                         t_int = F(lon_int, lat_int);
                        
                        dist_vec = sqrt((WP(:,2)-WP(1,2)).^2+(WP(:,3)-WP(1,3)).^2);
                        dist_int = sqrt((lon_int-WP(1,2)).^2+(lat_int-WP(1,3)).^2);
                        t_int = interp1(dist_vec, WP(:,1), dist_int);
                        h_int = interp1(dist_vec, WP(:,4), dist_int);
                        
                        % Identify entry and exit point
                        [time, index_t] = sort(t_int);
                        altitude        = h_int(index_t);
                        longitude       = lon_int(index_t);
                        latitude        = lat_int(index_t);
                        
                        % Two options: a) p_in in the edge b) p_in in the
                        % inside of the polygon
                        in_ini = inpolygon(WP(1,2),WP(1,3),sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                        
                        if in_ini
                            if (WP(1,4)>(a_band_div(k,2)+5)) % p_in in the inside of the polygon
                            h_top = a_band_div(k,2)+5; % top altitude of the sector division
                            int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_top);
                            time_int = int_values(1); lon_int = int_values(2); lat_int = int_values(3);
                            in_aux = inpolygon(lon_int,lat_int,sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));  
                            if in_aux
                                p_in_found = true;
                                p_in = [time_int, lon_int, lat_int, h_top];
                        
                                break % break for loop
                            end
                            end
                            
                        else % p_in in the edge
                            
                            if (altitude(1)>=(a_band_div(k,1)-5))&&(altitude(1)<(a_band_div(k,2)+5))
                                p_in_found = true;
                                p_in = [time(1), longitude(1), latitude(1), altitude(1)];
                                
                                break
                                
                            elseif altitude(1)>(a_band_div(k,2)+5)
                                h_top = a_band_div(k,1)-5; % top altitude of the sector division
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_top);
                                time_int = int_values(1); lon_int = int_values(2); lat_int = int_values(3);
                                in_aux = inpolygon(lon_int,lat_int,sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                                if in_aux
                                    p_in_found = true;
                                    p_in = [time_int, lon_int, lat_int, h_top];
                                    
                                    break % break for loop
                                end              
                            end                           
                        end 
                        
                        
                    end
                    
                end
                
            end
            
        end
end

end

function [p_out_found, p_out] = function_find_p_out(WP, phase, sector_ab, a_band)
% Is the exit point in segment?
% Compute exit point
p_out_found = false;
p_out = [];

lineseg = [WP(:,2),WP(:,3)]; % horizontal trajectory

switch phase
    case 0  % level
        h = WP(1,4); % Trajectory altitude
        ab_idx = ((a_band-5)<=h)&((a_band+5)>h);
        
        if any(ab_idx) % if the segment is in the altitude range (otherwise is_p_in = flase)
            sec_poly = sector_ab{ab_idx};
            [in,~] = intersect(sec_poly,lineseg);
            if ~isempty(in)
                
                % Only if the last point is outside the sector
                in_end = inpolygon(WP(end,2),WP(end,3),sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                if ~in_end
                    p_out_found = true;
                    
                    % Interpolate time
                    
                    lon_int = in([1,end],1);
                    lat_int = in([1,end],2);
                    
                    %                     F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,1));
                    %                     t_int = F(lon_int, lat_int);
                    
                    
                    dist_vec = sqrt((WP(:,2)-WP(1,2)).^2+(WP(:,3)-WP(1,3)).^2);
                    dist_int = sqrt((lon_int-WP(1,2)).^2+(lat_int-WP(1,3)).^2);
                    t_int = interp1(dist_vec, WP(:,1), dist_int);
                    
                    % Identify entry and exit point
                    [time, index_t] = sort(t_int);
                    longitude       = lon_int(index_t);
                    latitude        = lat_int(index_t);
                    
                    p_out = [time(2), longitude(2), latitude(2), h];
                end
            end
        end
    case 1  % Climb
        h_ini = WP(1,4); % Initial trajectory altitude
        h_end = WP(end,4); % Final trajectory altitude
        
        if (h_end>=a_band(1))&&(h_ini<=a_band(end))
            
            % Sector division into different geometries
            [sector_div, a_band_div] = function_divide_sector(sector_ab, a_band);
            n_div = numel(sector_div); % number of divisions
            
            for k = n_div:-1:1
                if h_end>a_band_div(k,1) 
                    sec_poly = sector_div{k};
                    [in,~] = intersect(sec_poly,lineseg);
                    if ~isempty(in)
                        % Interpolate altitude and time
                        
                        lon_int = in([1,end],1);
                        lat_int = in([1,end],2);
                        
%                         F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,4));
%                         h_int = F(lon_int, lat_int);
%                         
%                         F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,1));
%                         t_int = F(lon_int, lat_int);
                        
                        
                        dist_vec = sqrt((WP(:,2)-WP(1,2)).^2+(WP(:,3)-WP(1,3)).^2);
                        dist_int = sqrt((lon_int-WP(1,2)).^2+(lat_int-WP(1,3)).^2);
                        t_int = interp1(dist_vec, WP(:,1), dist_int);
                        h_int = interp1(dist_vec, WP(:,4), dist_int);
                        
                        % Identify entry and exit point
                        [time, index_t] = sort(t_int);
                        altitude        = h_int(index_t);
                        longitude       = lon_int(index_t);
                        latitude        = lat_int(index_t);
                        
                        % Two options: a) p_out in the edge b) p_out in the
                        % inside of the polygon
                        in_end = inpolygon(WP(end,2),WP(end,3),sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                        
                        if in_end
                            if(WP(end,4)>(a_band_div(k,2)+5)) % p_out in the inside of the polygon
                            h_top = a_band_div(k,2)+5; % top altitude of the sector division
                            int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_top);
                            time_int = int_values(1); lon_int = int_values(2); lat_int = int_values(3);
                            in_aux = inpolygon(lon_int,lat_int,sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));  
                            if in_aux
                                p_out_found = true;
                                p_out = [time_int, lon_int, lat_int, h_top];
                        
                                break % break for loop
                            end
                            end
                        else % p_in in the edge
                            
                            if (altitude(2)>=(a_band_div(k,1)-5))&&(altitude(2)<(a_band_div(k,2)+5))
                                p_out_found = true;
                                p_out = [time(2), longitude(2), latitude(2), altitude(2)];
                                
                                break
                                
                            elseif altitude(2)>(a_band_div(k,2)+5)
                                h_top = a_band_div(k,1)-5; % bottom altitude of the sector division
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_top);
                                time_int = int_values(1); lon_int = int_values(2); lat_int = int_values(3);
                                in_aux = inpolygon(lon_int,lat_int,sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                                if in_aux
                                    p_out_found = true;
                                    p_out = [time_int, lon_int, lat_int, h_top];
                                    
                                    break % break for loop
                                end              
                            end                           
                        end
                    end
                end
            end
            
        end
    case -1  % Descent
        h_ini = WP(1,4); % Initial trajectory altitude
        h_end = WP(end,4); % Final trajectory altitude
        
        if (h_end<=a_band(end))||(h_ini>=a_band(1))
            
            % Sector division into different geometries
            [sector_div, a_band_div] = function_divide_sector(sector_ab, a_band);
            n_div = numel(sector_div); % number of divisions
            
            for k = 1:n_div
                if h_end<a_band_div(k,2) 
                    sec_poly = sector_div{k};
                    [in,~] = intersect(sec_poly,lineseg);
                    if ~isempty(in)
                        
                        % Interpolate altitude and time
                        
                        lon_int = in([1,end],1);
                        lat_int = in([1,end],2);
                        
%                         F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,4));
%                         h_int = F(lon_int, lat_int);
%                         
%                         F = scatteredInterpolant(WP(:,2),WP(:,3),WP(:,1));
%                         t_int = F(lon_int, lat_int);
                        
                        dist_vec = sqrt((WP(:,2)-WP(1,2)).^2+(WP(:,3)-WP(1,3)).^2);
                        dist_int = sqrt((lon_int-WP(1,2)).^2+(lat_int-WP(1,3)).^2);
                        t_int = interp1(dist_vec, WP(:,1), dist_int);
                        h_int = interp1(dist_vec, WP(:,4), dist_int);
                        
                        % Identify entry and exit point
                        [time, index_t] = sort(t_int);
                        altitude        = h_int(index_t);
                        longitude       = lon_int(index_t);
                        latitude        = lat_int(index_t);
                        
                        % Two options: a) p_out in the edge b) p_out in the
                        % inside of the polygon
                        in_end = inpolygon(WP(end,2),WP(end,3),sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                        
                        if in_end
                            if(WP(end,4)<(a_band_div(k,1)-5)) % p_in in the inside of the polygon
                            h_bottom = a_band_div(k,1)-5; % bottom altitude of the sector division
                            int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_bottom);
                            time_int = int_values(1); lon_int = int_values(2); lat_int = int_values(3);
                            in_aux = inpolygon(lon_int,lat_int,sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));  
                            if in_aux
                                p_out_found = true;
                                p_out = [time_int, lon_int, lat_int, h_bottom];
                        
                                break % break for loop
                            end
                            end
                            
                        else % p_in in the edge
                            
                            if (altitude(2)>=(a_band_div(k,1)-5))&&(altitude(2)<(a_band_div(k,2)+5))
                                p_out_found = true;
                                p_out = [time(2), longitude(2), latitude(2), altitude(2)];
                                
                                break
                                
                            elseif altitude(2)<(a_band_div(k,1)-5)
                                h_bottom = a_band_div(k,1)-5; % bottom altitude of the sector division
                                int_values = interp1(WP(:,4), [WP(:,1), WP(:,2), WP(:,3)], h_bottom);
                                time_int = int_values(1); lon_int = int_values(2); lat_int = int_values(3);
                                in_aux = inpolygon(lon_int,lat_int,sec_poly.Vertices(:,1),sec_poly.Vertices(:,2));
                                if in_aux
                                    p_out_found = true;
                                    p_out = [time_int, lon_int, lat_int, h_bottom];
                                    
                                    break % break for loop
                                end              
                            end                           
                        end                
                    end
                end
            end    
            
        end
end


end

function segments = function_divide_trajectory(WP)

% This function divides the aircraft trajectory WP(t, lon, lat, h) into segments:
% 0: Level flight
% 1: Climb
% -1: Descent

[nwp,~] = size(WP); % number of waypoints

alt_dif = diff(WP(:,4));
nseg = numel(find(diff(sign(alt_dif))))+1; % Number of segments

segments(nseg) = struct;

i_vec = 1;
for i = 1:nseg
    type = sign(alt_dif(i_vec));
    segments(i).phase = type;
    i_ini = i_vec;
    while sign(alt_dif(i_vec))==type
        i_vec = i_vec+1;
        if i_vec == nwp
            break
        end
    end
    segments(i).WP = WP(i_ini:i_vec, :);
end

end