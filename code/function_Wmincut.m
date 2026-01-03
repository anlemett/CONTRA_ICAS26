
function Wmincut = function_Wmincut(sector_pgon, T, B, weather_polygons)

global nowcast

% Wmincut = zeros(size(weather_polygons));
% [nT,nM] = size(weather_polygons); % nT: forecast times - nM: total members
        
% sector_lon = [S(:,1);T(:,1);D(:,1);B(:,1)]; 
% sector_lat = [S(:,2);T(:,2);D(:,2);B(:,2)];
% sector_pgon = polyshape(sector_lon, sector_lat);

[lon_c, lat_c] = centroid(sector_pgon);
[sector_y, sector_x] = function_spherical_to_eq_azimuth(sector_pgon.Vertices(:,2), sector_pgon.Vertices(:,1), lat_c, lon_c);
sector_pgon_xy = polyshape(sector_x, sector_y);

[T_y, T_x] = function_spherical_to_eq_azimuth(T(:,2), T(:,1), lat_c, lon_c);
[B_y, B_x] = function_spherical_to_eq_azimuth(B(:,2), B(:,1), lat_c, lon_c);

% figure, hold on
% plot(sector_pgon_xy)
% plot(B_x, B_y, 'Linewidth', 2, 'Color', 'r', 'Marker', 'x')
% plot(T_x, T_y, 'Linewidth', 2, 'Color', 'b', 'Marker', 'o')

% name_alph = 'a':'z';

% for t = 1:nT
%     for m = 1:nM
        
        nw = numel(weather_polygons); % number of weather polygons
        
        % Select weather cells within the sector
        count = 0;
        weather_in = [];
        
        for i = 1:nw
            
            polyout = intersect(sector_pgon, weather_polygons(i));
            
            if polyout.NumRegions>0 
                
                w_vertices = weather_polygons(i).Vertices;
                [y, x] = function_spherical_to_eq_azimuth(w_vertices(:,2), w_vertices(:,1), lat_c, lon_c);
                w_pgon_xy = polyshape(x,y);
                polyout = intersect(sector_pgon_xy, w_pgon_xy);
                
                if polyout.NumRegions>0 
                    weather_in = [weather_in, i];
                    count = count+1;
                    x_w{count} = polyout.Vertices(:,1);
                    y_w{count} = polyout.Vertices(:,2);
                end
%                 plot(polyout, 'EdgeColor', [199, 0, 57 ]/255, 'FaceColor', [199, 0, 57 ]/255, 'FaceAlpha', 0.2)
  
            end
        end
        
        % Create graph
        NW = length(weather_in); % Number of weather cells
        NG = 2 + NW; % Number of nodes in the graph
        
        % 1) nodes and edges
%         nodenames = ['T', graph_names, 'B'];
        nodes_comb = nchoosek(1:NG,2);
        sg = nodes_comb(:,1)'; tg = nodes_comb(:,2)';
        
        % 2) weights
        weights = zeros(size(sg));
        
        [weights(NG-1), x1, x2] = boundary_to_boundary_fun(B_x, B_y, T_x, T_y);
%         plot([x1(1), x2(1)], [x1(2), x2(2)], 'b--')

        for i = 1:NW
            
            x = x_w{i};
            y = y_w{i};
            
            if isempty(polyxpoly(T_x, T_y, x, y))
                [weights((sg==1)&(tg==(i+1))), x1, x2] = boundary_to_boundary_fun(T_x, T_y, x, y); % T
%                 plot([x1(1), x2(1)], [x1(2), x2(2)],'b--')
                if ~nowcast
                    % 13.5 NM = 25002 meters
                    weights((sg==1)&(tg==(i+1))) = weights((sg==1)&(tg==(i+1))) - 25002;
                    if weights((sg==1)&(tg==(i+1))) < 0
                        weights((sg==1)&(tg==(i+1))) = 0;
                    end
                end
            end
            if isempty(polyxpoly(B_x, B_y, x, y))
                [weights((sg==(i+1))&(tg==NG)), x1, x2] = boundary_to_boundary_fun(B_x, B_y, x, y); % B
%                 plot([x1(1), x2(1)], [x1(2), x2(2)], 'b--')
                if ~nowcast
                    % 13.5 NM = 25002 meters
                    weights((sg==(i+1))&(tg==NG)) = weights((sg==(i+1))&(tg==NG)) - 25002;
                    if weights((sg==(i+1))&(tg==NG)) < 0
                        weights((sg==(i+1))&(tg==NG)) = 0;
                    end
                end
            end
            
            % for each weather cell
            
            for j = i+1:NW
                x2 = x_w{j};
                y2 = y_w{j};
                
                [weights((sg==(i+1))&(tg==(j+1))), x1, x2] = boundary_to_boundary_fun(x2, y2, x, y);
%                 plot([x1(1), x2(1)], [x1(2), x2(2)],'b:')
                if ~nowcast
                    % 27 NM = 50004 meters
                    weights((sg==(i+1))&(tg==(j+1))) = weights((sg==(i+1))&(tg==(j+1))) - 50004;
                    if weights((sg==(i+1))&(tg==(j+1))) < 0
                        weights((sg==(i+1))&(tg==(j+1))) = 0;
                    end
                end
            end
            
        end

        G = graph(sg,tg,weights); % Create graph
        [short_path, Wmincut] = shortestpath(G, 1, NG); % Compute shortest path
        
%         figure
%         H = plot(G,'EdgeLabel',G.Edges.Weight,'Layout','layered');
%         highlight(H,short_path,'EdgeColor','r')
% %     end
% end

end


function [d, p1, p2] = boundary_to_boundary_fun(x1, y1, x2, y2)

% Function that computes the distance between two boundaries 1 and 2

% x1, y1: boundary 1
% x2, y2: boundary 2

% d: distance between two boundaties
% p1, p2, points over the boundaries where they are closer together

N1 = length(x1); %length of boundary 1
N2 = length(x2); %length of boundary 2

d = inf;

if N1 == 1
    
    pt = [x1, y1];
    p1 = pt;
    %         v2 = [x1(i+1), y1(i+1)];
    for j = 1:N2-1
        w1 = [x2(j), y2(j)];
        w2 = [x2(j+1), y2(j+1)];
        
        %             [d_new, p1_new, p2_new] = line_to_line(v1, v2, w1, w2);
        [d_new, x] = point_to_line(pt, w1, w2);
        if d_new<d
            d = d_new;
            p2 = x;
        end
    end
    
    
elseif N2 == 1
    
    pt = [x2, y2];
    p2 = pt;
    
    for i = 1:N1-1
        v1 = [x1(i), y1(i)];
        v2 = [x1(i+1), y1(i+1)];
        
        %             [d_new, p1_new, p2_new] = line_to_line(v1, v2, w1, w2);
        [d_new, x] = point_to_line(pt, v1, v2);
        
        if d_new<d
            d = d_new;
            p1 = x;
        end
        
    end
    
else
    
    for i = 1:N1-1
        v1 = [x1(i), y1(i)];
        v2 = [x1(i+1), y1(i+1)];
        for j = 1:N2-1
            w1 = [x2(j), y2(j)];
            w2 = [x2(j+1), y2(j+1)];
            
            [d_new, p1_new, p2_new] = line_to_line(v1, v2, w1, w2);
            if d_new<d
                d = d_new;
                p1 = p1_new;
                p2 = p2_new;
            end
        end
    end
    
end

end

function [d, x1, x2] = line_to_line(v1, v2, w1, w2)

% Function that computes the distance between 2 segments
% v1, v2, w1, w2 are the three-dimensional coordinates of the segments

% x1, x2, points over the lines where the two are closer together

x1v = zeros(4,2); x1v(1,:) = v1; x1v(2,:) = v2;
x2v = zeros(4,2); x2v(3,:) = w1; x2v(4,:) = w2;

[dv(1), x2v(1,:)] = point_to_line(v1, w1, w2);
[dv(2), x2v(2,:)] = point_to_line(v2, w1, w2);
[dv(3), x1v(3,:)] = point_to_line(w1, v1, v2);
[dv(4), x1v(4,:)] = point_to_line(w2, v1, v2);

[d, idx] = min(dv);

x1 = x1v(idx,:);
x2 = x2v(idx,:);

end

function [d, x] = point_to_line(pt, v1, v2)

% Function that computes the distance between a point and a segment
% pt, v1, and v2 are the three-dimensional coordinates of the point,
% one vertex on the segment, and a second vertex on the segment, respectively

% x: point that is closer to pt

    a = pt - v1;
    b = v2 - v1;
    dotp = a*b';
    len = b*b';

    param = dotp/len;

    if param<0
        x = v1; 
    elseif param>1
        x = v2;
    else
        x = v1 + param*b;
    end

    d = sqrt((x-pt)*(x-pt)');
    
end