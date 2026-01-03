
function Omincut = function_Omincut(sector_pgon, T, B)

% [nT,nM] = size(weather_polygons); % nT: forecast times - nM: total members
        
% sector_lon = [S(:,1);T(:,1);D(:,1);B(:,1)]; 
% sector_lat = [S(:,2);T(:,2);D(:,2);B(:,2)];
% sector_pgon = polyshape(sector_lon, sector_lat);

[lon_c, lat_c] = centroid(sector_pgon);

[T_y, T_x] = function_spherical_to_eq_azimuth(T(:,2), T(:,1), lat_c, lon_c);
[B_y, B_x] = function_spherical_to_eq_azimuth(B(:,2), B(:,1), lat_c, lon_c);

% figure, hold on
% plot(sector_pgon_xy)
% plot(B_x, B_y, 'Linewidth', 2, 'Color', 'r')
% plot(T_x, T_y, 'Linewidth', 2, 'Color', 'r')

% name_alph = 'a':'z';


weights = boundary_to_boundary_fun(B_x, B_y, T_x, T_y);

G = graph(1,2,weights); % Create graph
[~, short_path] = shortestpath(G, 1, 2); % Compute shortest path

% Omincut = short_path*ones(size(weather_polygons));
Omincut = short_path;
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