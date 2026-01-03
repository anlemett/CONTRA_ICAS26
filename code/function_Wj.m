function [Wj, Wij, total_ac] = function_Wj(t_ini, t_fin, p_in, p_out, a_band, flows_j, AC_in)

% All flows j in sector k:
aux = [flows_j{:}]; aux2 = vertcat(aux.triplet); all_triplets = unique(aux2, 'rows');
[ntri, ~] = size(all_triplets); % number of triplets

nab = numel(a_band); % number of altitude bands
Wij = zeros(nab, ntri);

% Aircraft that fly sector k in the time interval [t_ini, t_fin]
ac_idx = (t_ini<=p_out(:,1))&(t_fin>=p_in(:,1));

total_ac = sum(ac_idx); % Total number of aircraft that fly over sector k 

p_in  = p_in(ac_idx,:);
p_out = p_out(ac_idx, :);
AC_in = AC_in(ac_idx);

% figure, hold on
%

for a = 1:total_ac
    [a_band_idx, tri] = function_find_ac(p_in(a,:), p_out(a,:), a_band, flows_j);
    if ~isempty(tri)
        tri_idx = find(ismember(all_triplets,tri,'rows'));
        Wij(a_band_idx, tri_idx) = Wij(a_band_idx, tri_idx) + 1;
        
    end
    
%    plot([p_in(a, 2), p_out(a, 2)], [p_in(a, 3), p_out(a, 3)])
%    plot(p_in(a, 2), p_in(a, 3), '*')
%    plot(p_out(a, 2), p_out(a, 3), 'o')

%    plot3(p_in(a,2), p_in(a,3), p_in(a,4), 'r*')
%    plot3(p_out(a,2), p_out(a,3),p_out(a,4), 'bo')
%    plot3(AC_in(a).WP(:,2),AC_in(a).WP(:,3),AC_in(a).WP(:,4))
%    plot(AC_in(a).WP(:,2),AC_in(a).WP(:,3))
end

total_ac = sum(sum(Wij));

Wj = sum(Wij)/total_ac;
Wij = Wij/total_ac;

end

function [a_band_idx, tri] = function_find_ac(p_in, p_out, a_band, flows_j)

a_band_idx = [];
tri = [];

min_lon = 0;
max_lon = 30;
% min_lat = 35;
% max_lat = 65;

% This function identifies the altitude band and the triplet associated
% with the aircraft

t_in = p_in(1); t_out = p_out(1);
lon_in = p_in(2); lon_out = p_out(2);
lat_in = p_in(3); lat_out = p_out(3);
h_in = p_in(4); h_out = p_out(4);

h = (h_in+h_out)/2;
% Find altitude band
a_band_idx = find((a_band>=(h-5))&(a_band<=(h+5)), 1, 'first');
% Find corresponding triplet
flow = flows_j{a_band_idx}; 

if ~isempty(flow)
all_triplets = vertcat(flow.triplet);
[ntri, ~] = size(all_triplets); % number of triplets

coefficients = polyfit([lon_in, lon_out], [lat_in, lat_out], 1);
a = coefficients (1);
b = coefficients (2);

x = [min_lon, max_lon]; y=a*x+b;

for j = 1:(ntri/2)
    [xS,yS] = polyxpoly(x,y, flow(j).S(:,1), flow(j).S(:,2)); % Intersect source edge
    if ~isempty(xS)
        [xD,yD] = polyxpoly(x,y, flow(j).D(:,1), flow(j).D(:,2)); % Intersect destination edge
        if ~isempty(xD)
            if norm([xS,yS]-[lon_in, lat_in])<norm([xD,yD]-[lon_in, lat_in])
                tri_idx = j; 
            else
                tri_idx = j + (ntri/2);
            end
            tri = flow(tri_idx).triplet;
            break
        end
    end
end

% tri = flow(tri_idx).triplet;

end

end
