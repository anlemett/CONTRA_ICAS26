function AFCR_i = function_AFCR_j_onlyweigths(weather_polygons, sector_ab, flows_j, Wij, a_band)

% [sector_ab, a_band, flows_j] = function_flows_sector_j(k, main_sectors, adjacent_sectors);

% All flows j in sector k:
aux = [flows_j{:}]; aux2 = vertcat(aux.triplet); all_triplets = unique(aux2, 'rows');
% [ntri, ~] = size(all_triplets); % number of triplets

% [nab,ntri] = size(Wij); 

AFCR_i = nan(size(Wij));


vec_ab = find(sum(Wij, 2)>0);
for i = vec_ab'
    
    h = a_band(i);
    vec_trip = find(Wij(i,:)>0);
    trip_j = vertcat(flows_j{i}.triplet);
    
    if numel(weather_polygons)>0
        aux = [weather_polygons{:}];
        w_poly = vertcat(aux.pgon);
        index = false(1,numel(weather_polygons));
    
        for j = vec_trip
            [~, j1] = ismember(all_triplets(j,:), trip_j, 'rows');

            %Wmincut = function_Wmincut(sector_ab{i}, flows_j{i}(j1).T, flows_j{i}(j1).B, w_poly); % Mincut at altitude band i with weather areas
            Wmincut = function_Wmincut_wrapper(sector_ab{i}, flows_j{i}(j1).T, flows_j{i}(j1).B, w_poly, h); % Mincut at altitude band i with weather areas
%             Omincut = function_Omincut(sector_ab{i}, flows_j{i}(j1).T, flows_j{i}(j1).B); % Mincut at altitude band i without weather areas
            Omincut = flows_j{i}(j1).Omincut;
            
            AFCR_i(i,j) = Wmincut./Omincut;
        end
    else
        for j = vec_trip
            AFCR_i(i,j) = 1;
        end
    end
end

end