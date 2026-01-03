function AFCR_i = function_AFCR_j(weather_polygons, sector_ab, a_band, flows_j)

% [sector_ab, a_band, flows_j] = function_flows_sector_j(k, main_sectors, adjacent_sectors);

% All flows j in sector k:
aux = [flows_j{:}]; aux2 = vertcat(aux.triplet); all_triplets = unique(aux2, 'rows');
[ntri, ~] = size(all_triplets); % number of triplets

[nab,~] = size(a_band); % number of altitude bands
AFCR_i = nan(nab, ntri);

for i = 1:nab
    
    [~, nj] = size(flows_j{i});
    for j = 1:(nj/2)
            Wmincut = function_Wmincut(sector_ab{i}, flows_j{i}(j).T, flows_j{i}(j).B, weather_polygons); % Mincut at altitude band i with weather areas
            Omincut = function_Omincut(sector_ab{i}, flows_j{i}(j).T, flows_j{i}(j).B); % Mincut at altitude band i without weather areas

            [~, j1] = ismember(flows_j{i}(j).triplet, all_triplets, 'rows');
            [~, j2] = ismember(flows_j{i}(j+nj/2).triplet, all_triplets, 'rows');
            
            AFCR_i(i,j1) = Wmincut./Omincut;
            AFCR_i(i,j2) = AFCR_i(i,j1); % inverse flows have the same mincut
    end
    
end

end