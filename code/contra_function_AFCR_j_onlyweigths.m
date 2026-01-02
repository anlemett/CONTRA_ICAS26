function AFCR_i = contra_function_AFCR_j_onlyweigths(cost_obstacles, sector_ab, flows_j, Wij, a_band)
% Compute AFCR_i for each altitude band i and triplet j.
% Expects cost_obstacles as a cell array of structs with field .pgon (polyshape).

% Robust collection of all triplets (avoid empty cells)
nonEmpty = ~cellfun(@isempty, flows_j);
if ~any(nonEmpty)
    AFCR_i = nan(size(Wij));
    return
end

aux = [flows_j{nonEmpty}];
aux2 = vertcat(aux.triplet);
all_triplets = unique(aux2, 'rows');

AFCR_i = nan(size(Wij));

vec_ab = find(sum(Wij, 2) > 0);

for i = vec_ab'

    h = a_band(i);

    vec_trip = find(Wij(i,:) > 0);
    trip_i = vertcat(flows_j{i}.triplet);

    if ~isempty(cost_obstacles)
        obs = [cost_obstacles{:}];
        obs_poly = vertcat(obs.pgon);

        for j = vec_trip
            [~, j1] = ismember(all_triplets(j,:), trip_i, 'rows');
            if j1 == 0
                continue
            end

            Wmincut = contra_function_Wmincut_wrapper(sector_ab{i}, flows_j{i}(j1).T, flows_j{i}(j1).B, obs_poly, h);
            Omincut = flows_j{i}(j1).Omincut;

            if isempty(Omincut) || ~isfinite(Omincut) || Omincut <= 0
                AFCR_i(i,j) = NaN;
            else
                AFCR_i(i,j) = Wmincut ./ Omincut;
            end
        end
    else
        AFCR_i(i, vec_trip) = 1;
    end
end

end
