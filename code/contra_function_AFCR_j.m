function AFCR_i = contra_function_AFCR_j(cost_obstacles, sector_ab, a_band, flows_j, adjacent_sectors)
% Compute AFCR_i for each altitude band i and each triplet (flow).
% Expects cost_obstacles as a cell array of structs with field .pgon (polyshape).

% Robustly gather all triplets across non-empty altitude bands
nonEmpty = ~cellfun(@isempty, flows_j);
if ~any(nonEmpty)
    AFCR_i = nan(numel(a_band), 0);
    return
end

aux = [flows_j{nonEmpty}];
aux2 = vertcat(aux.triplet);
all_triplets = unique(aux2, 'rows');
ntri = size(all_triplets, 1);

nab = numel(a_band);
AFCR_i = nan(nab, ntri);

% Prepare obstacle polyshapes (can be empty)
obs_poly = polyshape.empty(0,1);
if ~isempty(cost_obstacles)
    obs = [cost_obstacles{:}];
    obs_poly = vertcat(obs.pgon);
end

for i = 1:nab

    h = a_band(i);

    if isempty(flows_j{i})
        continue
    end

    nj = numel(flows_j{i});
    if mod(nj, 2) ~= 0
        error("flows_j{%d} must contain an even number of flows (forward + inverse).", i);
    end

    for j = 1:(nj/2)

        Wmincut = contra_function_Wmincut_wrapper(sector_ab{i}, flows_j{i}(j).T, flows_j{i}(j).B, obs_poly, h, adjacent_sectors);

        % Prefer stored baseline if present, otherwise compute it
        if isfield(flows_j{i}(j), 'Omincut') && ~isempty(flows_j{i}(j).Omincut)
            Omincut = flows_j{i}(j).Omincut;
        else
            Omincut = contra_function_Omincut(sector_ab{i}, flows_j{i}(j).T, flows_j{i}(j).B);
        end

        [~, j1] = ismember(flows_j{i}(j).triplet, all_triplets, 'rows');
        [~, j2] = ismember(flows_j{i}(j + nj/2).triplet, all_triplets, 'rows');

        if j1 == 0 || j2 == 0
            continue
        end

        if isempty(Omincut) || ~isfinite(Omincut) || Omincut <= 0
            AFCR_i(i, j1) = NaN;
            AFCR_i(i, j2) = NaN;
        else
            AFCR_i(i, j1) = Wmincut ./ Omincut;
            AFCR_i(i, j2) = AFCR_i(i, j1);
        end
    end
end

end
