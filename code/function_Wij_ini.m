% Initiate variable Wij_m1

function Wij_m1 = function_Wij_ini(flows_j)

aux = [flows_j{:}]; aux2 = vertcat(aux.triplet); all_triplets = unique(aux2, 'rows');
[ntri, ~] = size(all_triplets); % number of triplets
nab = numel(flows_j); % number of altitude bands
Wij_m1 = ones(nab, ntri);

Wij_m1(cellfun(@isempty,flows_j),:) = 0; % 0 in altitudes with no possible flows (no adjacent sectors)

for i = find(~cellfun(@isempty,flows_j))'
    trip_j = vertcat(flows_j{i}.triplet);
    for j = 1:ntri
        [~, j1] = ismember(all_triplets(j,:), trip_j, 'rows');
        if j1 == 0
            Wij_m1(i,j) = 0;
        end
    end
end

Wij_m1 = Wij_m1/(sum(sum(Wij_m1)));

end