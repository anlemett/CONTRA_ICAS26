function Wij_m1 = contra_function_Wij_ini(flows_j)

nonEmpty = ~cellfun(@isempty, flows_j);

if ~any(nonEmpty)
    Wij_m1 = 0;
    return
end

aux = [flows_j{nonEmpty}];
aux2 = vertcat(aux.triplet);
all_triplets = unique(aux2, 'rows');

[ntri, ~] = size(all_triplets);
nab = numel(flows_j);

Wij_m1 = ones(nab, ntri);
Wij_m1(~nonEmpty, :) = 0;

for i = find(nonEmpty)'
    trip_i = vertcat(flows_j{i}.triplet);
    for j = 1:ntri
        if ~ismember(all_triplets(j,:), trip_i, 'rows')
            Wij_m1(i,j) = 0;
        end
    end
end

s = sum(Wij_m1, 'all');
if s > 0
    Wij_m1 = Wij_m1 / s;
else
    % Fallback: uniform over altitude bands that have flows
    Wij_m1 = double(nonEmpty);
    Wij_m1 = Wij_m1 / sum(Wij_m1);
end

end
