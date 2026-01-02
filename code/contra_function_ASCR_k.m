function [ASCR_k, M] = contra_function_ASCR_k(sector_ab, flows_j, cost_obstacles, Wij, a_band)

AFCR_i = contra_function_AFCR_j_onlyweigths(cost_obstacles, sector_ab, flows_j, Wij, a_band);

if isempty(AFCR_i) || isempty(Wij)
    M = [];
    ASCR_k = NaN;
    return
end

if ~isequal(size(AFCR_i), size(Wij))
    error("Size mismatch: AFCR_i is %s but Wij is %s.", mat2str(size(AFCR_i)), mat2str(size(Wij)));
end

M = AFCR_i .* Wij;
ASCR_k = sum(M, 'all', 'omitnan');

end
