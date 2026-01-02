function [ASCR_k, M] = contra_function_ASCR_k(sector_ab, flows_j, cost_obstacles, Wij, a_band, adjacent_sectors)
% Computes ASCR for one sector as sum_i,sum_j AFCR(i,j) * Wij(i,j).
%
% If global FORCE_UNIT_WEIGHTS is true:
%   - AFCR is computed for ALL flows (to force plotting / diagnostics) using
%     contra_function_AFCR_j.
% Otherwise:
%   - AFCR is computed only where Wij>0 using contra_function_AFCR_j_onlyweigths.

global FORCE_UNIT_WEIGHTS

if isempty(FORCE_UNIT_WEIGHTS)
    FORCE_UNIT_WEIGHTS = false;
end

if isempty(Wij)
    M = [];
    ASCR_k = NaN;
    return
end

% Compute AFCR matrix
if FORCE_UNIT_WEIGHTS
    % Compute AFCR for all flows (independent of Wij), mainly for testing/plotting
    AFCR_i = contra_function_AFCR_j(cost_obstacles, sector_ab, a_band, flows_j, adjacent_sectors);
else
    % Compute AFCR only where Wij>0
    AFCR_i = contra_function_AFCR_j_onlyweigths(cost_obstacles, sector_ab, flows_j, Wij, a_band, adjacent_sectors);
end

if isempty(AFCR_i)
    M = [];
    ASCR_k = NaN;
    return
end

% Safety: AFCR_i and Wij must match dimensions to multiply elementwise
if ~isequal(size(AFCR_i), size(Wij))
    error("Size mismatch: AFCR_i is %s but Wij is %s.", mat2str(size(AFCR_i)), mat2str(size(Wij)));
end

if ~isempty(FORCE_UNIT_WEIGHTS) && FORCE_UNIT_WEIGHTS
    Wij_eff = ones(size(Wij));              % unit weights where weights are used
else
    Wij_eff = Wij;
end

M = AFCR_i .* Wij_eff;
ASCR_k = sum(M, 'all', 'omitnan');


end
