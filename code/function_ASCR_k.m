function [ASCR_k, M] = function_ASCR_k(sector_ab, flows_j, weather_polygons, Wij, a_band)

% Compute ASCR (Available Sector Capacity Ratio) of sector k

AFCR_i = function_AFCR_j_onlyweigths(weather_polygons, sector_ab, flows_j, Wij, a_band);

% a) Most restrictive altitude
% AFCR = min(AFCR_i);
% ASCR_k = AFCR.*Wj;

% b) Altitude-dependent weights
M = AFCR_i.*Wij;
ASCR_k = sum(sum(M,'omitnan'));

end