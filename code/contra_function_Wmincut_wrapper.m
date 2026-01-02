function Wmincut = contra_function_Wmincut_wrapper(sector_pgon, T, B, cost_polygons, h)
% Wrapper for mincut computation:
% - If global PLOT_MINCUT is true, it uses contra_function_Wmincut_draw
%   but ONLY for the selected altitude band PLOT_MINCUT_FL.
% - Otherwise, it uses contra_function_Wmincut.
%
% Inputs:
%   sector_pgon   : polyshape in lon/lat
%   T, B          : [N x 2] edges in lon/lat (columns: lon, lat)
%   cost_polygons : polyshape array (lon/lat) of obstacles (can be empty)
%   h             : altitude band (FL) currently being processed
%
% Globals:
%   PLOT_MINCUT    : logical
%   PLOT_MINCUT_FL : scalar FL (e.g., 370). If empty, plot for all altitudes.

global PLOT_MINCUT PLOT_MINCUT_FL

if isempty(PLOT_MINCUT)
    PLOT_MINCUT = false;
end

doPlot = false;
if PLOT_MINCUT
    if isempty(PLOT_MINCUT_FL)
        doPlot = true;        % plot for all altitude bands
    else
        doPlot = (h == PLOT_MINCUT_FL); % plot only for selected altitude band
    end
end

if doPlot
    Wmincut = contra_function_Wmincut_draw(sector_pgon, T, B, cost_polygons);
    drawnow; % ensure figures render immediately
else
    Wmincut = contra_function_Wmincut(sector_pgon, T, B, cost_polygons);
end

end
