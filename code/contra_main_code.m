clearvars; clc;

warnState = warning;
warning('off', 'MATLAB:polyshape:repairedBySimplify');

set(0,'DefaultFigureVisible','on');   % make sure figures can appear

global PLOT_MINCUT PLOT_MINCUT_FL FORCE_UNIT_WEIGHTS
PLOT_MINCUT    = true;
PLOT_MINCUT_FL = 370;
FORCE_UNIT_WEIGHTS = false; % if true - plots for all flows appear

tic
[ASCR, sector_name] = contra_main_function();
toc

warning(warnState);  % restore warnings
