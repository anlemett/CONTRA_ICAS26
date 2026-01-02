clearvars; clc;

warnState = warning;
warning('off', 'MATLAB:polyshape:repairedBySimplify');

set(0,'DefaultFigureVisible','on');   % make sure figures can appear

global PLOT_MINCUT PLOT_MINCUT_FL
PLOT_MINCUT    = true;
PLOT_MINCUT_FL = 370;

tic
[ASCR, sector_name] = contra_main_function();
toc

warning(warnState);  % restore warnings
