clear; clc; clear function_Wmincut_draw;
warning('off');

global PLOT_MINCUT PLOT_MINCUT_FL TARGET_FL
PLOT_MINCUT    = false;
PLOT_MINCUT_FL = 370;
TARGET_FL        = 370;

tic
[ASCR, sector_names, sector_time, sector_data] = function_main();

disp(array2table(ASCR{1}(:).', 'VariableNames', cellstr(string(sector_time)), 'RowNames', cellstr(string(sector_names))))

toc