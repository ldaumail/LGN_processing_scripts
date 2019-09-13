
function [ppnev] = create_ppnev()
dir = 'C:\Users\maier\Documents\MATLAB\KiloSortUtils-Master/offlineBRAutoSort';
directory  = 'C:\Users\maier\Documents\LGNinfo_4LD-20190826T172747Z-001\LGN_19019_B_cinterocdrft002_data\';
BRdatafile = '190119_B_cinterocdrft002';
filename   = [directory BRdatafile]; 
ext = 'ns6';
addpath(genpath(directory))
addpath(genpath(dir))

ppnev = offlineBRAutoSort(strcat(filename,'.', ext));
save(strcat(filename, '.ppnev'))
end