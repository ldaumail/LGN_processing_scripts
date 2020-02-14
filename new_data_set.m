gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
newdatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';


oldchannelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
old_data_file = load(oldchannelfilename);

%Reject 160517, K cell, as no well triggered (1)
%Reject 180806 p1 uclust17, M cell, as doesn't seem well triggered (46)
%Reject 181207 (B) uclust22, M cell, as doesn't seem well triggered (55)
new_data = old_data_file.good_data([2:45,47:54,56:end]);

allfilename = [newdatadir 'refined_dataset'];
save(strcat(allfilename, '.mat'), 'new_data', '-v7.3');
 