gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels\';
data_file = load(channelfilename);

keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
     31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 for i = 1:length(keepidx)
 data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(600:1901,:,:));
 
 filename = [data_file.good_data(i).channel_data.filename, f{2}];
 filename(strfind(filename, 'mat')) = [];
 filename(strfind(filename, '.')) = []; 
 
 all_locs = nan(4,length(data(1,:)));
 all_pks = nan(4,length(data(1,:)));
 all_lpsu= nan(length(data(:,1)),length(data(1,:)));
 nyq = 15000;
 %all_locs =nan(4,length(channel_data(1,1,:)));
 for n =1:length(data(1,:))
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpsu      = filtfilt(bwb,bwa, data(:,n));
   %all_lpsu(:,n) = lpsu;
   %plot(lpsu)
   locs = findpeaks(lpsu(30:1150));
    if length(locs.loc) >= 4
   all_pks(:,n) = lpsu(locs.loc(1:4));
   end 
 end
 
 channelfilename = [channeldir filename];
 save(strcat(channelfilename, '.mat'), 'all_pks');

 end