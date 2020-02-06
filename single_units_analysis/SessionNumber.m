gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);


channum = 1: length(data_file.good_data); 
all_filenames = cell(1,length(channum));
numcontacts = cell(length(size(data_file.good_data(1).channel_data.wf.waveForms)), length(channum));
filenamecont = cell(2,length(channum));
for i = 1:length(channum)
    filename = data_file.good_data(i).channel_data.filename;
    filename = erase(filename, "cinterocdrft_stab_fft_sigmat.mat") ;    
    all_filenames(i) = cellstr(filename(1:6));
    numcontacts(:,i) = cellstr(string(size(data_file.good_data(i).channel_data.wf.waveForms)));
    filenamecont(1,i) = all_filenames(i);
    filenamecont(2,i) = numcontacts(3,i);
end
filenamecont = str2double(filenamecont);
sessions = unique(filenamecont(1,:));
%Find common values between arrays (C) (C() is the same vector as sessions()) and corresponding indexes (ia). 
[C,ia,~] = intersect(filenamecont(1,:),sessions);

%only keep the first occurrences of each filename = each session, together with the electrode contact
%type
uniquefilenamecont = filenamecont(:,ia);
%determine the electrode contact types (3,24,32 or 24+32 =56)
conttypes = unique(uniquefilenamecont(2,:));
%determine the number of sessions with each electrode contact type
out = [conttypes;histc(uniquefilenamecont(2,:),conttypes)];

