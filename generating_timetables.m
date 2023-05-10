
uiwait(msgbox("Selecteer EDF bestand van Capnografie metingen","modal"));   % Dialog box for opening EDF file
[EDF_filename,EDF_path] = uigetfile('*.edf')   % Select EDF file

uiwait(msgbox("Selecteer GDT bestand van Ergospirometrie metingen","modal"));   % Dialog box for opening GDT file
[GDT_filename,GDT_path] = uigetfile('*GDT')    % Select GDT file

EDF_file = edfread(strcat(EDF_path,EDF_filename))


% ------------- Paste gdt file parser here ------------------
options = detectImportOptions(strcat(GDT_path,GDT_filename), 'FileType','text');
T = readtable(strcat(GDT_path,GDT_filename), options); 
T.PaCO2(isnan(T.PaCO2)) = 0; %replace Nan values with 0 so they dont get removed in next step
gdtfile = T(~any(ismissing(T),2),:);%remove all rows that have missing value or NaN. now we only have the table remaining.
gdtSamplingTime = seconds(duration(gdtfile.Time, 'InputFormat','mm:ss'));%duration of recording
gdtSamplePeriod = round(mean(diff(gdtSamplingTime)));%difference between the samples and we take the mean of that
gdtSampleRate = 1/gdtSamplePeriod;
gdtfile.Time=[];%remove Time collumn
Time=transpose(seconds(0:gdtSamplePeriod:(height(gdtfile)-1)*gdtSamplePeriod ));%create new time collumn
gdtTimeTable=table2timetable(gdtfile,'RowTimes',Time);

%GDT_file = ones(97,12, 'uint8'); % example of gdt timetable, delete this line if code is pasted in
