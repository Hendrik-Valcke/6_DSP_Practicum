
uiwait(msgbox("Selecteer EDF bestand van Capnografie metingen","modal"));   % Dialog box for opening EDF file
[EDF_filename,EDF_path] = uigetfile('*.edf')   % Select EDF file

uiwait(msgbox("Selecteer GDT bestand van Ergospirometrie metingen","modal"));   % Dialog box for opening GDT file
[GDT_filename,GDT_path] = uigetfile('*GDT')    % Select GDT file

EDF_file = edfread(strcat(EDF_path,EDF_filename))


% ------------- Paste gdt file parser here ------------------

GDT_file = ones(97,12, 'uint8'); % example of gdt timetable, delete this line if code is pasted in
