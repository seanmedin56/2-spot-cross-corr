FolderPath = '../dat/steven_2colorelongation';
project = 'steven_2colorelongation'; %Project Identifier

trace_struct_final = struct;

%loads compiled particles from every folder
dirinfo = dir(FolderPath);
dirinfo(~[dirinfo.isdir]) = []; %remove non-directories
load([FolderPath filesep 
