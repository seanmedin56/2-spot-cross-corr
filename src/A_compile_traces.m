% Nick's script for compiling pipeline data
% Script to Compile Data Sets and Find Stripe Centers
close all
clear 
%--------------------------Set Path Specs, ID Vars------------------------%
FolderPath = '../dat/steven_2colorelongation';
% FolderPath = 'D:/Data/Nick/LivemRNA/LivemRNAFISH/Dropbox (Garcia Lab)/mHMM/weka/';
project = 'steven_2colorelongation'; %Project Identifier
%%% folders
fig_path = ['../fig/experimental_system/' project '/preprocessing/'];
data_path = ['../dat/' project '/']; % data mat directory
%%% fig subfolders
ap_pos_path = [fig_path 'ap_positioning/'];
fluo_path = [fig_path 'fluo_stats/'];
%%% assign save names
trace_name = [data_path 'raw_traces_01' project]; % names for compiled trace struct
nucleus_name = [data_path 'ellipse_info_01' project]; % names for compiled elipse struct

%%% make filepaths
mkdir(data_path);
mkdir([ap_pos_path '/stripe_fits']);
mkdir(fluo_path);
%%% cleaning params
keyword = 'elongation'; % Keyword to ensure only sets from current project are pulled
include_vec = [4]; %data set numbers to include
% show_ap_fit_figs = 0;
snippet_size = 15; % particles within snippet/2+1 are at risk for tracking issues
% pre_post_padding = 10; % max mun frames for which nucleus can be MIA at start or end
%--------------------------Obtain Relevant Filepaths----------------------%
% store set names
dirinfo = dir(FolderPath);
dirinfo(~[dirinfo.isdir]) = []; %remove non-directories
cp_filenames = {}; % particles
set_nums = [];
for d = 1 : length(dirinfo)
    thisdir = dirinfo(d).name;
    % Skip files lacking project keyword 
    if isempty(strfind(thisdir,keyword)) 
        continue
    end
    % Remove sets not in include_vec
    set_num_start_ind = strfind(thisdir,'_');
    set_num_start_ind = set_num_start_ind(end);
    set_num = str2num(thisdir(set_num_start_ind+1:end));    
    if sum(set_num==include_vec) ~= 1 
        continue
    end    
    set_nums = [set_nums set_num];
    % append file paths
    cp_filenames = [cp_filenames {[thisdir '/CompiledParticles.mat']}];            
    fov_filenames = [fov_filenames {[thisdir '/FrameInfo.mat']}];           
end

trace_struct = struct; % Generate data structure to store extracted trace sets

%%% compile traces and nuclei from across experimental sets
j_pass = 0; % indexing variable to keep track of iterations
total_matched = 0;
for i = 1:length(cp_filenames) % Loop through filenames    
    % read in raw files
    load([FolderPath fov_filenames{i}]) % FrameInfo Info
    yDim = FrameInfo(1).LinesPerFrame;
    xDim = FrameInfo(1).PixelsPerLine;
    raw_data = load([FolderPath cp_filenames{i}]); % Particles    
    setID = set_nums(i);    
    
    % pull trace and nuclei variables
    time_raw = raw_data.ElapsedTime*60; % time vector            
    traces_raw1 = raw_data.AllTracesVector{1}; % Array with a column for each trace
    traces_raw2 = raw_data.AllTracesVector{2};
    frames_raw = 1:length(time_raw); % Frame list    
    first_frame = raw_data.nc14; % Get frame that marks start of nc14
    last_frame = frames_raw(end); % last frame in set
    % filter trace mat and time
    traces_clean1 = traces_raw1(first_frame:end,:);
    traces_clean2 = traces_raw2(first_frame:end,:);
    time_clean = time_raw(first_frame:end);    
    time_clean = time_clean - min(time_clean); % Normalize to start of nc14
    frames_clean = frames_raw(first_frame:end);    

    %%% Now Particles        
    fn = cp_filenames{i}; % Get filename to store in struct           
    % iterate through traces
    particles = raw_data.CompiledParticles; % extract particle set
    j_init = j_pass;
    for j = 1:size(traces_clean,2)        
        raw_trace = traces_clean(:,j);     
        trace_start = find(~isnan(raw_trace),1);
        trace_stop = find(~isnan(raw_trace),1,'last');        
        %Creat versions with all intervening frames present (missing frames
        %appear as NaNs
        trace_full = raw_trace(trace_start:trace_stop)';               
        time_full = time_clean(trace_start:trace_stop);                
        frames_full = frames_clean(trace_start:trace_stop);
        % skip small fragments
        short_flag = 0;
        if length(raw_trace(~isnan(raw_trace))) < 3            
            continue
        elseif length(raw_trace(~isnan(raw_trace))) < 5            
            short_flag = 1;
        end
        j_pass = j_pass + 1;
        
        % Pull variables from particle structure        
        fov_xPos_raw = particles(j).xPos;
        fov_yPos_raw = particles(j).yPos;
        fluo3_raw = particles(j).Fluo3;
        fluo5_raw = particles(j).Fluo5;
        pt_frames = particles(j).Frame;
        schnitz = particles(j).schnitz;
        % Only take values from CP struct corresponding to frames in
        % filtered frames
        NaN_vec = NaN(1,length(frames_full));             
        xPos = fov_xPos_raw(ismember(pt_frames,frames_full));        
        yPos = fov_yPos_raw(ismember(pt_frames,frames_full));
        fluo3 = NaN_vec;
        fluo3(ismember(frames_full,pt_frames)) = fluo3_raw(ismember(pt_frames,frames_full));
        fluo5 = NaN_vec;
        fluo5(ismember(frames_full,pt_frames)) = fluo5_raw(ismember(pt_frames,frames_full));
        pt_frames = pt_frames(ismember(pt_frames,frames_full));
        % look for edge issues
        edge_frames = (((xDim-xPos) <= 1.5*snippet_size/2)|(xPos <= 1.5*snippet_size/2)|...
                          ((yDim-yPos) <= 1.5*snippet_size/2)|(yPos <= 1.5*snippet_size/2))>0;
        trace_struct(j_pass).edge_flag = max(edge_frames)&&~short_flag;        
        % Record info in trace struct
        trace_struct(j_pass).cp_frames = pt_frames; % Keep track of frame correspondence
        trace_struct(j_pass).all_frames = frames_full;
        trace_struct(j_pass).nc14 = first_frame;
        trace_struct(j_pass).last_frame = last_frame;        
        trace_struct(j_pass).xPos = xPos;
        trace_struct(j_pass).yPos = yPos;
        trace_struct(j_pass).fluo = trace_full;
        trace_struct(j_pass).fluo3 = fluo3;
        trace_struct(j_pass).fluo5 = fluo5;
        trace_struct(j_pass).time = time_full;
        trace_struct(j_pass).FluoError = particles(j).FluoError; %Estimated error in bkg subtraction     

        % Identifier variables        
        particle = particles(j).OriginalParticle;        
        pID = eval([num2str(setID) '.' sprintf('%04d',particle)]);                    
        trace_struct(j_pass).ParticleID = pID;
        trace_struct(j_pass).setID = setID;
        trace_struct(j_pass).source_path = fn;                        
    end      
end

%%% Look for trace fragments that belong together. Stitch them up
%%% This is premised on the trace start times being sorted in ascending
%%% order!

set_index = [trace_struct.setID];
set_vec = unique(set_index);
% stitch together overlaps
remove_indices = [];
match_indices = [];
dupe_indices = [];
for s = 1:length(set_vec)
    set_trace_struct = trace_struct(set_index==set_vec(s));     
    base_indices = find(set_index==set_vec(s));
    % vectors used to asses fragment proximity
    start_t = [];
    stop_t = [];
    start_x = [];
    stop_x = [];
    start_y = [];
    stop_y = [];
    for i = 1:length(set_trace_struct)
        time = set_trace_struct(i).time;
        xPos = set_trace_struct(i).xPos;
        yPos = set_trace_struct(i).yPos;
        % add start and stop info
        start_t = [start_t time(1)];
        stop_t = [stop_t time(end)];
        start_x = [start_x xPos(1)];
        stop_x = [stop_x xPos(end)];
        start_y = [start_y yPos(1)];
        stop_y = [stop_y yPos(end)];
    end
    %%% enforce ascending sort order
    [start_t, si] = sort(start_t);
    stop_t = stop_t(si);
    start_x = start_x(si);
    stop_x = stop_x(si);
    start_y = start_y(si);
    stop_y = stop_y(si);
    base_indices = base_indices(si);   
    t_mat = repmat(start_t,length(start_t),1) - repmat(stop_t',1,length(stop_t));
    x_mat = repmat(start_x,length(start_x),1) - repmat(stop_x',1,length(stop_x));
    y_mat = repmat(start_y,length(start_y),1) - repmat(stop_y',1,length(stop_y));
    logic_mat = (sqrt((x_mat.^2 + y_mat.^2))<=5)&(t_mat>0)&(t_mat<120); % look for spatially and temporally proximate fragments
    logic_mat(eye(length(stop_y))==1) = 0; % remove diagonals
    overlap_ids = find(logic_mat) ;    
    tr_col = floor((overlap_ids-1)/length(stop_x))+1; % convert from linear to 2D
    tr_row = overlap_ids - (tr_col-1)*length(stop_x);
    tr_col = base_indices(tr_col);
    tr_row = base_indices(tr_row);        
    if length(unique(tr_col))~=length(tr_col) || length(unique(tr_row))~=length(tr_row)           
        warning('Duplicate Trace Fragments Detected. Removing.')
        col_mat = repmat(tr_col,length(tr_col),1)==repmat(tr_col',1,length(tr_col));
        row_mat = repmat(tr_row,length(tr_row),1)==repmat(tr_row',1,length(tr_row));        
        row_mat(eye(size(row_mat,1))==1)=0;
        col_mat(eye(size(col_mat,1))==1)=0;
        [~, row_ind] = find(max(row_mat));
        [~, col_ind] = find(max(col_mat));
        rm_vec = 1:length(tr_col);        
        dupe_indices = [dupe_indices tr_col(row_ind) tr_row(col_ind)];        
        tr_row = tr_row(~ismember(rm_vec,[row_ind col_ind]));
        tr_col = tr_col(~ismember(rm_vec,[row_ind col_ind]));                
    end    
    cat_fields = {'fluo','fluo3','fluo5','time','xPos','yPos'...
                    'cp_frames','all_frames'};
    for j = 1:length(tr_col)
        ID1 = min(tr_col(j),tr_row(j)); % keep earlier as base
        ID2 = max(tr_col(j),tr_row(j));
        if ismember(ID1,remove_indices)            
            ID1 = match_indices(ID1==remove_indices);            
        end
        % take ID info from earlier trace (smaller ID)
        base = trace_struct(ID1);
        extra = trace_struct(ID2);        
        for f = 1:length(cat_fields)
            base.(cat_fields{f}) = [base.(cat_fields{f}) extra.(cat_fields{f})];
        end
        base.edge_flag = max(base.edge_flag,extra.edge_flag);        
        % assign nucleus that is better fit
%         xp1 = nucleus_struct(nc_id_vec == base.ncID).xPos;
%         yp1 = nucleus_struct(nc_id_vec == base.ncID).yPos;
%         t1 = nucleus_struct(nc_id_vec == base.ncID).time;
%         xp2 = nucleus_struct(nc_id_vec == extra.ncID).xPos;
%         yp2 = nucleus_struct(nc_id_vec == extra.ncID).yPos;
%         t2 = nucleus_struct(nc_id_vec == extra.ncID).time;
%         f1 = ismember(round(t1),round(base.time));
%         if sum(f1)==length(base.time)             
%             r1 = sqrt((xp1(f1)-base.xPos).^2+(yp1(f1)-base.yPos).^2);
%         else
%             r1 = Inf;
%         end
%         f2 = ismember(round(t2),round(base.time));
%         if sum(f2) < length(base.time)   
%             error('asfa')
%         end
        % Add particle info to nucleus struct
        fn = fieldnames(trace_struct);
        for f = 1:length(fn)
            trace_struct(ID1).(fn{f}) = base.(fn{f});
        end        
        remove_indices = [remove_indices ID2];
        match_indices = [match_indices ID1];
    end        
end

% remove extra entries
index_vector = 1:length(trace_struct);
tr_particles = [trace_struct.ParticleID];
trace_struct = trace_struct(~ismember(index_vector,[remove_indices dupe_indices]));
rm_particles = tr_particles([remove_indices dupe_indices]);


%%

for i = 1:length(include_vec)   
    set_trace_struct = trace_struct([trace_struct.setID] == include_vec(i));
    new_trace_struct = [new_trace_struct set_trace_struct]; 
    disp(['Completed ' num2str(i) ' of ' num2str(length(include_vec))])
end         

save([trace_name '.mat'],'new_trace_struct') 
%%
