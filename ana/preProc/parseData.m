function dat = parseData(header,paths,varargin)

% dat = parseData(header,paths,varargin)
%
% Builds data structure for all projects.
%
% HEADER file-specific and project-specific information
% PATHS project path structure containing at minimum paths.rare field
% VARARGIN options string specifying what fields to create.
% valid VARARGIN values
% 'k' use spike-times & event-markers from Dr. Koyano
%
% last modified 6-28-16
% apj

option                                          = lower(varargin);

% choose event-markers from Dr. Koyano
if ~isempty(strfind(option{1},'k'));
    ML_Data                                     = 1;
else
    ML_Data                                     = 0;
end

if ~isempty(strfind(option{1},'m'));
    ML_Tank                                     = 1;
else
    ML_Tank                                     = 0;
end

% build DAT
dat.h                                           = header;
dat.h.ver                                       = 1.05;
dat.h.parseDate                                 = datestr(now);

% set format of spike files in DAT
dat.h.fileType                                  = header.src.spike_format;

% pre-allocate in DAT
dat.c                                           = [];
dat.s                                           = [];

%% load spike-times from kenji
if ML_Data
    
    % load corresponding neurostruct file
    filename                                    = fullfile(paths.raw,dat.h.filename,...
        'SD_NeuroStruct',['SD' num2str(header.sdNum) '.mat']);
    load(filename)
    
    % get cell list
    ML_chans                                    = cell2mat(NeuroStruct(1).cells(:,1));
    
    % get position of block in data structure
    for i = 1:length(NeuroStruct)
        if ML_Data
            %             regexp(NeuroStruct(i).block,dat.h.src.kenji)
            blockNum                            = i;
            break
        end
    end
    
end

%% load header (TDT-mat) file
if ML_Tank
    header_name                                 = [dat.h.filename header.src.block '-tdt' header.suffix.load '.mat'];
    header_file                                 = fullfile(paths.rare,header_name);
else
    header_name                                 = [dat.h.filename header.src.block '-tdt.mat'];
    header_file                                 = fullfile(paths.rare,header_name);
end

if exist(header_file, 'file')
    load(header_file); % load tdt struct
    fprintf('%s\n',['Loaded tank file for header: ' header_name])
    % else
    %     error('%s\n',['Can''t find any tank files for header: ' header_name]);
    % end
    
    dat.h.date                                  = tdt.dateS;
    % dat.h.date                                  = datestr(tdt.dateN);
    if ML_Tank
        date_vector                             = datevec(dat.h.date,'yyyy-mmm-dd');
    else
        date_vector                             = datevec(dat.h.date,'dd-mmm-yyyy');
    end
    % define stimulus LUT
    [ind, set, file]                            = stimLUT(header.proj,date_vector);
    dat.h.stim.all                              = ind;
    dat.h.stim.img                              = file;
    dat.h.stim.set                              = set;
    % define digital signals
    if ML_Tank
        timestamp.time                          = tdt.times*1000;    % save in msec
    else
        timestamp.time                          = tdt.times;    % save in msec
    end
    timestamp.code                              = tdt.codes;
    %                                     tdt.duration
    
    %% decode events
    
    if ML_Tank
        stimID_shift                            = 100;
        % add monkey logic event codes
        dat.h.stim.ML                           = header.stim.ML;
    else
        stimID_shift                            = 1000;
        
        %% figure out these codes
        elevens                                 = find(timestamp.code==11); % remove 11's- why? outdated?
        if ~isempty(elevens)
            keyboard
            %         timestamp.code(elevens)             = [];
            %         timestamp.time(elevens)             = [];
        end
    end
    
    % shift stim id to correspond with event codes
    stimID                                      = dat.h.stim.all + stimID_shift;
    
    %% find event-markers
    if ML_Tank
        [ts_ind, ~]                             = ismember(timestamp.code, stimID);
        trial_stim                              = find(ts_ind);
        % % trial_started                           = find(timestamp.code==header.stim.ML.isi_code);
        % % trial_end                               = find(timestamp.code==header.stim.ML.stim_code);
        trial_started                           = find(timestamp.code==header.stim.ML.stim_code);
        trial_end                               = find(timestamp.code==header.stim.ML.isi_code);
        
        % select completed trials
        % % wholeTrials                             = intersect(trial_stim-1, trial_started);
        % % GoodTrials                              = intersect(wholeTrials, trial_end-2);
        wholeTrials                             = intersect(trial_stim+1, trial_started);
        GoodTrials                              = intersect(wholeTrials, trial_end-1);
    else
        GoodTrials                              = find(ismember(tdt.codes(1:end-1),stimID) ...
            & tdt.codes(2:end)==21);
    end
    
    % unpack header to current workspace
    eval(unpack_header(dat));
    
    %% get spikes
    abc                                         = 'a':'z';
    name_tick                                   = 0;
    
    % channel loop
    channels                                    = dat.h.chan;
    for c = 1:length(channels)
        
        if ML_Data
            % create file name
            tmpTitle                                = ['*' dat.h.src.tank ...
                dat.h.src.block '_' num2str(channels(c)) '_*.' header.src.spike_format];
            
            %% select channel in pre-sorted file from Kenji
            chan_ind                            = find(ML_chans==channels(c));
            
            if ~isempty(chan_ind)
                spike_times                     = NeuroStruct(blockNum).cells{chan_ind(1),3};
                
                spike_nums                      = ones(size(spike_times));
                tsc.cluster_class               = [spike_nums spike_times];
            else
                tsc                             = [];
                fprintf('%s\n',[tmpTitle ' not found in Neurostruct data.'])
            end
            
        else
            % create file name
            tmpTitle                                = ['times_' dat.h.src.tank ...
                'a_' num2str(channels(c)) '_*.' header.src.spike_format];
            
            %% load waveClus spike files
            tmpDir                              = dir(fullfile(paths.rare, tmpTitle));
            spike_time_file                     = fullfile(paths.rare, tmpDir.name);
            try
                tsc                             = load(spike_time_file);
                fprintf(1,'.');  % indicate progress
            catch
                tsc                             = [];
                fprintf('%s\n',[tmpTitle ' not read.'])
            end
        end
        
        %% save spike-times for dat struct
        if ~isempty(tsc)
            %         Nspikes                    = length(unique(tsc.cluster_class(:,1)));
            %         for s = 1:Nspikes
            s = 1; % until I figure out how to add second neurons from the same channel
            % define channel
            name_tick                           = name_tick+1;
            TS{name_tick}                       = tsc.cluster_class(:,2); % in milliseconds
            
            % create spike-name
            snames(name_tick,:)     = ['sig' sprintf('%03d',channels(c)) abc(s)];
            %         end
        else
            fprintf('%s\n',['Spike-file empty for channel# ' num2str(c) ...
                ' in tank: ' dat.h.src.tank])
        end
    end
    
    % check quality
    if length(GoodTrials)<=1
        error('No GoodTrials found')
    end
    
    try
        dat.h.snames                            = snames;       % names of spikes contained in dat
    catch
        error(['No spikes in /rare for: ' dat.h.src.tank])
    end
    
    %% create dat array
    dat.s                                       = nan(length(GoodTrials),SPIKE_DIMS,length(TS)); % Good_crossing
    
    for trialNum = 1:length(GoodTrials)
        
        dat.c(trialNum,SESSION)                 = dat.h.session;  % session number
        dat.c(trialNum,DATE)                    = datenum(date_vector(1,1:3));  % session date - yearMonthDay
        dat.c(trialNum,TR)                      = trialNum;  % trial number
        % %     stim                                    = timestamp.code(GoodTrials(trialNum)+1) - stimID_shift;  % stimulus id
        stim                                    = timestamp.code(GoodTrials(trialNum)+(ML_Tank*-1)) - stimID_shift;  % stimulus id
        
        % if kenji data - decode accordingly
        if ML_Tank
            stim                                = dat.h.stim.ML.decode(stim);
        end
        
        % set stimulus descriptors
        if header.proj==1 % category selectivity
            
            [category]                          = decodeToroidCategorySelectivity(stim,dat.h);
            dat.c(trialNum,CATEGORY)            = category;
            
        elseif header.proj==2 % face-learning project
            
            dat.c(trialNum,STIM)                = stim;
            [type, face, step]                  = decodeFaceLearning(stim,dat.h);
            
            dat.c(trialNum,TYPE)                = type;
            dat.c(trialNum,FACE)                = face;
            dat.c(trialNum,STEP)                = step;
            
        elseif header.proj==3 % anti-car project
            
            dat.c(trialNum,STIM)                = stim;
            [type, face, step]                  = decodeRadTan(stim,dat.h);
            
            dat.c(trialNum,TYPE)                = type;
            dat.c(trialNum,FACE)                = face;
            dat.c(trialNum,STEP)                = step;
            
        end
        
        dat.c(trialNum,TREL)                    = timestamp.time(GoodTrials(trialNum));
        dat.c(trialNum,T_STIMON)                = timestamp.time(GoodTrials(trialNum));
        dat.c(trialNum,T_STIMOFF)               = timestamp.time(GoodTrials(trialNum)+1);
        dat.c(trialNum,T_PRE)                   = dat.c(trialNum,T_STIMON)-300;
        dat.c(trialNum,T_POST)                  = dat.c(trialNum,T_STIMOFF)+300;
        
        % look for spikes, relative to each stimulus
        for chan = 1:length(TS)
            Tstart                              = dat.c(trialNum,T_PRE);
            Tstop                               = dat.c(trialNum,T_POST);
            spikes                              = find(TS{chan}>=Tstart&TS{chan}<=Tstop);
            if ~isempty(spikes)
                dat.s(trialNum,1:length(spikes),chan) = TS{chan}(spikes);
            end
        end
        
        if mod(trialNum,100)==0
            fprintf('%d\n',trialNum)
        end
        
        if trialNum==1
            fprintf('\n');
        end
        
    end
    
    %% adjust times to be relative to the stim onset
    
    %% ?
    dat.s                                       = dat.s(:,1:SPIKE_DIMS,:); % unnecessary?
    
    % define relative time
    dat.h.Taligned                              = TREL;
    Talign                                      = dat.c(:,dat.h.Taligned);
    
    % shift event timestamps
    T_align_events                              = repmat(Talign,[1 length(dat.h.timestamps)]);
    dat.c(:,dat.h.timestamps)                   = dat.c(:,dat.h.timestamps)-T_align_events;
    
    % shift spike timestamps
    [~,c,s]                                     = size(dat.s);
    T_align_spikes                              = repmat(Talign,[1,c,s]);
    dat.s                                       = dat.s - T_align_spikes;
    
end

end

