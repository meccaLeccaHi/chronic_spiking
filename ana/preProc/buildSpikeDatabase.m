function snames = buildSpikeDatabase(header,paths)
%
% Usage: snames = buildSpikeDatabase(header,paths)
% Script for making datasets that contain a single spike across multiple days.
%
% HEADER file-specific and project-specific information
% PATHS project path structure containing at minimum paths.rare field
%
% last modified 6-17-16
% apj

% read experiment descriptors
monkey                      = header{1}.filename(1:end-8);
% block                       = header{1}.src.block;
exp_code                    = [monkey(1) paths.code];

% get list of dat-struct files for each day
search_name                 = [header{1}.filename header{1}.src.block header{1}.suffix.load '.mat'];
% search_name                = [monkey '*' block header{1}.suffix '.mat'];
dat_files                   = dir(fullfile(paths.mas,search_name));

% display error if empty
if isempty(dat_files)
    error(['No spikes files found for ' search_name])
end

%% create spike-list
% loop through each file in header
SS                          = [];
for n = 1:length(header)
    
    % load file
    load(fullfile(paths.mas,[header{n}.src.tank header{n}.src.block header{n}.suffix.load '.mat']));
    
    SS                      = [SS; dat.h.snames];
    
    % get number of spikes
    nSpikes                 = size(SS,1);
    
    % loop through spikes to get list
    for n = 1:nSpikes
        
        if n>nSpikes
            break
        else
            
            % find instances of each spike in SS
            str             = SS(n,:);
            reps            = strmatch(str,SS);
            
            % remove redundant instances
            if length(reps)>1
                SS(reps(2:end),:) = [];
            end
            
            nSpikes         = size(SS,1);
            
        end
    end
end

% sort spike list by cell and print
snames                      = sortrows(SS,[4 5 6]);
disp(snames)


%% build new file for each spike
% spike loop
for s = 1:size(snames,1)

    day                     = cell(length(header),1);
    observedIn              = [];
        
    % loop through days and build dat array
    for n = 1:length(header);

        % load each file
        foo                 = load(fullfile(paths.mas,[header{n}.src.tank header{n}.src.block header{n}.suffix.load '.mat']));
        
        if ~isempty(foo.dat.c)
            % add to cell array
            day{n}          = foo.dat;
            % select only particular spike
            day{n}          = select_spikename(day{n},snames(s,:));
            
            %         %%?
            %         day{n}.s            = day{n}.s(:,1:80);
        end
    end
    
    % delete empty cells
    day                     = day(~cellfun('isempty',day));
    
    % concatenate dat array elements
    dat                     = day{1};
    for n = 2:length(day)
        if ~isempty(day{n})
            dat             = catdat(dat,day{n});
        end
    end
    
    % set filename and path
%     filename                = [header{n}.filename header{n}.src.block header{n}.suffix '.mat'];
    filename                = [exp_code '_' snames(s,:) header{1}.suffix.save '.mat'];
    dat.h.filename          = filename;
    dat.h.observedIn        = observedIn;
    savefile                = fullfile(paths.mas,filename);
    
    keyboard
    
    %% save file
    save(savefile,'dat');
    disp(['Built spikefile: ' savefile]);
end
