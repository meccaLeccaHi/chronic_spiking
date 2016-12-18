function waveClus_shell
%
% Script for executing the entire data-production "pipeline".
%
% last modified 6-28-16
% apj

addpath('/home/lab/Cloud3/projects/reference_code')
addpath('/home/lab/Cloud3/projects/ana/preProc')

% set experimental variables
runSpikes                                           = 0;
parseSpikes                                         = 0;
analyzeMonkeys                                      = 1;
analyzePop                                          = 1;
spikeForm                                           = 'mat';
keepOriginal                                        = 'o';
getUnique                                           = 'u';
suffix                                              = '_test';
sdNum                                               = 40; % for parseData -- try 40, to start

% set paths
proj                                                = 3;
[paths,monkeys]                                     = setFilePath(proj);

monkeys([4]) = [];

for mo = 3:length(monkeys)
    
    % use original data for toroid
    if ~isempty(strfind('tr',monkeys{mo}(1)))
        tankVers                                    = 'a';
        dataVers                                    = 'a';
    else
        tankVers                                    = 'm';
        dataVers                                    = 'k';
    end
            dataVers                                    = 'a';

    
    % load header info
    header                                          = fileList(proj,monkeys{mo},spikeForm, ...
                                        [dataVers 's' keepOriginal],0,suffix,sdNum);
    
    dates                                           = [];
    for i = 1:length(header)
        dates                                       = [dates; header{i}.filename];
    end
    [~,date_idx]                                    = unique(dates,'rows');
    
    %% extract spike-times from analog data
    
    if runSpikes
        
        fprintf('%s\n','Extracting spike-times from analog data')
        
        % generate filenames
        ind                                         = 0;
        for i = 1:length(header)
            
            % create directory, if none exists
            tankName                                = fullfile(paths.scratch,header{i}.expt,header{i}.src.tank);
            if ~exist(tankName,'dir');
                mkdir(tankName);
            end
            
            % create file-list
            channels                                = header{i}.chan;
            for ii = 1:length(channels)
                file_name                           = [header{i}.filename header{i}.src.spike '_' num2str(channels(ii))];
                full_file                           = fullfile(tankName,file_name);
                if exist(full_file,'file')
                    ind                             = ind + 1;
                    files{ind}                      = [full_file '.mat'];
                    clust_files{ind}                = [full_file '_spikes.mat'];
                else
                    fprintf('%s\n',['Channel not found: ' [file_name '.mat']]);
                end
            end
            [~,indA]                                = unique(files);
            files                                   = files(indA);
            clust_files                             = clust_files(indA);
            
            % get spikes
            parfor k = 1:length(files)
                Get_spikes(files(k));
            end
            
            % get clusters
            parfor k = 1:length(files)
                Do_clustering(clust_files(k),paths.rare);
            end
        end
        
        %     % cleanup any leftover temp data
        %     waveClus_cleanup(paths);
    end
    
    %% parse spike-time data
    
    if parseSpikes
        
        fprintf('%s\n','Parsing spike-time data')
        
                % write file for each recording date (all neurons)
                buildDatabase(header,paths,[tankVers dataVers])
        
                % write file for each neuron (all dates)
                for d = 1:length(date_idx)
                    %             temp_header                             = fileList(proj,monkeys{mo},spikeForm,[dataVers keepOriginal],[],suffix,sdNum);
                    buildDailyDatabase(header{d},paths)
                end
        
        spike_header                                = fileList(proj,monkeys{mo},spikeForm,[dataVers keepOriginal getUnique],0,suffix,sdNum);
        buildSpikeDatabase(spike_header,paths)
        
    end
    
    %% analyze spike data
    
    if analyzeMonkeys
        if proj==2
            %             % plot raw responses
            %             insp_dat_fine(paths,monkeys{mo},header{mo}.suffix)
            %             % make dodecahedrons with embedded rasters
            %             rastDodec(paths,monkeys{mo},header{mo}.suffix)
            % make dodecahedrons
            insp_dat_identTraj(paths,monkeys{mo},header{mo}.suffix);
        elseif proj==3
            % plot raw responses
            insp_dat_rtTraj(paths,monkeys{mo},header{mo}.suffix);
        end
    end
end

%% population analyses

if analyzePop
    if proj==2
        insp_dat_identTrajPop(paths,monkeys,header{1}.suffix)
    elseif proj==3
        insp_dat_rtPop(paths,monkeys,header{1}.suffix)
    end
end

end


%%
function waveClus_cleanup(paths)
%
% Function to cleanup any leftover temp data created during waveClus sort process
%
% last modified 12-09-15
% apj

close all

dump_tags                   = ['exe' 'ata' 'run' '500' '_01' 'mag' 'lab'];
% ship_tags                 = ['times'];
dir_list                    = dir;
for i = 1:length(dir_list)
    if dir_list(i).isdir==0;  % if it's a directory, ignore
        if ~strncmp(dir_list(i).name,'.',1);  % if hidden, ignore
            if ~isempty(regexp(dir_list(i).name,'fig2print','once'))
                figName = dir_list(i).name;
                outName = fullfile(paths.sortFigs,[dir_list(i).name(1:end-3) 'png']);
                hT = imread(figName);
                imwrite(hT,outName)
                clear hT;
                %                 movefile(dir_list(i).name,fullfile(paths.sortFigs,dir_list(i).name))
                
            elseif ~isempty(regexp(dump_tags,dir_list(i).name((end-2):end),'once'))...
                    ||~isempty(regexp('spikes',dir_list(i).name((end-9):(end-4)),'once'))
                % if it matches with dump_tags,delete it
                if ~exist('dibs1','var')
                    fprintf('%s\n','Cleaning up files...');
                    dibs1       = 1;
                end
                delete(dir_list(i).name)
                fprintf('%s\n',[dir_list(i).name ' deleted'])
                %             elseif ~isempty(regexp(ship_tags,dir_list(i).name(1:5),'once'))
                %                 % if it matches with ship_tags,move it
                %                 if ~exist('dibs2','var')
                %                     fprintf('%s\n','Moving spike files...');
                %                     dibs2       = 1;
                %                 end
                %                 movefile(dir_list(i).name,fullfile(paths.rare,dir_list(i).name));
                %                 fprintf('%s\n',[dir_list(i).name ' moved to /rare'])
            end
        end
    end
end
end