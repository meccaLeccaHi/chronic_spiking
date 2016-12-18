function insp_dat_identTraj(paths,monkey,suffix)
%   loop through relevant files in 'paths.mas'
%
% Usage: insp_dat_identTraj(paths,monkey,suffix)
% MONKEY is the monkey's name
% PATHS is the project path structure containing at minimum paths.rare field
% SUFFIX is a string appended to the end of each filename
%
% last modified 6-24-16
% apj

% setup
figsOn                                                          = 0; % switch for figures (on vs. off)
visTog                                                          = 'off';
Bkgwin                                                          = 300; % length of prestimulus period (msec)
Plotwin                                                         = [-300 600]; % plot display size (msec)
Rspwin                                                          = [100 400]; % time window for evoked rate
prs                                                             = 150; % print resolution
% fullscr                                                         = get(0,'ScreenSize'); % get screen size
recomb                                                          = {1;2:4;5:6;7;8;9}; % values for averaging
dodecPlot_cutoff                                                = .5; % minimum rate necessary to plot dodec
pop_low_crit                                                    = 2.5; % minimum rate (low) for pop. analysis
pop_hi_crit                                                     = 5; % minimum rate (high) for pop. analysis

% set face numbers
humans                                                          = 1:12;
monkeys                                                         = 101:112;
faces                                                           = [humans monkeys];
% objects                                                         = 201:212;

dodecStr                                                        = fullfile(paths.results,'dodecs'); % directory for saving plots
monkeyDataDir                                                   = fullfile(paths.results,'fine',monkey); % directory for saving plots

% get file list
exp_code                                                        = [monkey(1) paths.code];
dirList                                                         = dir(fullfile(paths.mas,[exp_code '*' suffix.load '.mat']));

% pre-allocate population variables
% humPopResp = cell(length(dirList),1);
humPopRecomb = cell(length(dirList),1);
humPopNormResp = cell(length(dirList),1);
humPopRankResp = cell(length(dirList),1);
% monPopResp = cell(length(dirList),1);
monPopRecomb = cell(length(dirList),1);
monPopNormResp = cell(length(dirList),1);
monPopRankResp = cell(length(dirList),1);
% objPopResp = cell(length(dirList),1);
% cellDateArray = cell(length(dirList),1);
inclCrit_rateLim                                                = zeros(length(dirList),1);
% inclCritRateLimLo = zeros(length(dirList),1);
% inclCritRateLimHi = zeros(length(dirList),1);
% inclCritRateT = zeros(length(dirList),1);
% inclCritRateD = zeros(length(dirList),1);
% inclCritRateL = zeros(length(dirList),1);

sNamesList                                                      = [];
% testStr1 = [];

% loop through file list
for i = 1:length(dirList);
    close all
    
    %load file
    filename                                                    = fullfile(paths.mas,dirList(i).name);
    load(filename)
    disp(['Loaded: ' dirList(i).name])
    unpack
    
    if strfind('t r',monkey(1))
        dates                                                   = datenum(dat.h.date,'dd-mmm-yy');
    else
        dates                                                   = datenum(dat.h.date,'yyyy-mmm-dd');
    end
    
    morphs                                                      = dat.h.stim.morphSteps;
    nMorphs                                                     = length(morphs);
    morphTraj                                                   = [0 morphs];
    sNamesList                                                  = [sNamesList; dat.h.snames];
    
    % interpolate responses to new length
    recombMorphs                                                = nan(1,length(recomb));
    for k = 1:length(recomb)
        recombMorphs(k)                                         = mean(morphTraj(recomb{k}));
    end
    
    % preallocate
    tetPile                                                     = nan(12,nMorphs+1); % for tetrahedron plot
    tetPile_imList                                              = nan(length(faces),nMorphs+1);
    tetPile_recomb                                              = nan(length(faces)/2,length(recomb));
    
    tcArray                                                     = nan(1,nMorphs+1); % for tuning curves
    %     stdArray = nan(1,nMorphs+1); % preallocate for tuning curves standard deviation
    %     steArray = nan(1,nMorphs+1); % preallocate for tuning curves standard deviation
    
    %     backRate = []; stimRate = []; rateT = [];
    
    inclCrit_rateLim_face                                       = nan(length(faces),1);
    
    % step through faces
    for f = 1:length(faces)
        %         f = 16;
        %         mI = 0;
        
        % define norm for current type
        if ismember(faces(f),humans)
            typeNum                                             = 1;
            %             normStim = 193;
            fNum                                                = 0;
        elseif ismember(faces(f),monkeys)
            typeNum                                             = 2;
            %             normStim = 194;
            fNum                                                = 100;
        end
        
        % look at norm response before stepping through morphs
        normdat                                                 = select_trials(dat,DATE,dates,TYPE,typeNum,FACE,fNum);
        imNum                                                   = unique(normdat.c(:,STIM)); % stimulus image index
        tetPile_imList(f,1)                                     = imNum(1);
        
        %         % loop through first seven days
        %         normDayArray = cell(1,numDays);
        %         for n = 1:numDays
        %             dayDat = select_trials(dat,FACE,fNum,TYPE,typeNum,DATE,dates(n));
        %
        %             [r,~] = find(dayDat.s>Rspwin(1)&dayDat.s<Rspwin(2));
        %             trialSpikes = histc(r,1:length(dayDat.s(:,1)));
        %             [r,~] = find(dayDat.s<0);
        %             trialBackGround = histc(r,1:length(dayDat.s(:,1)));
        %             eResp = mean(trialSpikes)*(1000/diff(Rspwin));
        %             bResp = mean(trialBackGround)*(1000/Bkgwin);
        %             normDayArray{n} = [normDayArray{n};eResp-bResp];
        %         end
        %         morphDateArray = normDayArray;
        
        % add to cumulative sdf matrix
        if ~isempty(normdat.s);
            [r,~]                                               = find(normdat.s>Rspwin(1)&normdat.s<Rspwin(2));
            trialSpikes                                         = histc(r,1:length(normdat.s(:,1)));
            [r,~]                                               = find(normdat.s<0);
            trialBackGround                                     = histc(r,1:length(normdat.s(:,1)));
            
            eResp                                               = mean(trialSpikes)*(1000/diff(Rspwin));
            bResp                                               = mean(trialBackGround)*(1000/Bkgwin);
            %             eResp = mean(y(Bkgwin+Rspwin(1):Bkgwin+Rspwin(2))); % mean of evoked response
            %             bResp = mean(y(Bkgwin+bg_win(1)+1:Bkgwin+bg_win(2))); % mean of background responses
            
            % subtract background
            tcArray(1)                                          = eResp-bResp;
            %             stdArray(1) = stdResp;
            %             steArray(1) = steResp;
            %             stdResp = std(trialSpikes)*(1000/diff(Rspwin));
            %             steResp = (std(trialSpikes)/sqrt(length(trialSpikes)))*(1000/diff(Rspwin));
            
        end
        
        % step through morphs
        for m = 1:nMorphs
            subdat                                              = select_trials(dat,STEP,morphs(m),FACE,faces(f),TYPE,typeNum,DATE,dates);
            
            imNum                                               = subdat.c(end,STIM); % stimulus image index
            tetPile_imList(f,m+1)                               = imNum;
            
            % add to cumulative sdf matrix
            if ~isempty(subdat.s);
                [r,~]                                           = find(subdat.s>Rspwin(1)&subdat.s<Rspwin(2));
                trialSpikes                                     = histc(r,1:length(subdat.s(:,1)));
                [r,~]                                           = find(subdat.s<0);
                trialBackGround                                 = reshape(histc(r,1:length(subdat.s(:,1))),size(trialSpikes));
                
                eResp                                           = mean(trialSpikes)*(1000/diff(Rspwin));
                bResp                                           = mean(trialBackGround)*(1000/Bkgwin);
                
                % subtract background
                tcArray(m+1)                                    = eResp-bResp;
                %                 stdArray(m+1) = stdResp;
                %                 steArray(m+1) = steResp;
                %                 stdResp = std(trialSpikes)*(1000/diff(Rspwin));
                %                 steResp = (std(trialSpikes)/sqrt(length(trialSpikes)))*(1000/diff(Rspwin));
                
                %                 % check significance of response over background
                %                 [rateT(end+1),p] = ttest(trialBackGround,trialSpikes,'Alpha',0.01);
                
            else
                disp('Make an alert regarding what was just thrown out.')
                keyboard
            end
        end
        
        % save rate limit inclusion criteria
        inclCrit_rateLim_face(f)                                = max(tcArray);
        
        %                 if tcArray(m+1)>=2.5;
        %                     inclCritRateLimLo(i) = inclCritRateLimLo(i)+1;
        %                 end
        %
        %                 if tcArray(m+1)>=5;
        %                     inclCritRateLimHi(i) = inclCritRateLimHi(i)+1;
        %                 end
        
        %         % normalize each trajectory in cumulative date array
        %         morphNormDateArray = [];
        %         for v = 1:length(morphDateArray(1,:))
        %
        %             normVec = cell2mat(morphDateArray(:,v)) ...
        %                +-(min(cell2mat(morphDateArray(:,v))));
        %             %             if max(normVec)~=0
        %             normVec = normVec./max(normVec);
        %             %             end
        %             morphNormDateArray = [morphNormDateArray normVec];
        %         end
        %         faceDateArray{f} = morphNormDateArray';
        
        %         % d-prime for all possible stimulus combinations in tuning curve
        %         pr = nchoosek(1:9,2);
        %         rateD = nan(length(pr(:,1)'),1);
        %         for w = 1:length(pr(:,1)')
        %             dP = (tcArray(pr(w,1))-tcArray(pr(w,2)))/sqrt(stdArray(pr(w,1))+stdArray(pr(w,2)));
        %             rateD(w) = dP^2;
        %         end
        
        %         % rate discriminability criterion
        %         if max(rateD)>=1;
        %             inclCritRateD(i) = inclCritRateD(i)+1;
        %         end
        %
        %         % t-test criterion
        %         [~,pVal] = ttest(backRate,stimRate);
        %         if pVal>.05;
        %             inclCritRateT(i) = inclCritRateT(i)+1;
        %         end
        
        %         % correlation criterion
        %         [~,Pvalue] = corrcoef(x,y);
        %         if Pvalue(2)<0.05;
        %             inclCritRateL(i) = inclCritRateL(i)+1;
        %         end
        
        % aggregate tuning curves for dodec plot
        foo                                                     = f-(floor(f/(length(faces)/2+.1))*12);
        tetPile(foo,:)                                          = tcArray;
        %         tetPile_imList(f) = imNum;
        
        % do some averaging for similar identity values
        for iii=1:length(recomb)
            tetPile_recomb(foo,iii)                             = mean(tcArray(recomb{iii}));
        end
        
        if f==length(faces)/2 || f==length(faces)
            
            %             % normalize accross all responses for a given cell
            %             norm_tetPile = tetPile/max(max(tetPile));
            % norm_recomb_tetPile = tetPile_recomb/max(max(tetPile_recomb));
            
            tempName                                            = dirList(i).name(1:end-4);
            
            %             pop_tetPile                                             = nan(size(tetPile));
            %             for p = 1:length(tetPile(:,1))
            %                 pop_tetPile(p,:) = tetPile(p,:);
            %                 % pop_tetPile(p,:) = tetPile(p,:)/tetPile(p,find(abs(tetPile(p,:))==max(abs(tetPile(p,:)))));
            %             end
            
            % add to population variables
            if f==length(faces)/2
                % humPopResp{i} = pop_tetPile;
                humPopRecomb{i}                                 = tetPile_recomb;
                
                humPopNormResp{i}                               = nan(size(tetPile_recomb));
                for a = 1:length(tetPile_recomb)
                    % subtract minimum resp from tuning curve
                    tcArrayMinusNorm                            = tetPile_recomb(a,:)+-(min(tetPile_recomb(a,:)));
                    % normalize tuning curve
                    humPopNormResp{i}(a,:)                      = tcArrayMinusNorm./max(tcArrayMinusNorm);
                end
                
                foo                                             = tetPile_recomb-tetPile_recomb(1,1); %+ -(min(min(tetPile_recomb)));
                humRankArray                                    = foo./max(max(foo));
                [~,sortInd]                                     = sort(max(humRankArray,[],2),'descend');
                if ~any(any(isnan(humRankArray)))
                    humPopRankResp{i}                               = humRankArray(sortInd,:);
                end
            elseif f==length(faces)
                %                 monPopResp{i} = pop_tetPile;
                monPopRecomb{i}                                 = tetPile_recomb;
                
                monPopNormResp{i}                               = nan(size(tetPile_recomb));
                for a = 1:length(tetPile_recomb)
                    % subtract minimum resp from tuning curve
                    tcArrayMinusNorm                            = tetPile_recomb(a,:)+-(min(tetPile_recomb(a,:)));
                    % normalize tuning curve
                    monPopNormResp{i}(a,:)                      = tcArrayMinusNorm./max(tcArrayMinusNorm);
                end
                
                foo                                             = tetPile_recomb-tetPile_recomb(1,1); %+ -(min(min(tetPile_recomb)));
                monRankArray                                    = foo./max(max(foo));
                [~,sortInd]                                     = sort(max(monRankArray,[],2),'descend');
                if ~any(any(isnan(monRankArray)))
                    monPopRankResp{i}                           = monRankArray(sortInd,:);
                end
            end
            
            %% plot dodecs for each neuron
            if f==length(faces)&&figsOn&&sum(sum([humPopRecomb{i}>0 monPopRecomb{i}>0]))>0
                
                %     %% normalize everything
                %     temp = min([min(humRecomb) min(monRecomb)]);
                %     humRecomb = humRecomb-temp./senorNorm;
                %     monRecomb = monRecomb-temp./senorNorm;
                
                
                % plot spheres
                %                     if sum(sum([humPopRecomb{i}>0 monPopRecomb{i}>0]))>0
                %                     if max(max([humPopRecomb{i} monPopRecomb{i}]))>0
                %                         if ishandle(37)
                %                             clf(37)
                %                         else
                %                             figure(37)
                %                         end
                %                         set(37,'Visible',visTog)
                
                if any(any([humPopRecomb{i}>dodecPlot_cutoff monPopRecomb{i}>dodecPlot_cutoff]))
                    %                             figure
                    %                             set(gcf,'Position',get(gcf,'Position').*[20 1 1 1])
                    data                                        = [humPopRecomb{i};monPopRecomb{i}];
                    temp_lims                                   = [min(min(data)) max(max(data))];
                    savename                                    = fullfile(dodecStr,monkey,[tempName '_dodec.png']);
                    plot_spheres(data,temp_lims,savename,prs,paths,1,visTog); % each arg is 2 cells of 12x6 arrays
                end
                
                %                         if any(any(humPopRecomb{i}>dodecPlot_cutoff))
                %                             plot_spheres_single(humPopRecomb{i},tempName,'h')
                %                             tempPos = get(37,'Position');
                %                             set(37,'PaperUnits','inches','PaperPosition',tempPos/prs);
                %                             savename = fullfile(dodecStr,monkey,[tempName '_dodecHum.png']);
                %                             export_fig(37,savename,['-r' num2str(prs)])
                %                             disp(['saved = ' filename])
                %                             clf(37)
                %                         end
                %
                
                data                                            = [humPopRecomb{i};monPopRecomb{i}];
                full_lims                                       = [min(min(data)) max(max(data))];
                
                if any(any(humPopRecomb{i}>dodecPlot_cutoff))
                    temp_lims                                   = [min(min(humPopRecomb{i})) max(max(humPopRecomb{i}))];
                    savename                                    = fullfile(dodecStr,monkey,[tempName '_dodecHumNoFace.png']);
                    plot_spheres(humPopRecomb{i},temp_lims,savename,prs,paths,0,visTog); % each arg is 2 cells of 12x6 arrays
                end
                
                if any(any(monPopRecomb{i}>dodecPlot_cutoff))
                    temp_lims                                   = [min(min(monPopRecomb{i})) max(max(monPopRecomb{i}))];
                    
                    if regexp(monkey(1),'m')
                        temp_lims                               = full_lims;
                    end
                    
                    savename                                    = fullfile(dodecStr,monkey,[tempName '_dodecMonk.png']);
                    plot_spheres(monPopRecomb{i},temp_lims,savename,prs,paths,1,visTog); % each arg is 2 cells of 12x6 arrays
                    
                    savename                                    = fullfile(dodecStr,monkey,[tempName '_dodecMonkNoFace.png']);
                    plot_spheres(monPopRecomb{i},temp_lims,savename,prs,paths,0,visTog); % each arg is 2 cells of 12x6 arrays
                end
                
                if any(any([humPopRecomb{i}>dodecPlot_cutoff monPopRecomb{i}>dodecPlot_cutoff]))
                    savename                                    = fullfile(dodecStr,monkey,[tempName '_dodecNoFace.png']);
                    plot_spheres(data,full_lims,savename,prs,paths,0,visTog); % each arg is 2 cells of 12x6 arrays
                end
                %                     end
                
            end
            clear tetPile;
        end
    end
    
    inclCrit_rateLim(i)                                         = max(inclCrit_rateLim_face);
    
    %     clear maxSdf minSdf
    % cellDateArray{i} = dateArray;
    % disp([num2str(inclCritRateL(i)) '/' num2str(24) 'trajectories included in inclCritRateL'])
end
% close all

% % handle inclusion criteria
% % inclCritList                                                    = inclCrit_rateLim;
% inclCritRateLimLo                                               = inclCrit_rateLim>low_crit;
% inclCritRateLimHi                                               = inclCrit_rateLim>hi_crit;
% inclCritList                                                    = {inclCritRateLimLo;inclCritRateLimHi;};

% inclCritList = {inclCritRateLimLo;inclCritRateLimHi;}; % inclCritRateD;inclCritRateL
% inclCritSpikeNames = {sNamesList(inclCritRateLimLo~=0,:);sNamesList(inclCritRateLimHi~=0,:);};
% inclCritRateLimLo                                               = inclCritList(inclCritList>low_crit);
% inclCritRateLimHi                                               = inclCritList(inclCritList>hi_crit);


% if figsOn
    pop_recombAve                                                   = cell(length(humPopRecomb),1);
    for i = 1:length(humPopRecomb)
        pop_recombAve{i}                                            = mean(cat(3,humPopRecomb{i},monPopRecomb{i}),3);
    end
    
    % plot population ave. dodec for low rate limit
    foo                                                             = pop_recombAve(inclCrit_rateLim>pop_low_crit);
    ave_dodec                                                       = mean(cat(3,foo{:}),3);
    full_lims                                                       = [min(min(ave_dodec)) max(max(ave_dodec))];
    filename                                                        = ['dodecPopAve_' exp_code '_rateLo' suffix.save '.png'];
    savename                                                        = fullfile(dodecStr,monkey,filename);
    plot_spheres(ave_dodec,full_lims,savename,prs,paths,0,visTog); % each arg is 2 cells of 12x6 arrays
    %     plot_spheres_single_noFace(ave_dodec,savename,prs,paths,suffix,visTog);
    
    % plot population ave. dodec for high rate limit
    % if regexp(monkey(1),'m')
    %     pop_hi_crit                                                 = 6;
    % end
    foo                                                             = pop_recombAve(inclCrit_rateLim>pop_hi_crit);
    ave_dodec                                                       = mean(cat(3,foo{:}),3);
    full_lims                                                       = [min(min(ave_dodec)) max(max(ave_dodec))];
    filename                                                        = ['dodecPopAve_'  exp_code '_rateHi' suffix.save '.png'];
    savename                                                        = fullfile(dodecStr,monkey,filename);
    plot_spheres(ave_dodec,full_lims,savename,prs,paths,0,visTog); % each arg is 2 cells of 12x6 arrays
    %     plot_spheres_single_noFace(ave_dodec,savename,prs,paths,suffix,visTog);
 
    close all
% end

critNames{1}                                                        = 'popRespLimLo';
critNames{2}                                                        = 'popRespLimHi';

savename                                                        = fullfile(monkeyDataDir,[exp_code 'Data' suffix.save '.mat']);
save(savename,'inclCrit_rateLim','recombMorphs','critNames',...
    'humPopRecomb','monPopRecomb','humPopRankResp','monPopRankResp');
disp({'Saved monkey data:'; savename})

% keyboard
%
% % nonContigCells = zeros(1,length(popTuningCurves));
% % for c = 1:length(popTuningCurves)
% %     if sum(cellfun(@isempty,popTuningCurves{c}))>=1
% %         nonContigCells(c) = 1;
% %     end
% % end
% % nonContigCellInds = find(nonContigCells);
%
% % save inclusion criteria
% savename = [dataStr  'inclCrit' monkey '.mat'];
% save(savename,'inclCritList','inclCritSpikeNames')
%
% % save variables
% savename = [dataStr  'identTraj' monkey '.mat'];
% save(savename)

end