function insp_dat_identTrajPop(paths,monkeys,suffix)
% loop through population data files for each monkey and make population plots
%
% Usage: insp_dat_identTrajPop(paths,monkey,suffix)
% MONKEY is the monkey's name
% PATHS is the project path structure containing at minimum paths.rare field
% SUFFIX is a string appended to the end of each filename
%
% NOTE: stimset includes 24 humans and 12 monkeys - 12 hum & 6 mon face pairs (3 trajectories[2 rad;1 tan])
%
% last modified 6-24-16
% apj


% % define stimuli
% [stimNames,~]                                   = getStimList(proj,fullfile(paths.physio,'ptb'));

% % load appropriate image files for display
% img                                                             = load(fullfile(paths.stim,'imgFaceLearningSet1OnWhite.mat'));
% stim_images                                                     = img;

% load(fullfile(paths.stim,'imgRadTanSetWhite.mat'));

% setup
figsOn                                                          = 0; % switch for figures (on vs. off)
visTog                                                          = 'off';
Bkgwin                                                          = 300; % length of prestimulus period (msec)
Plotwin                                                         = [-300 600]; % plot display size (msec)
Rspwin                                                          = [100 400]; % time window for evoked rate
prs                                                             = 150; % print resolution
fullscr                                                         = get(0,'ScreenSize'); % get screen size
ax_Fsz                                                          = 12; % axes fontsize
tTl_Fsz                                                         = 11; % title fontsize
sPtTl_Fsz                                                       = 13.5; % super-title fontsize
ax_Lw                                                           = 2; % axes linewidth
erBr_Lw                                                         = 3.5; % error-bar linewidth
low_crit                                                        = 2.5; % minimum rate (low) for pop. analysis
hi_crit                                                         = 5; % minimum rate (high) for pop. analysis
crit_lists                                                      = {low_crit; hi_crit};
% figJit                                          = round(rand*10);           % avoid using the same fig number in parallel
% title_Fsz                                       = 10;                       % title fontsize
% stim_Lw                                         = 4.5;                      % Stimulus presentation linewidth
% sdf_Lw                                          = 2;                        % lineWidth for sdf's
% popTC_Lw                                        = 4;                        % lineWidth for sdf's in population plot
% figPos1                                         = [fullscr(3)/8   fullscr(4)/10   fullscr(3)/2   fullscr(4)/4]; % [x y width height]
% if regexp(visTog,'off')
%     figPos1(1)                                  = figPos1(1)+5000;
% end

% % face numbers
% humFacePairs                                    = 1:12;
% monFacePairs                                    = 1:6;
% % faceList                                        = [humFacePairs monFacePairs];
% humNorm                                         = 535;
% monNorm                                         = 536;

popResults_dir                                                  = fullfile(paths.results,'population');


% % setup
% % % Plotwin         = [-300 600]; % plot display size (msec)
% % % figPos          = [340   120   1300   400]; % [x y width height]
% %
% % bkground_len    = 300;          % length of prestimulus period (msec)
% % resp_win        = [100 400];    % time window for evoked rate
% % recomb          = {1;2:4;5:6;7;8;9}; % values for averaging
%
% % % face numbers
% % human_facenums  = 1:12;
% % monkey_facenums = 101:112;
% % faces           = [human_facenum monkey_facenums];
% % object = 201:212;

% define subplot positions
positions                                                       = repmat(1:length(monkeys),length(monkeys),1)+...
    repmat((length(monkeys)+1).*[0:length(monkeys)-1]',1,length(monkeys));

% % load all image files for disp (stimscreen 1 & 2)
% img = load(fullfile(paths.stim,'imgFaceLearningSet1OnWhite.mat'));
% % img = load([paths.stim 'imgFaceLearningSet1.mat']);

% step through inclusion criteria lists
for g = 1:length(crit_lists)
    
    %% COLLECT RESPONSES
    
    cellCount                                                   = cell(length(monkeys),1);
    aveMonNormResps                                             = cell(length(monkeys),1);
    aveHumNormResps                                             = cell(length(monkeys),1);
    aveSpeciesNormResps                                         = cell(length(monkeys),1);
    errSpeciesNormResps                                         = cell(length(monkeys),1);
    % setup color table for each monkey in study
    monkey_colors                                               = parula(length(monkeys));
    for mo = 1:length(monkeys)
        
        monkey                                                  = monkeys{mo};
        %         monkeyCode      = decodeMonkeyName(monkey);
        exp_code                                                = [monkey(1) paths.code];
        
        monkeyDataDir                                           = fullfile(paths.results,'fine',monkey); % directory for saved data
        dataname                                                = fullfile(monkeyDataDir,[exp_code  'Data'  suffix.load '.mat']);
        monkey_data                                             = load(dataname);
        
        %         saveDir = fullfile(paths.results,'population',monkey);
        % load('/projects/apjones/results/faceLearning/heatSpk/identTrajtoroid.mat');
        
        %         humPopTemp = monkey_data.humPopResp;
        %         monPopTemp = monkey_data.monPopResp;
        
        %% from identTraj
        %          if tcArray(m+1)>=2.5;
        % %                     inclCritRateLimLo(i) = inclCritRateLimLo(i)+1;
        % %                 end
        % %
        % %                 if tcArray(m+1)>=5;
        % %                     inclCritRateLimHi(i) = inclCritRateLimHi(i)+1;
        % %                 end
        %%
        
        % collect inclusion criteria
        if g==1
            inclCrit                                            = monkey_data.inclCrit_rateLim>low_crit;
        else
            
            if mo==3
                hi_crit = 6
            else
                hi_crit = 5
            end
            
            inclCrit                                            = monkey_data.inclCrit_rateLim>hi_crit;
        end
        %         inclCritList                                            = {inclCritRateLimLo;inclCritRateLimHi;};
        
        indLocs                                                 = find(inclCrit);
        inclCritName                                            = monkey_data.critNames{g};
        
        if mo==3&&g==2
            inclCrit([3 16 19 20 24 32 36 41]) = 0;
        end
        
        %         inclCritSpikeNames = monkey_data.inclCritNames{g};
        %         if mo==3&g==1
        % %             keyboard
        %             inclCrit([10 11 12 17 18 19 40 41 42 ...
        %                 2 3 4 5 13 16 25 26 27 43]) = 0;
        % %                     % 14 15 16 17 18 19 20 21 48
        %             inclCritSpikeNames([10 11 12 17 18 19 40 41 42 ...
        %                 2 3 4 5 13 16 25 26 27 43],:) = [];
        %         end
        
        % aggregate neuron-counts for each monkey
        cellCount{mo}                                           = num2str(sum(inclCrit));
        
        %     cutout = find(inclCrit==0);
        %     humPopTemp(cutout) = [];
        %     monPopTemp(cutout) = [];
        
        %     % remove inhibitory responses
        %     co = zeros(length(humPopTemp),1);
        %     for r = 1:length(humPopTemp)
        %         for rr = 1:length(humPopTemp{r}(:,1))
        %             if humPopTemp{r}(rr,:)<0ad
        %                 co(r)= 1;
        %             end
        %         end
        %     end
        %     cutout = find(co==1);
        %     humPopTemp(cutout) = [];
        %     monPopTemp(cutout) = [];
        
        %         inclIndices = (inclCrit>0);
        %     inclIndices(cutout) = 0;
        
        
        %% SPECIES AVERAGES
        
        monNormResps                                            = [];
        humNormResps                                            = [];
        for r = 1:length(indLocs)
            
            % normalize responses across species
            monkeyFace_resp                                     = monkey_data.monPopRecomb{indLocs(r)};
            humanFace_resp                                      = monkey_data.humPopRecomb{indLocs(r)};
            %             monResps = [monResps; monkeyFace_resp];
            %             humResps = [humResps; humanFace_resp];
            
            % normalize to minimum/maximum responses
            foo                                                 = min(min([humanFace_resp; monkeyFace_resp]));
            maxNorm                                             = max(max([humanFace_resp; monkeyFace_resp]))-foo;
            
            % aggregate
            monNormResps                                        = [monNormResps; (monkeyFace_resp-foo)./maxNorm];
            humNormResps                                        = [humNormResps; (humanFace_resp-foo)./maxNorm];
        end
        
        %         if g==2
        %             savename = [dataStr  'popSpeciesrhombus.mat'];
        %             rhombPopSpecies = load(savename);
        
        % aggregate normalized responses
        aveMonNormResps{mo}                                     = monNormResps;
        aveHumNormResps{mo}                                     = humNormResps;
        
        aveSpeciesNormResps{mo}                                 = [mean(monNormResps(:,6)) mean(humNormResps(:,6))];
        errSpeciesNormResps{mo}                                 = [std(monNormResps(:,6))/sqrt(length(monNormResps(:,6))) ...
            std(humNormResps(:,6))/sqrt(length(humNormResps(:,6)))];
        
        %         % statistical test
        %         [H,m1Pnorm] = ttest2(monNormResps(:,6),humNormResps(:,6));
        %         [H,m2Pnorm] = ttest2(rhombPopSpecies.monNormResps(:,6),rhombPopSpecies.humNormResps(:,6));
        
        % aggregate raw responses
        humResps                                                = cell2mat(monkey_data.humPopRecomb(indLocs));
        monResps                                                = cell2mat(monkey_data.monPopRecomb(indLocs));
        
        aveMonResps3d{mo}                                       = cat(3,monkey_data.monPopRecomb{indLocs}); % cell2mat(monkey_data.monPopRankResp(indLocs));
        aveHumResps3d{mo}                                       = cat(3,monkey_data.humPopRecomb{indLocs});
        
        %         m2Hums = cell2mat(rhombPopSpecies.humRawResps);
        %         m2Mons = cell2mat(rhombPopSpecies.monRawResps);
        aveMonResps{mo}                                         = monResps;
        aveHumResps{mo}                                         = humResps;
        
        aveSpeciesResps{mo}                                     = [mean(monResps(:,6)) mean(humResps(:,6))];
        errSpeciesResps{mo}                                     = [std(monResps(:,6))/sqrt(length(monResps(:,6))) ...
            std(humResps(:,6))/sqrt(length(monResps(:,6)))];
        
        %         % statistical test
        %         [H,m1Praw] = ttest2(m1Mons(:,6),m1Hums(:,6));
        %         [H,m2Praw] = ttest2(m2Mons(:,6),m2Hums(:,6));
        
        %         figure; hold on
        %         for w = 1:length(indLocs)
        %             foo = ceil(sqrt(length(indLocs)));
        %             subplot(foo,foo,w)
        %             plot(monkey_data.monPopRankResp{indLocs(w)}')
        %             title(num2str(w))
        %         end
        %         keyboard
        aveMonRankResps{mo}                                     = cell2mat(monkey_data.monPopRankResp(indLocs));
        aveHumRankResps{mo}                                     = cell2mat(monkey_data.humPopRankResp(indLocs));
        
        foo                                                     = monkey_data.monPopRecomb(indLocs);
        dim                                                     = ndims(foo{1}); % Get the number of dimensions for your arrays
        aveMonIdentResps{mo}                                    = mean(cat(dim+1,foo{:}),dim+1); % Get the mean across arrays
        foo                                                     = monkey_data.humPopRecomb(indLocs);
        aveHumIdentResps{mo}                                    = mean(cat(dim+1,foo{:}),dim+1); % Get the mean across arrays
        
        aveMonRankResps3d{mo}                                   = cat(3,monkey_data.monPopRankResp{indLocs}); % cell2mat(monkey_data.monPopRankResp(indLocs));
        aveHumRankResps3d{mo}                                   = cat(3,monkey_data.humPopRankResp{indLocs});
    end
    
    
    %% SPECIES PLOTS
    
    % loop through monkeys and setup legend text
    legend_text1                                                 = cell(length(monkeys),1);
    for mo = 1:length(monkeys)
        string1                                                 = ['M' num2str(mo)]; % monkeys{mo}
        string2                                                 = ['(n = ' cellCount{mo} ')'];
        %         blank_piece                                             = blanks(diff([length(string1) length(string2)]));
        legend_text1{mo}                                        = [string1 ' ' string2];
    end
    legend_text2                                                = {'monkey faces';'human faces'};
    
    %% plot species-averaged,normalized responses
    model_series                                                = cell2mat(aveSpeciesNormResps)';
    model_error                                                 = cell2mat(errSpeciesNormResps)';
    
    figure;
    bH                                                          = bar(model_series); hold on
    for mc = 1:length(monkey_colors)
        set(bH(mc),'FaceColor',monkey_colors(mc,:))
    end
    
    set(bH,'BarWidth',1);
    set(gca,'YGrid','on','GridLineStyle','-','XTick',[1 2],...
        'XTickLabel',legend_text2,'LineWidth',ax_Lw,...
        'Fontsize',ax_Fsz,'Fontweight','bold')
    lh                                                          = legend(legend_text1,'Box','off');
    set(lh,'Location','BestOutside','Orientation','horizontal')
    numgroups                                                   = size(model_series,1);
    numbars                                                     = size(model_series,2);
    groupwidth                                                  = min(0.8,numbars/(numbars+1.5));
    % align error bar with individual bar
    for i = 1:numbars
        foo                                                     = (1:numgroups)-groupwidth/2+(2*i-1)*groupwidth/(2*numbars);
        errorbar(foo,model_series(:,i),model_error(:,i),'k','Linestyle','none','LineWidth',ax_Lw);
    end
    line(xlim,[min(ylim) min(ylim)],'Color','k','LineWidth',ax_Lw);
    line([min(xlim) min(xlim)],ylim,'Color','k','LineWidth',ax_Lw);
    set(gca,'Position',get(gca,'Position').*[1 .85 1 1])
    title(['incl. crit.: ' inclCritName],'Fontsize',ax_Fsz,'Fontweight','bold');
    ylabel('ave. norm. response');
    
    tempPos                                                     = get(gcf,'Position'); % [550   530   550   330];
    set(gcf,'Color','w','Position',tempPos);
    
    savename                                                    = fullfile(popResults_dir,['popNormAveSpecies_crit' num2str(g) suffix.save]);
    saveas(gcf,[savename '.eps'],'epsc');
    export_fig(gcf,[savename '.png'],'-nocrop',['-r' num2str(prs)])
    disp(['saved: ' savename])
    close(gcf)
    
    % normalized, species-ave responses, grouped by subject
    model_series                                                = cell2mat(aveSpeciesNormResps);
    model_error                                                 = cell2mat(errSpeciesNormResps);
    
    figure;
    bH                                                          = bar(model_series); hold on
    %     for mc = 1:length(monkey_colors)
    %         set(bH(mc),'FaceColor',monkey_colors(mc,:))
    %     end
    set(bH,'BarWidth',1);
    set(gca,'YGrid','on','GridLineStyle','-','XTick',[1 2 3],...
        'XTickLabel',legend_text1,'LineWidth',ax_Lw,...
        'Fontsize',ax_Fsz,'Fontweight','bold')
    
    lh                                                          = legend(legend_text2,'Box','off');
    set(lh,'Location','BestOutside','Orientation','horizontal')
    numgroups                                                   = size(model_series,1);
    numbars                                                     = size(model_series,2);
    groupwidth                                                  = min(0.8,numbars/(numbars+1.5));
    % align error bar with individual bar
    for i = 1:numbars
        foo                                                     = (1:numgroups)-groupwidth/2+(2*i-1)*groupwidth/(2*numbars);
        errorbar(foo,model_series(:,i),model_error(:,i),'k','Linestyle','none','LineWidth',ax_Lw);
    end
    line(xlim,[min(ylim) min(ylim)],'Color','k','LineWidth',ax_Lw);
    line([min(xlim) min(xlim)],ylim,'Color','k','LineWidth',ax_Lw);
    set(gca,'Position',get(gca,'Position').*[1 .85 1 1])
    title(['incl. crit.: ' inclCritName],'Fontsize',ax_Fsz,'Fontweight','bold');
    ylabel('ave. norm. response');
    set(gcf,'Color','w','Position',tempPos);
    
    savename                                                    = fullfile(popResults_dir,['popNormAveSpeciesXSubj_crit' num2str(g) suffix.save]);
    saveas(gcf,[savename '.eps'],'epsc');
    export_fig(gcf,[savename '.png'],'-nocrop',['-r' num2str(prs)])
    disp(['saved: ' savename])
    close(gcf)
    
    
    %% plot species-averaged responses
    model_series                                                = cell2mat(aveSpeciesResps')';
    model_error                                                 = cell2mat(errSpeciesResps')';
    
    figure; hold on
    bH                                                          = bar(model_series); hold on
    for mc = 1:length(monkey_colors)
        set(bH(mc),'FaceColor',monkey_colors(mc,:))
    end
    set(bH,'BarWidth',1);    % The bars will now touch each other
    set(gca,'YGrid','on','GridLineStyle','-','XTick',[1 2],...
        'XTickLabel',legend_text2,'LineWidth',ax_Lw,...
        'Fontsize',ax_Fsz,'Fontweight','bold')
    
    lh                                                          = legend(legend_text1,'Box','off');
    set(lh,'Location','BestOutside','Orientation','horizontal')
    numgroups                                                   = size(model_series,1);
    numbars                                                     = size(model_series,2);
    groupwidth                                                  = min(0.8,numbars/(numbars+1.5));
    % align error bar with individual bar
    for i = 1:numbars
        foo                                                     = (1:numgroups)-groupwidth/2+(2*i-1)*groupwidth/(2*numbars);
        errorbar(foo,model_series(:,i),model_error(:,i),'k','Linestyle','none','LineWidth',ax_Lw);
    end
    line(xlim,[min(ylim) min(ylim)],'Color','k','LineWidth',ax_Lw);
    line([min(xlim) min(xlim)],ylim,'Color','k','LineWidth',ax_Lw);
    set(gca,'Position',get(gca,'Position').*[1 .85 1 1])
    title(['incl. crit.: ' inclCritName],'Fontsize',ax_Fsz,'Fontweight','bold');
    ylabel('ave. response (Hz)')
    set(gcf,'Color','w','Position',tempPos);
    
    savename                                                    = fullfile(popResults_dir,['popAveSpecies_crit' num2str(g) suffix.save]);
    saveas(gcf,[savename '.eps'],'epsc');
    export_fig(gcf,[savename '.png'],'-nocrop',['-r' num2str(prs)])
    disp(['saved: ' savename])
    close(gcf)
    
    % species-ave responses, grouped by subject
    model_series                                                = cell2mat(aveSpeciesResps');
    model_error                                                 = cell2mat(errSpeciesResps');
    
    figure;
    bH                                                          = bar(model_series); hold on
    %     for mc = 1:length(monkey_colors)
    %         set(bH(mc),'FaceColor',monkey_colors(mc,:))
    %     end
    set(bH,'BarWidth',1);
    set(gca,'YGrid','on','GridLineStyle','-','XTick',[1 2 3],...
        'XTickLabel',legend_text1,'LineWidth',ax_Lw,...
        'Fontsize',ax_Fsz,'Fontweight','bold')
    
    lh                                                          = legend(legend_text2,'Box','off');
    set(lh,'Location','BestOutside','Orientation','horizontal')
    numgroups                                                   = size(model_series,1);
    numbars                                                     = size(model_series,2);
    groupwidth                                                  = min(0.8,numbars/(numbars+1.5));
    % align error bar with individual bar
    for i = 1:numbars
        foo                                                     = (1:numgroups)-groupwidth/2+(2*i-1)*groupwidth/(2*numbars);
        errorbar(foo,model_series(:,i),model_error(:,i),'k','Linestyle','none','LineWidth',ax_Lw);
    end
    line(xlim,[min(ylim) min(ylim)],'Color','k','LineWidth',ax_Lw);
    line([min(xlim) min(xlim)],ylim,'Color','k','LineWidth',ax_Lw);
    set(gca,'Position',get(gca,'Position').*[1 .85 1 1])
    title(['incl. crit.: ' inclCritName],'Fontsize',ax_Fsz,'Fontweight','bold');
    ylabel('ave. response (Hz)')
    set(gcf,'Color','w','Position',tempPos);
    
    savename                                                    = fullfile(popResults_dir,['popAveSpeciesXSubj_crit' num2str(g) suffix.save]);
    saveas(gcf,[savename '.eps'],'epsc');
    export_fig(gcf,[savename '.png'],'-nocrop',['-r' num2str(prs)])
    disp(['saved: ' savename])
    
    close all
    
    %% plot morphLevel-averaged responses for each monkey
    hLine                                                       = [];
    figure; hold on
    for mo = 1:length(monkeys)
        model_series_mon                                        = aveMonResps{mo};
        model_series_hum                                        = aveHumResps{mo};
        
        %     nmbReps = length(humNormResps(:,1)); % number of replications
        %     [~,~,stats] = anova2([humNormResps; monNormResps],nmbReps);
        %     c = multcompare(stats);
        
        tempErr                                                 = std(model_series_mon)./sqrt(length(model_series_mon(:,1)));
        hL                                                      = errorbar(mean(model_series_mon),tempErr, ...
            'Color',monkey_colors(mo,:),'LineWidth',erBr_Lw);
        hLine                                                   = [hLine; hL];
        tempErr                                                 = std(model_series_hum)./sqrt(length(model_series_hum(:,1)));
        errorbar(mean(model_series_hum),tempErr,'Color',monkey_colors(mo,:),'LineWidth',erBr_Lw,'LineStyle','--')
    end
    lh                                                          = legend(hLine,legend_text1,'Box','off');
    set(lh,'Location','BestOutside','Orientation','horizontal');
    
    yLims                                                       = ylim;
    %         text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) 'cells)'],...
    %             'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
    set(gca,'XTick',1:6,'XTickLabel',monkey_data.recombMorphs,'LineWidth',ax_Lw,...
        'YTick',yLims(1):diff(yLims)/5:yLims(2),'Fontsize',ax_Fsz,'Fontweight','bold')
    xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
    ylabel('pop. ave. response','Fontsize',ax_Fsz,'Fontweight','bold')
    set(gca,'YLim',yLims,'Position',get(gca,'Position').*[1 .85 1 1])
    title('pop. ave. morph-level responses','Fontsize',tTl_Fsz,'Fontweight','bold')
    tempPos                                                     = get(gcf,'Position').*[1 1 1 1.5]; %[550   530   550   330];
    set(gcf,'Color','w','Position',tempPos);
    %    ,'PaperUnits','inches','PaperPosition',tempPos/prs
    
    savename                                                    = fullfile(popResults_dir,['popAveMorphLevel_crit' num2str(g) suffix.save]);
    saveas(gcf,[savename '.eps'],'epsc');
    export_fig(gcf,[savename '.png'],'-nocrop',['-r' num2str(prs)]);
    disp(['saved: ' savename])
    close all
    
    %% plot morphLevel-averaged normalized responses for each monkey
    hLine                                                       = [];
    figure; hold on
    for mo = 1:length(monkeys)
        model_series_mon                                        = aveMonNormResps{mo};
        model_series_hum                                        = aveHumNormResps{mo};
        
        %         if mo==3
        %             keyboard
        %         end
        
        % get indices
        %         [~,sortI] = sort(model_series_hum(:,2),'descend');
        %         % get top 50 response amplitudes
        %         top50 = model_series_hum(sortI(end-500:end),2);
        %         figure; plot(mean(model_series_hum(sortI(end-300:end),:)))
        
        %     nmbReps = length(humNormResps(:,1)); % number of replications
        %     [~,~,stats] = anova2([humNormResps; monNormResps],nmbReps);
        %     c = multcompare(stats);
        
        tempErr                                                 = std(model_series_mon)./sqrt(length(model_series_mon(:,1)));
        hL                                                      = errorbar(mean(model_series_mon),tempErr, ...
            'Color',monkey_colors(mo,:),'LineWidth',erBr_Lw);
        hLine                                                   = [hLine; hL];
        tempErr                                                 = std(model_series_hum)./sqrt(length(model_series_hum(:,1)));
        errorbar(mean(model_series_hum),tempErr,'Color',monkey_colors(mo,:),'LineWidth',erBr_Lw,'LineStyle','--')
    end
    %          line(xLims,[0 0],'Color','k');
    %          line([1 1],yLims,'Color','k');
    lh                                                          = legend(hLine,legend_text1,'Box','off');
    set(lh,'Location','BestOutside','Orientation','horizontal');
    
    yLims                                                       = ylim;
    %         text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) 'cells)'],...
    %             'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
    set(gca,'XTick',1:6,'XTickLabel',monkey_data.recombMorphs,'LineWidth',ax_Lw,...
        'YTick',[yLims(1):diff(yLims)/5:yLims(2)],'Fontsize',ax_Fsz,'Fontweight','bold')
    xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
    ylabel('pop. ave. normalized resp.','Fontsize',ax_Fsz,'Fontweight','bold')
    set(gca,'YLim',yLims,'Position',get(gca,'Position').*[1 .85 1 1])
    title('pop. ave. norm. morph-level resp.','Fontsize',tTl_Fsz,'Fontweight','bold')
    set(gcf,'Color','w','Position',tempPos);
    
    savename                                                    = fullfile(popResults_dir,['popNormAveMorphLevel_crit' num2str(g) suffix.save]);
    saveas(gcf,[savename '.eps'],'epsc');
    export_fig(gcf,[savename '.png'],'-nocrop',['-r' num2str(prs)]);
    disp(['saved: ' savename])
    close all
    
    %     if g==2
    %         % both animals
    %         figure; hold on
    %         % toroid
    %         tempMean = mean(monNormResps);
    %         tempErr = std(monNormResps)./sqrt(length(monNormResps(:,1)));
    %         errorbar(tempMean,tempErr,'g','LineWidth',erBr_Lw)
    %         text(length(tempMean)+.1,tempMean(end),'M1',...
    %             'Fontsize',ax_Fsz,'Fontweight','bold','LineWidth',ax_Lw)
    %
    %         tempMean = mean(humNormResps);
    %         tempErr = std(humNormResps)./sqrt(length(humNormResps(:,1)));
    %         errorbar(tempMean,tempErr,'LineWidth',erBr_Lw)
    %         text(length(tempMean)+.1,tempMean(end),'M1',...
    %             'Fontsize',ax_Fsz,'Fontweight','bold','LineWidth',ax_Lw)
    %
    %         % rhombus
    %         tempMean = mean(rhombPopSpecies.monNormResps);
    %         tempErr = std(rhombPopSpecies.monNormResps)./sqrt(length(rhombPopSpecies.monNormResps(:,1)));
    %         errorbar(tempMean,tempErr,'g','LineWidth',erBr_Lw)
    %         text(length(tempMean)+.1,tempMean(end),'M2',...
    %             'Fontsize',ax_Fsz,'Fontweight','bold','LineWidth',ax_Lw)
    %
    %         tempMean = mean(rhombPopSpecies.humNormResps);
    %         tempErr = std(rhombPopSpecies.humNormResps)./sqrt(length(rhombPopSpecies.humNormResps(:,1)));
    %         errorbar(tempMean,tempErr,'LineWidth',erBr_Lw)
    %         text(length(tempMean)+.1,tempMean(end),'M2',...
    %             'Fontsize',ax_Fsz,'Fontweight','bold','LineWidth',ax_Lw)
    %
    %         yLims = ylim;
    %         text(2.5,yLims(1)+diff(yLims)*.75,{['M1 (n = ' num2str(sum(inclIndices)) 'cells)'];...
    %             ['M2 (n = ' num2str(length(rhombPopSpecies.monNormResps(:,1))/12) 'cells)']},...
    %             'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw,...
    %             'HorizontalAlignment','center')
    %         set(gca,'XTick',1:6,'XTickLabel',recombMorphs,'LineWidth',ax_Lw,...
    %             'YTick',[yLims(1):diff(yLims)/5:yLims(2)],'Fontsize',ax_Fsz,'Fontweight','bold')
    %         xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
    %         ylabel('normalized response','Fontsize',ax_Fsz,'Fontweight','bold')
    %         title('population ave species responses','Fontsize',tTl_Fsz,'Fontweight','bold')
    %         lH = legend('monkey faces','human faces');
    %         yLims = ylim;
    %         tempPosLeg = get(lH,'Position');
    %         set(lH,'Position',tempPosLeg.*[.45 .95 1 1],'Box','off','LineWidth',ax_Lw);
    %         line([1 1],yLims,'Color','k');
    %         set(gca,'YLim',yLims)
    %
    %         %     tempPos = [550   530   550   330]; % get(gcf,'Position');
    %         set(gcf,'Color','w')
    %         %     ,'Position',tempPos,'PaperUnits','inches',...
    %         %         'PaperPosition',tempPos/prs);
    %
    %         savename = [popResults_dir 'bothPopAveSpecies' num2str(g) suffix '.eps'];
    %         saveas(gcf,savename,'epsc');
    %         savename = [popResults_dir 'bothPopAveSpecies' num2str(g) suffix '.png'];
    %         export_fig(gcf,savename,'-nocrop',['-r' num2str(prs)])
    %         close(gcf)
    %     end
    
    %     % human face population plot
    %     dim = ndims(humPopNormResp{1});          % get number of dimensions for array
    %     M = cat(dim+1,humPopNormResp{indLocs});        % convert to (dim+1)-dimensional matrix
    %     meanArray = mean(M,dim+1);  % get mean across arrays
    %     steArray = std(M,[],dim+1)./sqrt(length(M));  % get mean across arrays
    %
    %     figure
    %     subplot(211); hold on
    %     errorbar(meanArray',steArray','LineWidth',erBr_Lw);
    %     yLims = ylim;
    %     xLims = xlim;
    %     line([1 1],yLims,'Color','k');
    %     text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) ' cells)'],...
    %         'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
    %     set(gca,'XTick',1:6,'XTickLabel',recombMorphs,'LineWidth',ax_Lw,...
    %         'YTick',[yLims(1):diff(yLims)/5:yLims(2)],'Fontsize',ax_Fsz,'Fontweight','bold')
    %     ylabel('normalized response','Fontsize',ax_Fsz,'Fontweight','bold')
    %     title('human face population responses','Fontsize',tTl_Fsz,'Fontweight','bold')
    %
    
    
    %% refine response data
    
    for mo = 1:length(monkeys)
        
        model_series_mon                                        = aveMonRankResps{mo};
        model_series_hum                                        = aveHumRankResps{mo};
        model_series_mon3                                       = aveMonRankResps3d{mo};
        model_series_hum3                                       = aveHumRankResps3d{mo};
        
        model_series_ave_mon3                                   = aveMonResps3d{mo};
        model_series_ave_hum3                                   = aveHumResps3d{mo};
        
        meanIdentArray_mon                                      = mean(model_series_mon3, 3);
        meanIdentArray_hum                                      = mean(model_series_hum3, 3);
        meanIdentArray_err_mon                                  = std(model_series_mon3, 0, 3)./length(model_series_mon3(1,1,:));
        meanIdentArray_err_hum                                  = std(model_series_hum3, 0, 3)./length(model_series_hum3(1,1,:));
        
        meanIdentArrayAve_mon                                   = mean(model_series_ave_mon3, 3);
        meanIdentArrayAve_hum                                   = mean(model_series_ave_hum3, 3);
        meanIdentArrayAve_err_mon                               = std(model_series_ave_mon3, 0, 3)./length(model_series_ave_mon3(1,1,:));
        meanIdentArrayAve_err_hum                               = std(model_series_ave_hum3, 0, 3)./length(model_series_ave_hum3(1,1,:));
        
        % get resp ste across arrays
        steArray                                                = std(model_series_mon)./sqrt(length(model_series_mon));
        
        % get slopes of tuning curves
        slopes                                                  = [];
        y0                                                      = monkey_data.recombMorphs;
        for e = 1:length(model_series_mon)
            P                                                   = polyfit(y0,model_series_mon(e,:),1); %[slope intercept]
            slopes                                              = [slopes; P(1)];
        end
        
        %% plot cell-by-cell tuning curves
        
        % ranked averages- monkey
        figure(2323); hold on
        set(2323,'Position',[500 200 550 900]);
        subplot(2,1,1)
        %         yLims = [-1 1];
        errorbar(meanIdentArray_mon',meanIdentArray_err_mon','LineWidth',erBr_Lw);
        line(xlim,[0 0],'Color','k');
        line([1 1],[-1 1],'Color','k');        
        %         ylim(yLims);
        % text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) ' cells)'],...
        %   'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
        yTicks                                                  = min([min(ylim) -.05]):diff([0 max(ylim)])/4:max(ylim);
        set(gca,'XLim',[0 7],'YLim',[-1 1],'XTick',1:6,'XTickLabel',monkey_data.recombMorphs,...
            'LineWidth',ax_Lw,'Fontsize',ax_Fsz,'Fontweight','bold',...
            'YTick',[yTicks(yTicks<0) 0 yTicks(yTicks>0)]); % diff(yLims)/5:yLims(2)
        
        xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
        ylabel('normalized response','Fontsize',ax_Fsz,'Fontweight','bold')
        title('rank ordered resp. - monkey faces','Fontsize',tTl_Fsz,'Fontweight','bold')

        % unranked averages- monkey
        figure(4242); hold on
        set(4242,'Position',[500 200 550 900]);
        subplot(2,1,1)
        %         yLims = [-1 1];
        errorbar(meanIdentArrayAve_mon',meanIdentArrayAve_err_mon','LineWidth',erBr_Lw);
        line([1 1],ylim,'Color','k');
        %         ylim(yLims);
        % text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) ' cells)'],...
        %   'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
%         yTicks                                                  = min(ylim -.05):diff([min(ylim -.05) max(ylim)])/4:max(ylim);
        set(gca,'XLim',[0 7],'YLim',ylim,'XTick',1:6,'XTickLabel',monkey_data.recombMorphs,...
            'LineWidth',ax_Lw,'Fontsize',ax_Fsz,'Fontweight','bold');
        xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
        ylabel('ave. response(Hz)','Fontsize',ax_Fsz,'Fontweight','bold')
        title('ave. resp. - monkey faces','Fontsize',tTl_Fsz,'Fontweight','bold')

        % plot population tuning curve
        figure(2525); hold on
        subplot(length(monkeys),length(monkeys)+1,positions)
        %         yLims = [-1 1];
        errorbar(mean(model_series_mon),steArray,'LineWidth',erBr_Lw,'Color',monkey_colors(mo,:)); hold on
        %         ylim(yLims);
        % text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) ' cells)'],...
        %   'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
        yTicks                                                  = min([min(ylim) -.05]):diff([0 max(ylim)])/4:max(ylim);
        set(gca,'XLim',[0 7],'XTick',1:6,'XTickLabel',monkey_data.recombMorphs,...
            'LineWidth',ax_Lw,'Fontsize',ax_Fsz,'Fontweight','bold',...
            'YTick',[yTicks(yTicks<0) 0 yTicks(yTicks>0)]); % diff(yLims)/5:yLims(2)
        
        xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
        ylabel('normalized response','Fontsize',ax_Fsz,'Fontweight','bold')
        title('rank ordered resp. - monkey faces','Fontsize',tTl_Fsz,'Fontweight','bold')
        
        % plot histogram of tuning curve slopes
        subplot(length(monkeys),length(monkeys)+1,mo*(length(monkeys)+1))
        % axes('Position',get(gca,'Position').*[2.25 1.25 .25 .2]); hold on; box on;
        hist(slopes,30);
        hP                                                      = findobj(gca,'Type','patch');
        set(hP,'FaceColor',monkey_colors(mo,:));
        axis tight;
        line([0 0],ylim,'Color','k');
        title([monkeys{mo} '(n = ' cellCount{mo} ')'],'Fontsize',8,'Fontweight','bold');
        xlim([-0.02 0.02]);
        %  maxTick = str2num(num2str(floor(max(xlim)*1000)/1000));
        set(gca,'YTick',max(ylim),'Fontsize',8,'Fontweight','bold','LineWidth',ax_Lw);
        %  ,'XTick',[min(xlim) 0 max(xlim)],'XTickLabel',num2str([min(xlim) 0 max(xlim)]),...
        if mo==length(monkeys)
            set(gca,'XTick',[-0.02 0 0.02]);
        else
            set(gca,'XTick',[]);
        end
        
        % get resp ste across arrays
        steArray                                                = std(model_series_hum)./sqrt(length(model_series_hum));
        
        % get slopes of tuning curves
        slopes                                                  = [];
        y0                                                      = monkey_data.recombMorphs;
        for e = 1:length(model_series_hum)
            P                                                   = polyfit(y0,model_series_hum(e,:),1); %[slope intercept]
            slopes                                              = [slopes; P(1)];
        end
        
        %% plot cell-by-cell tuning curves
        
        % ranked averages- human
        figure(2323); hold on
        set(2323,'Color','w');
        subplot(2,1,2)
        %         yLims = [-1 1];
        errorbar(meanIdentArray_hum',meanIdentArray_err_hum','LineWidth',erBr_Lw);
        line(xlim,[0 0],'Color','k');
        line([1 1],[-1 1],'Color','k');
        %         ylim(yLims);
        % text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) ' cells)'],...
        %   'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
        yTicks                                                  = min([min(ylim) -.05]):diff([0 max(ylim)])/4:max(ylim);
        set(gca,'XLim',[0 7],'YLim',[-1 1],'XTick',1:6,'XTickLabel',monkey_data.recombMorphs,...
            'LineWidth',ax_Lw,'Fontsize',ax_Fsz,'Fontweight','bold',...
            'YTick',[yTicks(yTicks<0) 0 yTicks(yTicks>0)]); % diff(yLims)/5:yLims(2)
        
        xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
        ylabel('normalized response','Fontsize',ax_Fsz,'Fontweight','bold')
        title('rank ordered resp. - human faces','Fontsize',tTl_Fsz,'Fontweight','bold')
        
        % plot super title
        foo                                                     = [monkeys{mo} ' (inclCrit#' num2str(g)  ': n=' cellCount{mo} ')'];
        sH                                                      = suptitle(foo);
        set(sH,'Fontsize',sPtTl_Fsz);

        savename                                                = fullfile(popResults_dir,['popRankAveIdentTraj_'...
            monkeys{mo} '_crit' num2str(g) suffix.save]);
        saveas(2323,[savename '.eps'],'epsc');
        export_fig(2323,[savename '.png'],'-nocrop',['-r' num2str(prs)]);
        disp(['saved: ' savename])
        close(2323);
        
        % unranked averages- monkey
        figure(4242); hold on
        set(4242,'Color','w');
        subplot(2,1,2)
        %         yLims = [-1 1];
        errorbar(meanIdentArrayAve_hum',meanIdentArrayAve_err_hum','LineWidth',erBr_Lw);
        line([1 1],ylim,'Color','k');
        %         ylim(yLims);
        % text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) ' cells)'],...
        %   'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
%         yTicks                                                  = min(ylim -.05):diff([min(ylim -.05) max(ylim)])/4:max(ylim);
        set(gca,'XLim',[0 7],'YLim',ylim,'XTick',1:6,'XTickLabel',monkey_data.recombMorphs,...
            'LineWidth',ax_Lw,'Fontsize',ax_Fsz,'Fontweight','bold');
        xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
        ylabel('ave. response(Hz)','Fontsize',ax_Fsz,'Fontweight','bold')
        title('ave. resp. - human faces','Fontsize',tTl_Fsz,'Fontweight','bold')

        % plot super title
        foo                                                     = [monkeys{mo} ' (inclCrit#' num2str(g)  ': n=' cellCount{mo} ')'];
        sH                                                      = suptitle(foo);
        set(sH,'Fontsize',sPtTl_Fsz);

        savename                                                = fullfile(popResults_dir,['popAveIdentTraj_'...
            monkeys{mo} '_crit' num2str(g) suffix.save]);
        saveas(4242,[savename '.eps'],'epsc');
        export_fig(4242,[savename '.png'],'-nocrop',['-r' num2str(prs)]);
        disp(['saved: ' savename])
        close(4242);
        
        % plot population tuning curve
        figure(3030); hold on
        subplot(length(monkeys),length(monkeys)+1,positions)
        errorbar(mean(model_series_hum),steArray,'LineWidth',erBr_Lw,'Color',monkey_colors(mo,:)); hold on
        line(xlim,[0 0],'Color','k');
        yTicks                                                  = [min([min(ylim) -.05]):diff([0 max(ylim)])/4:max(ylim)];
        set(gca,'XLim',[0 7],'XTick',1:6,'XTickLabel',monkey_data.recombMorphs,...
            'LineWidth',ax_Lw,'Fontsize',ax_Fsz,'Fontweight','bold',...
            'YTick',[yTicks(yTicks<0) 0 yTicks(yTicks>0)]); % diff(yLims)/5:yLims(2)
        xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
        ylabel('normalized response','Fontsize',ax_Fsz,'Fontweight','bold')
        title('rank ordered resp. - human faces','Fontsize',tTl_Fsz,'Fontweight','bold')
        
        % plot histogram of tuning curve slopes
        subplot(length(monkeys),length(monkeys)+1,mo*(length(monkeys)+1))
        hist(slopes,30);
        hP                                                      = findobj(gca,'Type','patch');
        set(hP,'FaceColor',monkey_colors(mo,:));
        axis tight;
        line([0 0],ylim,'Color','k');
        title([monkeys{mo} '(n = ' cellCount{mo} ')'],'Fontsize',8,'Fontweight','bold');
        xlim([-0.012 0.012]);
        set(gca,'YTick',max(ylim),'Fontsize',8,'Fontweight','bold','LineWidth',ax_Lw);
        if mo==length(monkeys)
            set(gca,'XTick',[-0.012 0 0.012]);
        else
            set(gca,'XTick',[]);
        end
        
        %         M = cat(dim+1,monPopNormResp{indLocs});        % convert to (dim+1)-dimensional matrix
        %         meanArray = mean(M,dim+1);  % get mean across arrays
        %         steArray = std(M,[],dim+1)./sqrt(length(M));  % get mean across arrays
        %
        %         figure
        %         subplot(211); hold on
        %         errorbar(meanArray',steArray','LineWidth',erBr_Lw);
        %         xLims = xlim;
        %         yLims = ylim;
        %         line([1 1],yLims,'Color','k');
        %         text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) ' cells)'],...
        %             'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
        %         set(gca,'XTick',1:6,'XTickLabel',recombMorphs,'LineWidth',ax_Lw,...
        %             'YTick',[yLims(1):diff(yLims)/5:yLims(2)],'Fontsize',ax_Fsz,'Fontweight','bold')
        %         ylabel('normalized response','Fontsize',ax_Fsz,'Fontweight','bold')
        %         title('monkey face population responses','Fontsize',tTl_Fsz,'Fontweight','bold')
        %
        %         M = cat(dim+1,monPopRankResp{indLocs});        % convert to (dim+1)-dimensional matrix
        %         meanArray = mean(M,dim+1);  % get mean across arrays
        %         steArray = std(M,[],dim+1)./sqrt(length(M));  % get mean across arrays
        %
        %         % get slopes
        %         slopes = [];
        %         y0 = recombMorphs;
        %         for s = 1:length(indLocs)
        %             for e = 1:length(monPopRankResp{indLocs(s)}(:,1))
        %                 P = polyfit(y0,monPopRankResp{indLocs(s)}(e,:),1); %[slope intercept]
        %                 slopes = [slopes; P(1)];
        %             end
        %         end
        %
        %         yLims = [-1 1];
        %         subplot(212); hold on
        %         line(xLims,[0 0],'Color','k');
        %         line([1 1],yLims,'Color','k');
        %         errorbar(meanArray',steArray','LineWidth',erBr_Lw); hold on
        %         ylim(yLims);
        %         text(1.5,yLims(1)+diff(yLims)*.75,['M1 (n = ' num2str(sum(inclIndices)) ' cells)'],...
        %             'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
        %         set(gca,'XTick',1:6,'XTickLabel',recombMorphs,'LineWidth',ax_Lw,...
        %             'YTick',[yLims(1):diff(yLims)/5:yLims(2)],'Fontsize',ax_Fsz,'Fontweight','bold')
        %         xlabel('identity level (%)','Fontsize',ax_Fsz,'Fontweight','bold')
        %         ylabel('normalized response','Fontsize',ax_Fsz,'Fontweight','bold')
        %         title('rank ordered monkey face responses','Fontsize',tTl_Fsz,'Fontweight','bold')
        %
        %         set(gcf,'Position',get(gcf,'Position').*[1 .25 1 2],'Color','w')
        %
        %         axes('Position',get(gca,'Position').*[2.25 1.25 .25 .2]); hold on
        %         box on
        %         hist(slopes,30); axis tight
        %         yTk = get(gca,'YTick');
        %         set(gca,'YTick',yTk(end),'Fontsize',8,'Fontweight','bold','LineWidth',ax_Lw)
        %         line([0 0],ylim,'Color','k')
        %         title('slopes','Fontsize',8,'Fontweight','bold')
        %         endPoint = ceil(max(abs(xlim))*1000)./1000;
        %         xlim([-endPoint endPoint])
        
    end
    
    set(2525,'Position',[500 600 500 500],'Color','w');
    set(3030,'Position',[500 100 500 500],'Color','w');
    
    % do some fancy saving here
    savename                                                    = fullfile(popResults_dir,['popAveIdentTrajMonkey'...
        '_crit' num2str(g) suffix.save]);
    
    saveas(2525,[savename '.eps'],'epsc');
    export_fig(2525,[savename '.png'],'-nocrop',['-r' num2str(prs)]);
    disp(['saved: ' savename])
    close(2525);
    
    savename                                                    = fullfile(popResults_dir,['popAveIdentTrajHuman'...
        '_crit' num2str(g) suffix.save]);
    saveas(3030,[savename '.eps'],'epsc');
    export_fig(3030,[savename '.png'],'-nocrop',['-r' num2str(prs)]);
    disp(['saved: ' savename])
    close(3030);
    
    %     % longitudinal tuning population plot
    %     indLongLocs = setdiff(indLocs,nonContigCellInds);
    %     % dailyMean = cell(length(popTuningCurves{1}),1);
    %
    %     days = cell(length(cellDateArray),1);
    %     for i = 1:length(cellDateArray) % days
    %         for ii = 1:length(cellDateArray{i}) % neurons
    %             days{i} = [days{i}; cellDateArray{i}{ii}];
    %         end
    %     end
    %
    % %     days = cell(length(cellDateArray),1);
    % %     for i = 1:length(cellDateArray) % days
    % %         for ii = 1:length(indLongLocs) % neurons
    % %             days{i} = [days{i}; cellDateArray{indLongLocs(ii)}{i}];
    % %         end
    % %     end
    %
    % %     figure
    % %     for i = 1:length(days)
    % %         tr = ceil(sqrt(length(days)));
    % %         subplot(tr,tr,i)
    % %         imagesc(days{i})
    % %     end
    %
    %     colmap = cool(numDays);
    %     figure; hold on
    %     for i = 1:numDays
    %         stagAx = [1:length(recomb)]+i*.005;
    %         dailyMean = mean(days{i});
    %         dailyErr = std(days{i})./sqrt(length(days{i}));
    %         recombDailyMean = nan(1,length(recomb));
    %         recombDailyErr = nan(1,length(recomb));
    %         for k = 1:length(recombDailyMean)
    %             recombDailyMean(k) = mean(dailyMean(recomb{k}));
    %             recombDailyErr(k) = mean(dailyErr(recomb{k}));
    %         end
    %
    %         errorbar(stagAx,recombDailyMean,recombDailyErr,...
    %             'Color',colmap(i,:))
    %     end
    %     legend('day 1','day 2','day 3','day 4','day 5','day 6','day 7',...
    %         'Location','Northeastoutside')
    %
    %     yLims = ylim;
    %     text(2.5,yLims(2)-diff(yLims)*.2,['M1 (n = ' num2str(length(indLongLocs)) ')'],'EdgeColor','k')
    %     set(gca,'XTick',1:nMorphs+1,'XTickLabel',[0 morphs])
    %     xlabel('identity level')
    %     ylabel('ave. normalized response')
    %     title('longitudinal selectivity for faces')
    %     line([1 1],yLims,'Color','k');
    %     set(gca,'YLim',yLims)
    %
    %     set(gcf,'Color','w')
    %
    %     savename = [popResults_dir exp_code 'PopAveIdentTrajLong' num2str(g) suffix '.png'];
    %     export_fig(gcf,savename,'-nocrop',['-r' num2str(prs)])
    %     close(gcf)
    
end













%%
%
% for g = 1:length(humNoms)
%
%     humPopTemp = humPopResp;
%     humPopTempAdapt = humPopAdapt;
%     humPopTempAdaptFirst = humPopAdaptFirst;
%     monPopTemp = monPopResp;
%     monPopTempAdapt = monPopAdapt;
%     monPopTempAdaptFirst = monPopAdaptFirst;
%
%     % impose each inclusion criteria
%     inclCrit = sum(inclCritList{g},2);
%     cutout = find(inclCrit==0);
%     humPopTemp(cutout) = [];
%     humPopTempAdapt(cutout) = [];
%     humPopTempAdaptFirst(cutout) = [];
%     monPopTemp(cutout) = [];
%     monPopTempAdapt(cutout) = [];
%     monPopTempAdaptFirst(cutout) = [];
%
%     % remove inhibitory responses
%     co = zeros(length(humPopTemp),1);
%     for r = 1:length(humPopTemp)
%         for rr = 1:length(humPopTemp{r}(:,1))
%
%             if humPopTemp{r}(rr,:)<0
%                 co(r)= 1;
%             end
%         end
%     end
%     cutout = find(co==1);
%     humPopTemp(cutout) = [];
%     humPopTempAdapt(cutout) = [];
%     humPopTempAdaptFirst(cutout) = [];
%     monPopTemp(cutout) = [];
%     monPopTempAdapt(cutout) = [];
%     monPopTempAdaptFirst(cutout) = [];
%
%     % normalize responses for each cell
%     for q = 1:length(humPopTemp);
%         nMin = min(min(humPopTemp{q}));
%         humPopTemp{q} = humPopTemp{q}-nMin;
%         nMax = max(max(humPopTemp{q}));
%         humPopTemp{q} = humPopTemp{q}./nMax;
%
%         if sum(sum(humPopTempAdapt{q}))~=0
%             nMin = min(min(humPopTempAdapt{q}));
%             humPopTempAdapt{q} = humPopTempAdapt{q}-nMin;
%             nMax = max(max(humPopTempAdapt{q}));
%             humPopTempAdapt{q} = humPopTempAdapt{q}./nMax;
%         end
%
%         if sum(sum(humPopTempAdaptFirst{q}))~=0
%             nMin = min(min(humPopTempAdaptFirst{q}));
%             humPopTempAdaptFirst{q} = humPopTempAdaptFirst{q}-nMin;
%             nMax = max(max(humPopTempAdaptFirst{q}));
%             humPopTempAdaptFirst{q} = humPopTempAdaptFirst{q}./nMax;
%         end
%     end
%
%     for q = 1:length(monPopTemp);
%         nMin = min(min(monPopTemp{q}));
%         monPopTemp{q} = monPopTemp{q}-nMin;
%         nMax = max(max(monPopTemp{q}));
%         monPopTemp{q} = monPopTemp{q}./nMax;
%
%         if sum(sum(monPopTempAdapt{q}))~=0
%             nMin = min(min(monPopTempAdapt{q}));
%             monPopTempAdapt{q} = monPopTempAdapt{q}-nMin;
%             nMax = max(max(monPopTempAdapt{q}));
%             monPopTempAdapt{q} = monPopTempAdapt{q}./nMax;
%         end
%
%         if sum(sum(monPopTempAdaptFirst{q}))~=0
%             nMin = min(min(monPopTempAdaptFirst{q}));
%             monPopTempAdaptFirst{q} = monPopTempAdaptFirst{q}-nMin;
%             nMax = max(max(monPopTempAdaptFirst{q}));
%             monPopTempAdaptFirst{q} = monPopTempAdaptFirst{q}./nMax;
%         end
%     end
%
%     % calc mean and std of population resp
%     dim = ndims(humPopTemp{1});
%     M = cat(dim+1,humPopTemp{:});
%     rawPop{1} = M;
%     meanPop{1} = mean(M,dim+1);
%     stdPop{1} = std(M,0,dim+1);
%     temp = std(M,0,dim+1);
%     temp2 = size(temp);
%     stePop{1} = temp/sqrt(temp2(2));
%     dim = ndims(humPopTempAdapt{1});
%     M = cat(dim+1,humPopTempAdapt{:});
%     meanPop{2} = mean(M,dim+1);
%     dim = ndims(humPopTempAdaptFirst{1});
%     M = cat(dim+1,humPopTempAdaptFirst{:});
%     meanPop{3} = mean(M,dim+1);
%
%     dim = ndims(monPopTemp{1});
%     M = cat(dim+1,monPopTemp{:});
%     rawPop{2} = M;
%     meanPop{4} = mean(M,dim+1);
%     stdPop{2} = std(M,0,dim+1);
%     temp = std(M,0,dim+1);
%     temp2 = size(temp);
%     stePop{2} = temp/sqrt(temp2(2));
%     dim = ndims(monPopTempAdapt{1});
%     M = cat(dim+1,monPopTempAdapt{:});
%     meanPop{5} = mean(M,dim+1);
%     dim = ndims(monPopTempAdaptFirst{1});
%     M = cat(dim+1,monPopTempAdaptFirst{:});
%     meanPop{6} = mean(M,dim+1);
%
%     dirStr = '/projects/apjones/results/faceLearning/population/identTrajAves/'; % directory for saving plots
%
%     %% plot population dodecs
%     aveWvFrm = squeeze(dat.wf.v(1,:,:));
%
%     humPopRecombDodec = nan(size(humRecomb));
%     monPopRecombDodec = nan(size(monRecomb));
%
%     humPopRecombDodec(:,1) = meanPop{1}(:,1);
%     humPopRecombDodec(:,2) = mean(meanPop{1}(:,2:4),2);
%     humPopRecombDodec(:,3) = mean(meanPop{1}(:,5:6),2);
%     humPopRecombDodec(:,4) = meanPop{1}(:,7);
%     humPopRecombDodec(:,5) = meanPop{1}(:,8);
%     humPopRecombDodec(:,6) = meanPop{1}(:,9);
%
%     humPopRecombAdapt(:,1) = meanPop{2}(:,1);
%     humPopRecombAdapt(:,2) = mean(meanPop{2}(:,2:4),2);
%     humPopRecombAdapt(:,3) = mean(meanPop{2}(:,5:6),2);
%     humPopRecombAdapt(:,4) = meanPop{2}(:,7);
%     humPopRecombAdapt(:,5) = meanPop{2}(:,8);
%     humPopRecombAdapt(:,6) = meanPop{2}(:,9);
%
%     humPopRecombAdaptFirst(:,1) = meanPop{3}(:,1);
%     humPopRecombAdaptFirst(:,2) = mean(meanPop{3}(:,2:4),2);
%     humPopRecombAdaptFirst(:,3) = mean(meanPop{3}(:,5:6),2);
%     humPopRecombAdaptFirst(:,4) = meanPop{3}(:,7);
%     humPopRecombAdaptFirst(:,5) = meanPop{3}(:,8);
%     humPopRecombAdaptFirst(:,6) = meanPop{3}(:,9);
%
%     monPopRecombDodec(:,1) = meanPop{4}(:,1);
%     monPopRecombDodec(:,2) = mean(meanPop{4}(:,2:4),2);
%     monPopRecombDodec(:,3) = mean(meanPop{4}(:,5:6),2);
%     monPopRecombDodec(:,4) = meanPop{4}(:,7);
%     monPopRecombDodec(:,5) = meanPop{4}(:,8);
%     monPopRecombDodec(:,6) = meanPop{4}(:,9);
%
%     monPopRecombAdapt(:,1) = meanPop{5}(:,1);
%     monPopRecombAdapt(:,2) = mean(meanPop{5}(:,2:4),2);
%     monPopRecombAdapt(:,3) = mean(meanPop{5}(:,5:6),2);
%     monPopRecombAdapt(:,4) = meanPop{5}(:,7);
%     monPopRecombAdapt(:,5) = meanPop{5}(:,8);
%     monPopRecombAdapt(:,6) = meanPop{5}(:,9);
%
%     monPopRecombAdaptFirst(:,1) = meanPop{6}(:,1);
%     monPopRecombAdaptFirst(:,2) = mean(meanPop{6}(:,2:4),2);
%     monPopRecombAdaptFirst(:,3) = mean(meanPop{6}(:,5:6),2);
%     monPopRecombAdaptFirst(:,4) = meanPop{6}(:,7);
%     monPopRecombAdaptFirst(:,5) = meanPop{6}(:,8);
%     monPopRecombAdaptFirst(:,6) = meanPop{6}(:,9);
%
%     plot_spheres_noFaces(humPopRecombDodec,monPopRecombDodec,visTog)
%     xLims = xlim;
%     yLims = ylim;
%     text(xLims(1)*2.75,yLims(1)*1.25,['M1 (n = ' num2str(length(humPopTemp)) 'cells)'],...
%         'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
%
%     figure(37)
%     tempPos = get(37,'Position');
%     set(37,'PaperUnits','inches','PaperPosition',tempPos/prs);
%     savename = [dodecStr exp_code 'popAveDodecTorCrit' num2str(g) suffix '.png'];
%     export_fig(37,savename,'-nocrop',['-r' num2str(prs)])
%     clf(37)
%
%     plot_spheres_single_noFace(monPopRecombDodec,visTog)
%     text(xLims(1)*2.75,yLims(1)*1.25,['M1 (n = ' num2str(length(humPopTemp)) 'cells)'],...
%         'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
%     %     tH = suplabel(['population ave-' num2str(length(humPopTemp)) ' cells'],'t');
%     %     set(tH,'Fontweight','Bold','Fontsize',10)
%
%     figure(37)
%     tempPos = get(37,'Position');
%     set(37,'PaperUnits','inches','PaperPosition',tempPos/prs);
%     savename = [dodecStr exp_code 'popAveDodecMonksTorCrit' num2str(g) suffix '.png'];
%     export_fig(37,savename,'-nocrop',['-r' num2str(prs)])
%     clf(37)
%
%     plot_spheres_noFaces(humPopRecombAdapt,monPopRecombAdapt,visTog)
%     xLims = xlim;
%     yLims = ylim;
%     text(xLims(1)*2.75,yLims(1)*1.25,['M1 (n = ' num2str(length(humPopTempAdapt)) 'cells)'],...
%         'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
%
%     figure(37)
%     tempPos = get(37,'Position');
%     set(37,'PaperUnits','inches','PaperPosition',tempPos/prs);
%     savename = [dodecStr exp_code 'popAveDodecTorCrit' num2str(g) 'Adapt.png'];
%     export_fig(37,savename,'-nocrop',['-r' num2str(prs)])
%     clf(37)
%
%     plot_spheres_single_noFace(monPopRecombAdapt,visTog)
%     text(xLims(1)*2.75,yLims(1)*1.25,['M1 (n = ' num2str(length(humPopTemp)) 'cells)'],...
%         'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
%
%     figure(37)
%     tempPos = get(37,'Position');
%     set(37,'PaperUnits','inches','PaperPosition',tempPos/prs);
%     savename = [dodecStr exp_code 'popAveDodecMonksTorCrit' num2str(g) 'Adapt.png'];
%     export_fig(37,savename,'-nocrop',['-r' num2str(prs)])
%     clf(37)
%
%     plot_spheres_noFaces(humPopRecombAdaptFirst,monPopRecombAdaptFirst,visTog)
%     xLims = xlim;
%     yLims = ylim;
%     text(xLims(1)*2.75,yLims(1)*1.25,['M1 (n = ' num2str(length(humPopTempAdaptFirst)) 'cells)'],...
%         'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
%
%     figure(37)
%     tempPos = get(37,'Position');
%     set(37,'PaperUnits','inches','PaperPosition',tempPos/prs);
%     savename = [dodecStr exp_code 'popAveDodecTorCrit' num2str(g) 'AdaptFirst.png'];
%     export_fig(37,savename,'-nocrop',['-r' num2str(prs)])
%     clf(37)
%
%     plot_spheres_single_noFace(monPopRecombAdaptFirst,visTog)
%     text(xLims(1)*2.75,yLims(1)*1.25,['M1 (n = ' num2str(length(humPopTempAdapt)) 'cells)'],...
%         'Fontsize',ax_Fsz,'Fontweight','bold','EdgeColor','k','LineWidth',ax_Lw)
%
%     figure(37)
%     tempPos = get(37,'Position');
%     set(37,'PaperUnits','inches','PaperPosition',tempPos/prs);
%     savename = [dodecStr exp_code 'popAveDodecMonksTorCrit' num2str(g) 'AdaptFirst.png'];
%     export_fig(37,savename,'-nocrop',['-r' num2str(prs)])
%     clf(37)
%
%     axLims = axis;
%     hold all
%
%     %     if g==1
%     %
%     %         % create movie for average dodec
%     %         winSize = get(37,'Position');
%     %         winSize(1:2) = [0 0];
%     %         numFrames = 180;
%     %         A = moviein(numFrames,37,winSize);
%     %         % set(37,'NextPlot','replacechildren')
%     %
%     %         for b = 1:numFrames
%     %             subplot(1,2,1)
%     %             view(b,15)
%     %             set(gca,'XLim',axLims(1:2),'YLim',axLims(3:4),'ZLim',axLims(5:6));
%     %
%     %             subplot(1,2,2)
%     %             view(b,15)
%     %             set(gca,'XLim',axLims(1:2),'YLim',axLims(3:4),'ZLim',axLims(5:6));
%     %
%     %             % add to movie matrix
%     %             A(:,b) = getframe(37,winSize);
%     %             % add to gif
%     %             im = frame2im(A(:,b));
%     %             [imind,cm] = rgb2ind(im,256);
%     %             if b==1;
%     %                 imwrite(imind,cm,[movDir 'FSgifPopMovieTor.gif'],'gif',...
%     %                     'Loopcount',inf,'DelayTime',.05,'TransparentColor',1);
%     %             else
%     %                 imwrite(imind,cm,[movDir 'FSgifPopMovieTor.gif'],'gif',...
%     %                     'WriteMode','Append','DelayTime',.05,'TransparentColor',1);
%     %             end
%     %         end
%     %
%     %         movie(37,A,3,9,winSize)
%     %         %                                                 save (fullfile(movDir,['/movie' tempName 'Movie.mat']) ,'A')
%     %         movie2avi(A,fullfile(movDir,'/popAveTorMov.avi'))
%     %     end
%     clf(37)
%
%     %      % load rhombus pop dodec variables
%     %     if g==1
%     %         savename = [dataStr  'popMonDodecrhombus.mat'];
%     %         rhombPopMonDodec = load(savename);
%     %         keyboard
%     %     end
%     %%
%
%     % {human monkey} identity X morph X cell
%     for v = 1:2
%
%         % plot simple identity trajectory-tuning function
%         figure(532+v);clf
%         set(gcf,'color','w','Visible','off'); % ,'renderer','zbuffer'
%         temp = get(gcf,'Position');
%         temp([1 2 4]) = [0 0 fullscr(4)*.9];
%         set(532+v,'Position',temp)
%         if v == 1
%             nom = humNoms{g};
%             rk = 0;
%         elseif v == 2
%             nom = monNoms{g};
%             rk = 12;
%         end
%
%         tetPileImList(12,end) = 96;
%
%         for q = 1:length(meanPop{v}(:,1))
%
%             subplot(6,4,4+q+floor(q/4.1)*4)
%             disp(4+q+floor(q/4.1)*4)
%             errorbar([0 morphs],meanPop{v}(q,:),stePop{v}(q,:),'LineWidth',tclw); hold on
%             axis tight
%             doop = [min(min(meanPop{v}))-max(max(stePop{v})) max(max(meanPop{v}))+max(max(stePop{v}))];
%             ylim(round(doop.*10)./10)
%
%             temp = get(gca,'Position');
%             temp(4) = temp(4)*.75;
%             set(gca,'Position',temp,'XTick',[0 morphs]);
%
%             if q==1
%                 set(gca,'XTickLabel',[],'YTick',round(doop.*10)./10);
%                 ylabel('normFR','FontWeight','bold','FontSize',ax_Fsz);
%             elseif q==12
%                 set(gca,'XTickLabel',{0,[],[],[],[],[],25,50,100},'YTick',[]);
%                 %set(gca,'YTickLabel',[]);
%                 xlabel('% identity','FontWeight','bold','FontSize',ax_Fsz);
%             else
%                 set(gca,'XTickLabel',[],'YTick',[]);
%                 %set(gca,'YTickLabel',[]);
%             end
%
%             % plot corresponding image
%             subplot(6,4,q+floor(q/4.1)*4)
%             b2 = [temp(1)+temp(3)/4 temp(2)+temp(4) 0.07 0.07];
%             %
%             if v == 1
%                 image(stim_images{tetPileImList(q+rk,end)})
%             elseif v == 2
%                 image(stim_images{96+(q*8)})
%             end
%
%             set(gca,'Position',b2);
%             axis off; axis image;
%             if q==1
%                 d = title([nom ' (n = ' num2str(sum(inclCrit)) '/' num2str(length(inclCrit)) ')'],'FontWeight','bold');
%                 set(d,'interpreter','none');
%             end
%
%             set(gcf,'Visible','off');
%         end
%
%         savename = [dirStr exp_code nom 'norm.png'];
%         export_fig(532+v,savename,'-nocrop',['-r' num2str(prs)])
%         close(532+v)
%
%         % plot simple identity trajectory-tuning function
%         figure(542+v);clf
%         temp = get(gcf,'Position');
%         temp([1 2 4]) = [0 0 fullscr(4)*.9];
%         set(gcf,'color','w','Position',temp,'Visible','off') % ,'renderer','zbuffer'
%         if v == 1
%             nom = humNoms{g};
%             rk = 0;
%         elseif v == 2
%             nom = monNoms{g};
%             rk = 12;
%         end
%
%         for q = 1:length(meanPop{v}(:,1))
%             subplot(6,4,4+q+floor(q/4.1)*4)
%             plot([0 morphs],squeeze(rawPop{v}(q,:,:)),'LineWidth',tclw); hold on
%             %                     errorbar([0 morphs],meanPop{v}(q,:),stePop{v}(q,:),'LineWidth',tclw); hold on
%             axis tight
%             ylim([0 1])
%             temp = get(gca,'Position');
%             temp(4) = temp(4)*.75;
%             set(gca,'Position',temp,'XTick',[0 morphs]);
%
%             if q==1
%                 set(gca,'XTickLabel',[]);
%                 ylabel('normFR','FontWeight','bold','FontSize',ax_Fsz);
%             elseif q==12
%                 set(gca,'XTickLabel',{0,[],[],[],[],[],25,50,100});
%                 %set(gca,'YTickLabel',[]);
%                 xlabel('% identity','FontWeight','bold','FontSize',ax_Fsz);
%             else
%                 set(gca,'XTickLabel',[]);
%                 %set(gca,'YTickLabel',[]);
%             end
%
%             % plot corresponding image
%             subplot(6,4,q+floor(q/4.1)*4)
%             b2 = [temp(1)+temp(3)/4 temp(2)+temp(4) 0.07 0.07];
%             image(stim_images{tetPileImList(q+rk,end)})
%             set(gca,'Position',b2);
%             axis off; axis image;
%             if q==1
%                 d = title([nom ' (n = ' num2str(sum(inclCrit)) '/' num2str(length(inclCrit)) ')'],'FontWeight','bold');
%                 set(d,'interpreter','none');
%             end
%
%             set(gcf,'Visible','off');
%         end
%
%         savename = [dirStr exp_code 'PopNorm.png'];
%         export_fig(542+v,savename,'-nocrop',['-r' num2str(prs)])
%         close(542+v)
%     end
% end

% close all
end