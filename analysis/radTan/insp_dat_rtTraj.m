function insp_dat_rtTraj(paths,monkey,suffix)
% loop through each '*rt' file in 'paths.mas'
%
% Usage: insp_dat_rtTraj(paths,monkey,suffix)
% MONKEY is the monkey's name
% PATHS project path structure containing at minimum paths.rare field
% SUFFIX is a string appended to the end of each filename
%
% NOTE: stimset includes 24 humans and 12 monkeys - 12 hum & 6 mon face pairs (3 trajectories[2 rad;1 tan])
%
% last modified 6-18-16
% apj


% set project identifier
proj                                            = 3;

% define stimuli
[stimNames,~]                                   = getStimList(proj,fullfile(paths.physio,'ptb'));

% load appropriate image files for display
load(fullfile(paths.stim,'imgRadTanSetWhite.mat'));

% set experimental variables
visTog                                          = 'on';
figsOn                                          = 0;                        % switch for figures (on vs. off)
prs                                             = 150;                      % print resolution
figJit                                          = round(rand*10);           % avoid using the same fig number in parallel
ax_Fsz                                          = 8;                        % axes fontsize
title_Fsz                                       = 10;                       % title fontsize
spk_Fsz                                         = 10;                       % spike-name fontsize
popTC_Fsz                                       = 10;                       % lineWidth for sdf's in population plot
stim_Lw                                         = 4.5;                      % Stimulus presentation linewidth
ax_Lw                                           = 1;                        % axes linewidth
sdf_Lw                                          = 2;                        % lineWidth for sdf's
popTC_Lw                                        = 4;                        % lineWidth for sdf's in population plot
fullscr                                         = get(0,'ScreenSize');      % get screen size
figPos1                                         = [fullscr(3)/3   fullscr(4)/2   fullscr(3)/2   fullscr(4)/4]; % [x y width height]
figPos2                                         = [fullscr(3)/5   fullscr(4)/6   fullscr(3)/5   fullscr(4)/2]; % [x y width height]
if regexp(visTog,'off')
    figPos1(1)                                  = figPos1(1)+5000;
    figPos2(1)                                  = figPos2(1)+5000;
end


Plotstart                                       = -300;                     % plot display size (msec)
Plotstop                                        = 600;
Respwin                                         = [100 400];                % time window for evoked rate

% face numbers
humFaceNums                                     = 1:12;
monFaceNums                                     = 1:6;
% faceList                                        = [humFaceNums monFaceNums];
humNorm                                         = 535;
monNorm                                         = 536;

% get file list
exp_code                                        = [monkey(1) paths.code];
fileList                                        = dir(fullfile(paths.mas,[exp_code '*' suffix.load '.mat']));

% get monkey index for saving within population data
monkey_positions                                = 'trms';
monkeyNum                                       = strfind(monkey_positions,monkey(1));

% pre-allocate
popRadHeat                                      = cell(1,length(fileList));
popTanHeat                                      = cell(1,length(fileList));
% popRadAve                                       = {cell(1,length(fileList)); cell(1,length(fileList))};
% popTanAve                                       = {cell(1,length(fileList)); cell(1,length(fileList))};
popRadSlopesHum                                    = cell(1,length(fileList));
popRadSlopesMon                                    = cell(1,length(fileList));
% popTanSlopes                                    = cell(1,length(fileList));
popRadSlopesNoiseHum                               = cell(1,length(fileList));
popRadSlopesNoiseMon                               = cell(1,length(fileList));
% popTanSlopesNoise                               = cell(1,length(fileList));

% loop through spike-files selected above
for i = 1:length(fileList);

    % load dat file
    filename                                    = fullfile(paths.mas,fileList(i).name);
%     filename                                    = fullfile(paths.mas,[fileList(i).name(1:end-9) '.mat']);
    load(filename);
    
    disp(filename);
    unpack;
    
    % read variables
    % dates                                       = datenum(dat.h.date,'dd-mmm-yy');
    radList                                     = [dat.h.stim.radSteps(dat.h.stim.radSteps<-10) 0 ...
        dat.h.stim.radSteps(dat.h.stim.radSteps>10&dat.h.stim.radSteps<=200)];
    tanList                                     = [0 dat.h.stim.tanSteps 100];
    stim_win                                    = round(dat.c(1,T_STIMON:T_STIMOFF)); % stimulus onset/offset (msec)
    
    if figsOn==1
        % open figure
        figure('Visible',visTog,'Position',figPos1,'Color','w');
    end

    
    %% species loop
    for sP = 1:2
        
        % set variables specific to each species presented
        if sP==1
            species.stimOffset                  = 0;
            species.parseOffset                 = 0;
            species.norm                        = humNorm;
            species.facePairs                   = humFaceNums;
        elseif sP==2
            species.stimOffset                  = 348;
            species.parseOffset                 = 12;
            species.norm                        = monNorm;
            species.facePairs                   = monFaceNums;
        end
        
        %% facepairs loop
        % pre-allocate for next loop
        radCellHeat                             = nan(length(species.facePairs)*2,length(radList));
        tanCellHeat                             = nan(length(species.facePairs),length(tanList));
        for fP = 1:length(species.facePairs);
            
            % define each stimulus number in this trajectory
            tempMorphs                          = ([1:length(radList)+1]+((fP-1)*29))+species.stimOffset;
            radTraj{1}                          = [tempMorphs(1:3) species.norm tempMorphs(6:end)]; % re-order so norm (0%) is in appropriate spot
            tempMorphs                          = ([1:length(radList)+1]+((fP-1)*29))+length(radList)+2+species.stimOffset;
            radTraj{2}                          = [tempMorphs(1:3) species.norm tempMorphs(6:end)];
            tanTraj                             = [radTraj{2}(8) [25:(25+length(tanList)-3)]+((fP-1)*29)+species.stimOffset radTraj{1}(8)];
            
            %% rad loop
            % pre-allocate for next loop
            radAxList                           = nan(2,length(radTraj{1}));
            for sW = 1:length(radTraj)
                stimNumTraj                     = radTraj{sW};
                
                for radsT = 1:length(stimNumTraj)
                    
                    subdat                      = select_trials(dat,STIM,stimNumTraj(radsT));
                    
                    % read sdf
                    if ~isempty(subdat.s)
                        [y,~]                   = get_sdf(subdat,1,[Plotstart Plotstop]);
                    else % provide some dummy variables
                        y                       = zeros(1,1+diff([Plotstart Plotstop])); % e1 = 0;
                    end
                    
                    % count spikes during response window
                    [r,~]                       = find(subdat.s>Respwin(1)&subdat.s<Respwin(2));
                    trialSpikes                 = histc(r,1:length(subdat.s(:,1)));
                    [r,~]                       = find(subdat.s>Plotstart&subdat.s<0);
                    trialBackGround             = histc(r,1:length(subdat.s(:,1)));
                    eResp                       = mean(trialSpikes)*(1000/diff(Respwin));
                    if monkeyNum==1
                        bResp                   = 0;
                    else
                        bResp                   = mean(trialBackGround)*(1000/abs(Plotstart));
                    end
                    
                    radCellHeat((fP-1)*2+sW,radsT)   = eResp-bResp;
                    % radCellHeat(fP,radsT) = mean(y(stim_win(1)+bkground_len+100:stim_win(2)+bkground_len+100));   % response
                    
                    if figsOn==1
                        
                        % plot image
                        subplot(3,length(stimNumTraj),radsT);
                        image(img{stimNumTraj(radsT)});
                        axis off image; hold on
                        
                        % print spike name
                        if radsT==1  
                            text(diff(xlim)*-2,diff(ylim)/2,dat.h.snames,'Color','k','FontSize',spk_Fsz)
                        end
                        
                        % plot sdf
                        subplot(3,length(stimNumTraj),radsT+length(stimNumTraj)); hold on
                        plot([1:length(y)]-abs(Plotstart),y,'-b','LineWidth',sdf_Lw);
                        set(gca,'LineWidth',ax_Lw);
                        axis tight
                        
                        % accumnulate info
                        radAxList(:,radsT)      = round(ylim'*10)/10;
                        sub_pos                 = get(gca,'Position');
                        plot_left               = sub_pos(1);
                        
                        % plot title
                        title([num2str(radList(radsT)) '%'],'Fontsize',title_Fsz)
                        
                        % plot raster
                        subplot(3,length(stimNumTraj),radsT+length(stimNumTraj)*2)
                        plot_rasters(subdat.s,[Plotstart Plotstop]);
                        templims = ylim;
                        yline(stim_win);
                        
                        if radsT==1
                            xtick               = round(stim_win);
                        else
                            xtick               = [];
                        end
                        sub_pos                 = get(gca,'Position');
                        sub_pos(1)              = plot_left;
                        set(gca,'XTick',xtick,'YTick',[],'Fontsize',ax_Fsz,...
                            'Position',sub_pos.*[1 .9 1 1],'Linewidth',ax_Lw);
                        if radsT==1
                            set(gca,'YTick',templims(1),'YTickLabel',length(subdat.s(:,1)));
                        end
                        box on
                        
                    end
                end
                
                if figsOn==1
                    
                    %% normalize axes
                    % get new limits
                    axMin                       = min(radAxList(1,:));
                    axMax                       = max(radAxList(2,:))*1.1;
                    
                    % loop through axes
                    for aX = 1:length(radAxList)
                        
                        subplot(3,length(stimNumTraj),aX+length(stimNumTraj))
                        sub_pos                 = get(gca,'position');
                        
                        % normalize y axis & mark stim presentation time
                        line(stim_win,[axMax axMax],'Color','k','LineWidth',stim_Lw);
                        set(gca,'Position',sub_pos.*[1 1.1 1 1],'Fontsize',ax_Fsz,...
                            'YLim',[axMin axMax],'Linewidth',ax_Lw)
                        box off
                        
                        if aX==1
                            set(gca,'XTick',stim_win,'YTick',ylim)
                            xlabel('time(ms)','FontWeight','bold','FontSize',ax_Fsz);
                            ylabel('fr(Hz)','FontWeight','bold','FontSize',ax_Fsz);
                        else
                            set(gca,'XTick',[],'YTick',[])
                        end
                    end
                    
                    % do some fancy saving here
                    temp                        = fullfile(paths.results,monkey);
                    if ~exist(temp,'dir')
                        mkdir(temp)
                    end
                    savename                    = fullfile(temp,[exp_code '_' dat.h.snames ...
                        '_radStim' stimNames{stimNumTraj(1)}(1:7) '' suffix.save '.png']);
%                 saveas(gcf,[savename '.eps'],'epsc');
                    export_fig(gcf,savename,'-nocrop',['-r' num2str(prs)]);
                    disp(['saved: ' savename])
                    clf(gcf)
                end
            end
            
            
            %% tan loop
            % pre-allocate for next loop
            tanAxList                           = nan(2,length(tanTraj));
            for tansT = 1:length(tanTraj)
                stimNumTraj                     = fliplr(tanTraj);
                
                subdat                          = select_trials(dat,STIM,stimNumTraj(tansT));
                if ~isempty(subdat.s)
                    % get sdf
                    [y,~]                       = get_sdf(subdat,1,[Plotstart Plotstop]);
                else % provide some dummy variables
                    y                           = zeros(1,1+diff([Plotstart Plotstop])); %  e1 = 0;
                end
                
                % count spikes during response window
                [r,~]                           = find(subdat.s>Respwin(1)&subdat.s<Respwin(2));
                trialSpikes                     = histc(r,1:length(subdat.s(:,1)));
                [r,~]                           = find(subdat.s>Plotstart&subdat.s<0);
                trialBackGround                 = histc(r,1:length(subdat.s(:,1)));
                eResp                           = mean(trialSpikes)*(1000/diff(Respwin));
                if monkeyNum==1
                    bResp                       = 0;
                else
                    bResp                       = mean(trialBackGround)*(1000/abs(Plotstart));
                end
                
                tanCellHeat(fP,tansT)           = eResp-bResp;
                % tanCellHeat(fP,tansT)         = mean(y(stim_win(1)+bkground_len+100:stim_win(2)+bkground_len+100));   % response
                
                if figsOn==1
                    
                    % plot image
                    subplot(3,length(stimNumTraj),tansT)
                    image(img{stimNumTraj(tansT)})
                    axis off image; hold on
                    if tansT==1  % print spike name
                        text(diff(xlim)*-2,diff(ylim)/2,dat.h.snames,'Color','k','FontSize',spk_Fsz)
                    end
                    
                    % plot sdf
                    subplot(3,length(stimNumTraj),tansT+length(stimNumTraj)); hold on
                    plot([1:length(y)]-abs(Plotstart),y,'-b','LineWidth',sdf_Lw);
                    set(gca,'LineWidth',ax_Lw);
                    axis tight
                    
                    % accumnulate info
                    tanAxList(:,tansT)          = round(ylim'*10)/10;
                    sub_pos                     = get(gca,'Position');
                    plot_left                   = sub_pos(1);
                    
                    % plot title
                    if tansT==1
                        title('face A','Fontsize',title_Fsz)
                    elseif tansT==length(tanTraj)
                        title('face B','Fontsize',title_Fsz)
                    else
                        title([num2str(tanList(tansT)) '%'],'Fontsize',title_Fsz)
                    end
                    
                    % plot raster
                    subplot(3,length(stimNumTraj),tansT+length(stimNumTraj)*2)
                    plot_rasters(subdat.s,[Plotstart Plotstop]);
                    templims                    = ylim;
                    yline(stim_win);
                    
                    if tansT==1
                        xtick                   = round(stim_win);
                    else
                        xtick                   = [];
                    end
                    sub_pos                     = get(gca,'Position');
                    sub_pos(1)                  = plot_left;
                    set(gca,'XTick',xtick,'YTick',[],'Fontsize',ax_Fsz,...
                        'Position',sub_pos.*[1 .9 1 1],'Linewidth',ax_Lw);
                    if tansT==1
                        set(gca,'YTick',templims(1),'YTickLabel',length(subdat.s(:,1)));
                    end
                    box on
                    
                end
            end
            
            if figsOn==1
                
                %% normalize axes
                % get new limits
                axMin                           = min(tanAxList(1,:));
                axMax                           = max(tanAxList(2,:))*1.1;
                
                % loop through axes
                for aX = 1:length(tanAxList)
                    %                 set(gcf,'Visible',visTog)
                    subplot(3,length(stimNumTraj),aX+length(stimNumTraj))
                    sub_pos                     = get(gca,'position');
                    ylim([axMin axMax])
                    hold on
                    
                    % normalize y axis & mark stim presentation time
                    line(stim_win,[axMax axMax],'Color','k','LineWidth',stim_Lw);
                    set(gca,'Position',sub_pos.*[1 1.1 1 1])
                    box off
                    if aX==1
                        set(gca,'XTick',stim_win,'YTick',ylim)
                        xlabel('time(ms)','FontWeight','bold','FontSize',ax_Fsz);
                        ylabel('fr(Hz)','FontWeight','bold','FontSize',ax_Fsz);
                    else
                        set(gca,'XTick',[],'YTick',[])
                    end
                end
                
                %             % print label
                %             suplabel([monkey dat.h.snames ' ' stimNames{stimNumTraj(tansT)}(1:end-7)],'t');
                
                %             tempPos = get(gcf,'Position');
                %             set(gcf,'Visible',visTog)
                
                % do some fancy saving here
                temp                            = fullfile(paths.results,monkey);
                if ~exist(temp,'dir')
                    mkdir(temp)
                end
                
                tempInd1                        = strfind(stimNames{stimNumTraj(1)},'to');
                tempInd2                        = strfind(stimNames{stimNumTraj(2)},'to');
                savename                        = fullfile(temp,[exp_code '_' dat.h.snames...
                    '_tanStim' stimNames{stimNumTraj(1)}(1:tempInd1-1) '-' stimNames{stimNumTraj(end)}(1:tempInd2-1) suffix.save]);
                
                export_fig(gcf,savename,'-nocrop',['-r' num2str(prs)]);
%                 saveas(gcf,[savename '.eps'],'epsc');
                disp(['saved: ' [savename '.png']])
                clf(gcf)
            end
        end
        
        %% save cell summary
        if sP==1
            popRadHeat{i}                       = radCellHeat; % size(radCellHeat) 12 10 double
            popTanHeat{i}                       = tanCellHeat;
        else
            popRadHeat{i}                       = [popRadHeat{i}; radCellHeat];
            popTanHeat{i}                       = [popTanHeat{i}; tanCellHeat];
        end
        
    end
    
%     if i==length(fileList)
%         keyboard % size(popRadHeat{sP}) 1 50 cells of 12 10 double
%     end
    
    %% radial trajectory plot
    % humans
    humRad                                      = popRadHeat{i}(1:species.parseOffset*2,:);
    humRadAve                                   = mean(humRad);
    humRadErr                                   = std(humRad)./sqrt(length(humRad(:,1)));
    [temp_max,temp_ind]                         = max(humRadAve);     % flip inverted radial TC's
    if temp_ind==find(radList==0)
        flipped                                 = humRadAve.*-1;
        humRadAve                               = flipped+diff([max(flipped) temp_max]);
        disp(['Human TC flipped for spike:' dat.h.snames])
    end
    
    % monkeys
    monRad                                      = popRadHeat{i}(species.parseOffset*2+1:end,:);
    monRadAve                                   = mean(monRad);
    monRadErr                                   = std(monRad)./sqrt(length(monRad(:,1)));
    [temp_max,temp_ind]                         = max(monRadAve);    % flip inverted radial TC's
    if temp_ind==find(radList==0)
        flipped                                 = monRadAve.*-1;
        monRadAve                               = flipped+diff([max(flipped) temp_max]);
        disp(['Monkey TC flipped for spike:' dat.h.snames])
    end
    %     popRadAve{i} = monRadAve;
    
%     if any(any(humRad>5))
%         humSlopeRatios = nan(2,length(humRad(:,1)));
%         humSlopeRatiosNoise = nan(2,length(humRad(:,1)));
%         foo = length(humRad(1,:));
%         for st = 1:length(humRad(:,1))
%             % measure slope in real TC
%             p1=polyfit(1:3,humRad(st,1:3),1);
%             p2=polyfit(5:foo,humRad(st,5:foo),1);
%             humSlopeRatios(1,st) = p1(1);
%             humSlopeRatios(2,st) = p2(1);
%             
%             % measure slope in randomized TC
%             rand_monRad = humRad(st,randperm(foo));
%             p1=polyfit(1:3,rand_monRad(1:3),1);
%             p2=polyfit(5:foo,rand_monRad(5:foo),1);
%             humSlopeRatiosNoise(1,st) = p1(1);
%             humSlopeRatiosNoise(2,st) = p2(1);
%             
%         end
% 
%         popRadSlopesHum{i} = humSlopeRatios;
%         popRadSlopesNoiseHum{i} = humSlopeRatiosNoise;
%     end
%     
%     if any(any(monRad>5))
%         monSlopeRatios = nan(2,length(monRad(:,1)));
%         monSlopeRatiosNoise = nan(2,length(monRad(:,1)));
%         foo = length(monRad(1,:));
%         for st = 1:length(monRad(:,1))
%             
%             % measure slope in real TC           
%             x1 = 1:4;
%             x2 = 4:foo;
%             y1 = monRad(st,x1);
%             y2 = monRad(st,x2);
%             p1 = polyfit(x1,y1,1);
% %             yFit1 = polyval(p1,x1);
%             p2 = polyfit(x2,y2,1);
%             yFit2 = polyval(p2,x2);
%             
%             if sqrt(mean((y2 - yFit2).^2))<1  % Root Mean Squared Error < 1
%                 monSlopeRatios(1,st) = p1(1);
%                 monSlopeRatios(2,st) = p2(1);
%                 
%                 % measure slope in randomized TC
%                 rand_monRad = monRad(st,randperm(foo));
%                 p1 = polyfit(x1,rand_monRad(x1),1);
%                 p2 = polyfit(x2,rand_monRad(x2),1);
%                 monSlopeRatiosNoise(1,st) = p1(1);
%                 monSlopeRatiosNoise(2,st) = p2(1);
%             end
%             
%         end
%         
%         popRadSlopesMon{i} = monSlopeRatios;
%         popRadSlopesNoiseMon{i} = monSlopeRatiosNoise;
%     end
    
    if figsOn==1
        
        figure(31+figJit)
        set(31+figJit,'Color','w','Position',figPos2,'Visible',visTog);
        
        subplot(2,1,1);
        plot(1:length(humRadAve),humRad,'Color',[.9 .9 .9],'LineWidth',popTC_Lw); hold on
        errorbar(humRadAve,humRadErr,'Color',[.1 .1 .1],'LineWidth',popTC_Lw+.5);
        
        ylim([min(min(humRad))-max(max(humRad))*.05 max(max(humRad))+max(max(humRad))*.05]);
        set(gca,'XLim',[0 length(humRadAve)+1],'XTick',1:length(humRadAve),'XTickLabel',radList,...
            'LineWidth',ax_Lw,'FontWeight','Bold','FontSize',popTC_Fsz,'YTick',ylim,'YTickLabel',round(ylim.*10)./10)
        ylabel('ave. FR (Hz)','FontSize',popTC_Fsz);
        title(['human ave.(n = ' num2str(length(humRad(:,1))) ' faces)'],'FontSize',popTC_Fsz);
        
        subplot(2,1,2);
        plot(1:length(monRadAve),monRad,'Color',[.9 .9 .9],'LineWidth',popTC_Lw); hold on
        errorbar(monRadAve,monRadErr,'Color',[.1 .1 .1],'LineWidth',popTC_Lw+.5);
        
        ylim([min(min(monRad))-max(max(monRad))*.05 max(max(monRad))+max(max(monRad))*.05]);
        set(gca,'XLim',[0 length(monRadAve)+1],'XTick',1:length(monRadAve),'XTickLabel',radList,...
            'LineWidth',ax_Lw,'FontWeight','Bold','FontSize',popTC_Fsz,'YTick',ylim,'YTickLabel',round(ylim.*10)./10);
        xlabel('identity level (%)','FontSize',popTC_Fsz);
        title(['monkey ave.(n = ' num2str(length(monRad(:,1))) ' faces)'],'FontSize',popTC_Fsz);
        
        temp_dir                                = fullfile(paths.results,monkey,'aves');
        if ~exist(temp_dir,'dir')
            mkdir(temp_dir)
        end
        savename = fullfile(temp_dir,[exp_code '_' dat.h.snames '_radTrajAve' suffix.save]);
        saveas(31+figJit,[savename '.eps'],'epsc');
        export_fig(31+figJit,[savename '.png'],'-nocrop',['-r' num2str(prs)]);
        disp(['saved: ' savename]);
        close(31+figJit);
    end
    
    %% tangential trajectory plot
    humTan                                      = popTanHeat{i}(1:species.parseOffset,:);
    humTanAve                                   = mean(humTan);
    humTanErr                                   = std(humTan)./sqrt(length((humTan(:,1))));
    monTan                                      = popTanHeat{i}(species.parseOffset+1:end,:);
    monTanAve                                   = mean(monTan);
    monTanErr                                   = std(monTan)./sqrt(length(monTan(:,1)));
    %     popTanAve{i} = monTanAve;
    
    if figsOn==1
        
        figure(29+figJit);
        set(29+figJit,'Color','w','Position',figPos2,'Visible',visTog);
        
        subplot(2,1,1);
        plot(1:length(humTanAve),humTan,'Color',[.9 .9 .9],'LineWidth',popTC_Lw); hold on
        errorbar(humTanAve,humTanErr,'Color',[.1 .1 .1],'LineWidth',popTC_Lw+.5);
        ylim([min(min(humTan))-max(max(humTan))*.05 max(max(humTan))+max(max(humTan))*.05]);
        set(gca,'XTick',1:length(tanList),'XTickLabel',tanList,'LineWidth',ax_Lw,'FontWeight','Bold',...
            'FontSize',popTC_Fsz,'YTick',ylim,'YTickLabel',round(ylim.*10)./10);
        ylabel('ave. FR (Hz)','FontSize',popTC_Fsz);
        title(['human ave.(n = ' num2str(length(humTan(:,1))) ' faces)'],'FontSize',popTC_Fsz);
        
        subplot(2,1,2);
        plot(1:length(monTanAve),monTan,'Color',[.9 .9 .9],'LineWidth',popTC_Lw); hold on
        errorbar(monTanAve,monTanErr,'Color',[.1 .1 .1],'LineWidth',popTC_Lw+1.5);
        set(gca,'XTick',1:length(tanList),'XTickLabel',tanList,'LineWidth',ax_Lw,'FontWeight','Bold',...
            'FontSize',popTC_Fsz,'YTick',ylim,'YTickLabel',round(ylim.*10)./10);
        xlabel('tang. morph level (%)','FontSize',popTC_Fsz);
        title(['monkey ave.(n = ' num2str(length(monTan(:,1))) ' faces)'],'FontSize',popTC_Fsz);
        
        savename                                = fullfile(temp_dir,....
                        [exp_code '_' dat.h.snames '_tanTrajAve' suffix.save '.png']);
        export_fig(29+figJit,savename,'-nocrop',['-r' num2str(prs)]);
        disp(['saved: ' savename]);
        close(29+figJit);
    end
    
    close all
    
end

% figure
% floop1 = cell2mat(popRadSlopesMon(~cellfun(@isempty,popRadSlopesMon)));
% floop1 = floop1(:,~isnan(floop1(1,:)));
% floop2 = cell2mat(popRadSlopesNoiseMon(~cellfun(@isempty,popRadSlopesNoiseMon)));
% floop2 = floop2(:,~isnan(floop2(1,:)));
% 
% % lims = [floor(min(min([floop1 floop2]))*100) ceil(max(max([floop1 floop2]))*100)]./100;
% % xaxis = max(abs(lims))*-1:(max(abs(lims))*2)/9:max(abs(lims));
% subplot(2,2,1)
% hist(floop1(1,:))
% line([0 0],ylim,'LineStyle','--')
% xlim([-2 2])
% title('anti-face slope')
% % xlim([max(abs(lims))*-1 max(abs(lims))])
% subplot(2,2,2)
% hist(floop1(2,:))
% line([0 0],ylim,'LineStyle','--')
% xlim([-1 1])
% title('ident traj slope')
% % xlim([max(abs(lims))*-1 max(abs(lims))])
% subplot(2,2,3)
% hist(floop2(1,:))
% line([0 0],ylim,'LineStyle','--')
% xlim([-2 2])
% title('shuffled anti-face')
% % xlim([max(abs(lims))*-1 max(abs(lims))])
% subplot(2,2,4)
% hist(floop2(2,:))
% line([0 0],ylim,'LineStyle','--')
% xlim([-1 1])
% title('shuffled ident traj')
% % xlim([max(abs(lims))*-1 max(abs(lims))])


%% save results for this monkey
if ~isdir(fullfile(paths.results,monkey));
    mkdir(fullfile(paths.results,monkey));
end

savename                                        = fullfile(paths.results,'pop',...
                                    [monkey(1) paths.code 'PopData' suffix.save '.mat']);
save(savename);
disp(['Saved population data: ' exp_code 'PopData.mat']);

end