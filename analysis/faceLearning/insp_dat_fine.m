function insp_dat_fine(paths,monkey,suffix)
%
% Usage: insp_dat_fine(paths,monkey,suffix)
%
% PATHS project path structure containing at minimum paths.rare field
% MONKEY is the monkey's name
% SUFFIX is a string appended to the end of each filename
%
% last modified 6-18-16
% apj

% % housekeeping
% close all
% warning('off','MATLAB:opengl:switchToSoftwareUnixNotSupported')
% opengl software

% % set experimental variables
% % msz                    = 6;        % MarkerSize
% % dirStr = [paths.dataAux 'results/resultsFaceLearning/fine/']; % directory for saving plots


visTog                          = 'off';                     % visibility of plots
plotFaces                       = 1;
plotObjects                     = 1;
prs                             = 100;                      % print resolution
figJit                          = round(rand*10);           % avoid using the same fig number in parallel
ax_Fsz                          = 8;                        % axes fontsize
fileName_Fsz                    = 13;                       % file label fontsize
title_Fsz                       = 8;                        % raster title fontsize
date_Fsz                        = 12;                       % date label fontsize
ax_Lw                           = 1;                        % axes linewidth
sdf_Lw                          = 2;                        % lineWidth for sdf's
fullscr                         = get(0,'ScreenSize');      % get screen size
figPos                          = [fullscr(3)/3   fullscr(4)/2   fullscr(3)/2   fullscr(4)/4]; % [x y width height]
if regexp(visTog,'off')
    figPos(1)                   = figPos(1)+5000;
end
RasterWin                       = [-300 600]; % raster plot display size (msec)

% face numbers
humans                          = 1:12;
monkeys                         = 101:112;
faces                           = [humans monkeys];
objects                         = 201:212;

% get list of data files
exp_code                        = [monkey(1) paths.code];
dirList                         = dir(fullfile(paths.mas,[exp_code '*' suffix.load '.mat']));

% set save path
temp_dir                        = fullfile(paths.results,'fine',monkey);
fprintf(['Saving figures to:\n' temp_dir '\n']);

% load stimuli
load(fullfile(paths.stim,'imgFaceLearningSet1OnWhite.mat')); % 'imgFaceLearningSet1.mat'
stim_images                     = img;

% open figure
figure(747+figJit);

% loop through spike files
for i = 1:length(dirList);
    
    % set file-name
    spike_name                  = dirList(i).name(1:end-4);
    fprintf(['Loading file: \n' spike_name '.mat\n']);
    filename                    = fullfile(paths.mas,[spike_name '.mat']);
    
    %% load spike file
    load(filename)
    unpack

    % set variables
    morphs                      = [0 dat.h.stim.morphSteps];   % [1,2,3,6,12,25,50,100]
    if str2double(dat.h.date(1,1:2))==20
        dates                   = unique(datenum(dat.h.date,'yyyy-mmm-dd'));
    else
        dates                   = unique(datenum(dat.h.date,'dd-mmm-yyyy'));
    %     dates                       = dates(min([length(dates(:,1)) 5]),:);
    end
    dates                       = dates(1:min([length(dates) 7]));
    
    stim_win                    = dat.c(1,T_STIMON:T_STIMOFF);
    colors                      = [0 0 1; jet(length(dates(:,1))-1)];
    
    if plotFaces
        %% plot faces
        % loop through stimuli
        for s = 1:length(faces)
            
            if ismember(faces(s),humans)
                norm_id             = 0;
            elseif ismember(faces(s),monkeys)
                norm_id             = 100;
            else
                norm_id             = nan;
            end
            
            %if strcmp(stim_id(1:3),'ref') % count morphs of current stimuli
            %    length(morphs) = 1;
            %else
            %    length(morphs) = length(strmatch(stim_id(1:3),stimlist));
            %end
            
            % setup figure
            clf(747+figJit);
            set(747+figJit,'Position',figPos,'color','w','Visible',visTog);
            rows                    = length(dates(:,1)) + 2;
            cols                    = length(morphs) + 1;
            ax                      = subplot(rows,cols,1);
            
            % plot name of spike-file
            text(-1.25,.5,spike_name,'HorizontalAlignment','left',...
                'Fontsize',fileName_Fsz,'Interpreter','none');
            set(ax,'Visible','off')
            
            cumuLsdf                = cell(1,length(morphs)); % preallocate sdf array over dates
            sdf_y                   = nan(2,length(morphs)); % preallocate for consistent ylims
            
            %% plot rasters for each day & morph-level
            
            % loop through morph-levels
            for m = 1:length(morphs)
                
                if m==1
                    mdat            = select_trials(dat,FACE,norm_id,DATE,dates);
                else
                    mdat            = select_trials(dat,FACE,faces(s),STEP,morphs(m),DATE,dates(:,1));
                end
                imNum               = unique(mdat.c(:,STIM));
                %  imNum = imNum(end);
                
                % plot face image
                h1                  = subplot(rows,cols,m+1);
                sub_pos             = get(h1,'Position');
                image(stim_images{imNum(1)});
                axis off; axis image;
                if m==length(morphs)
                    title(['face# ' num2str(faces(s))],'FontSize',title_Fsz)
                end
                set(h1,'Position',sub_pos.*[1 1 1 1]);  % [x y width height]
                plot_left           = sub_pos(1);  % get position for later
                     
                % loop through days
                for d = 1:length(dates(:,1));
                    
                    % select data for each day
                    subdat          = select_trials(mdat,DATE,dates(d,:));
                    
                    % plot raster
                    subplot(rows,cols,m+d*cols+1);
                    plot_rasters(subdat.s,RasterWin); %hold all;
                    templims = ylim;
                    yline(stim_win);
                    sub_pos         = get(gca,'Position');
                    sub_pos(1)      = plot_left;   % align with subplots above
                    if m==1&&d==length(dates(:,1));
                        xtick       = round(stim_win);
                    else
                        xtick       = [];
                    end
                    
                    xlims = xlim;
                    set(gca,'XTick',xtick,'YTick',[],'Fontsize',ax_Fsz,...
                        'Position',sub_pos,'Linewidth',ax_Lw);
                    if m==1
                        set(gca, 'YTick',templims(1),'YTickLabel',length(subdat.s(:,1)));
                    end
                    if d==1
                        % plot title
                        title([num2str(unique(subdat.c(:,STEP))) '%'],'FontSize',title_Fsz)
                    end
                    box on
                    
                    % get sdf
                    %                 l = size(subdat.s);
                    if numel(subdat.s)~=0; % if subdat.s contains responses,then extract sdf
                        cumuLsdf{m}(d,:) = get_sdf(subdat,1,RasterWin);
                    end
                    
                    % date label
                    if m == 1;
                        ax          = subplot(rows,cols,m+d*cols);
                        foo         = datestr(dates(d,:));
                        text(-1,.5,foo(1:end-5),'Color',colors(d,:),...
                            'FontSize',date_Fsz)
                        set(ax,'Visible','off')
                    end
                    
                end
                
                % get range of data for setting plot limits
                sdf_y(1,m)          = min(min(cumuLsdf{m}));
                sdf_y(2,m)          = max(max(cumuLsdf{m}));
                
            end
            
            %% plot cumulative sdf's
            
            % get overall range of data for setting plot limits
            plot_lims               = [ceil(min(sdf_y(1,:))*10)/10 floor(max(sdf_y(2,:))*10)/10];
            
            % loop through morphs
            for m = 1:length(morphs)
                
                if ~isempty(cumuLsdf{m})
                    
                    % plot sdf
                    subplot(rows,cols,m+(1+length(dates(:,1)))*cols+1);
                    hold on
                    for w = 1:length(dates(:,1))
                        if max(cumuLsdf{m}(w,:))~=0
                            plot(xlims(1):xlims(2),cumuLsdf{m}(w,:),'Color',colors(w,:),'LineWidth',sdf_Lw);
                        end
                    end
                    
                    % set limits
                    xlim([subdat.c(1,T_PRE) subdat.c(1,T_POST)]);
                    if sum(plot_lims)~=0
                        y_ax        = plot_lims;
                        set(gca,'YLim',y_ax)
                    else
                        y_ax        = [];
                    end
                    
                    if m==1
                        ylabel('FR (hZ)','Fontsize',ax_Fsz,'FontWeight','bold');
                    else
                        y_ax        = [];
                    end
                    
                    % align with subplots above
                    cur_pos         = get(gca,'Position');
                    new_pos         = [cur_pos([1 2]) sub_pos(3) cur_pos(4)];
                    set(gca,'Position',new_pos,'XTick',[],'YTick',y_ax,...
                        'Linewidth',ax_Lw,'Fontsize',ax_Fsz,'XLim',xlims);
                    hold on
                    xline(0);
                    yline(subdat.c(1,T_STIMON:T_STIMOFF));
                    box off
                end
            end
            
            % do some fancy saving here
            filename = [exp_code '_' dat.h.snames '_stim' num2str(s) '_fine' suffix.save '.png'];
            savename = fullfile(temp_dir,filename);
            export_fig(747+figJit,savename,'-nocrop',['-r' num2str(prs)])
            disp(['Saved: '  filename]);
        end
    end
    
    if plotObjects
        %% plot objects stim
        
        cumuLsdf                    = cell(1,length(objects)); % preallocate array of sdfs over all dates
        sdf_y                       = nan(2,length(objects)); % preallocate for ylims
        
        clf(747+figJit);
        set(747+figJit,'Position',figPos,'color','w','Visible',visTog);
        rows                        = length(dates(:,1)) + 2;
        cols                        = length(objects)+1;
        
        % display spike name
        ax                          = subplot(rows,cols,1);
        text(-1.25,.5,spike_name,'HorizontalAlignment','left','Fontsize',fileName_Fsz,'Interpreter','none');
        set(ax,'Visible','off')
        
        for o = 1:length(objects)  % step through morphs
            
            mdat                    = select_trials(dat,FACE,objects(o));
            imNum                   = unique(mdat.c(:,STIM));
            
            % plot object image
            h1                      = subplot(rows,cols,o+1);
            sub_pos                 = get(h1,'Position');
            image(stim_images{max(imNum)});
            axis off; axis image;
            if o==length(objects)
                title(['objects #' num2str(objects(1)) ':' num2str(objects(end))],'FontSize',title_Fsz)
            end
            set(h1,'Position',sub_pos.*[1 1 1 1]); % [x y width height]
            plot_left = sub_pos(1);
            
            % loop through days
            for d = 1:length(dates(:,1));
                
                subdat              = select_trials(mdat,DATE,dates(d,:));
                
                % plot raster
                subplot(rows,cols,o+d*cols+1);
                plot_rasters(subdat.s,RasterWin); %hold all;
                yline(stim_win);
                sub_pos             = get(gca,'Position');
                sub_pos(1)          = plot_left;
                if o==1;
                    xtick           = round(stim_win);
                else
                    xtick           = [];
                end
                xlims               = xlim;
                set(gca,'XTick',xtick,'YTick',[],'Fontsize',ax_Fsz,...
                    'Position',sub_pos.*[1 1 1 1],'Linewidth',ax_Lw);
                
                % get sdf
                if numel(subdat.s)~=0; % if subdat.s contains responses,then extract sdf
                    cumuLsdf{o}(d,:) = get_sdf(subdat,1,RasterWin);
                end
                
                % date label
                if o == 1;
                    ax              = subplot(rows,cols,o+d*cols);
                    foo             = datestr(dates(d,:));
                    text(-1,.5,foo(1:end-5),'Color',colors(d,:),'FontSize',date_Fsz)
                    set(ax,'Visible','off')
                end                  
                
            end
            
            % get range of data for setting plot limits
            sdf_y(1,o)              = min(min(cumuLsdf{o}));
            sdf_y(2,o)              = max(max(cumuLsdf{o}));
        end
        
        plot_lims = [ceil(min(sdf_y(1,:))*10)/10 floor(max(sdf_y(2,:))*10)/10];
        
        for o = 1:length(objects)  % step through morphs
            
            if ~isempty(cumuLsdf{o})
                
                % plot sdf
                subplot(rows,cols,o+(1+length(dates(:,1)))*cols+1);
                hold on
                for w = 1:length(dates(:,1))
                    if max(cumuLsdf{o}(w,:))~=0
                        plot(xlims(1):xlims(2),cumuLsdf{o}(w,:),'Color',colors(w,:),...
                            'LineWidth',sdf_Lw);
                    end
                end
                xlim([subdat.c(1,T_PRE) subdat.c(1,T_POST)]);
                cur_pos = get(gca,'Position');
                new_pos = [cur_pos([1 2]) sub_pos(3) cur_pos(4)];
                
                if o==1&&sum(plot_lims)~=0
                    y_ax = plot_lims;
                    ylabel('FR (hZ)','Fontsize',ax_Fsz,'FontWeight','bold');
                    set(gca,'YLim',y_ax);
                else
                    y_ax = [];
                end
                set(gca,'Position',new_pos,'XTick',[],'YTick',y_ax,...
                    'Linewidth',ax_Lw,'Fontsize',ax_Fsz,'XLim',xlims);
                hold on
                xline(0);
                yline(subdat.c(1,T_STIMON:T_STIMOFF));
                box off
            end
        end
        
        % do some fancy saving here
        savename                = fullfile(temp_dir, ...
                            [exp_code '_' dat.h.snames '_stim' num2str(objects(o)) '_fine' suffix.save '.png']);
        export_fig(747+figJit,savename,'-nocrop',['-r' num2str(prs)]);
        disp(['Saved: '  savename]);
    end
end

end