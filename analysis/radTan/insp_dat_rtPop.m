function insp_dat_rtPop(paths,monkeys,suffix)
%   loop through relevant files in 'paths.mas'
%
% Usage: insp_dat_rtPop(paths,monkey,suffix)
% MONKEY is the monkey's name
% PATHS is the project path structure containing at minimum paths.rare field
% SUFFIX is a string appended to the end of each filename
%
% last modified 6-18-16
% apj


% set experimental variables
visTog                                                  = 'off';
figsOn                                                  = 0; % switch for figures (on vs. off)
ax_Lw                                                   = 1; % axes linewidth
eB_Lw                                                   = 4.5; % error-bar linewidth
prs                                                     = 150;
fullscr                                                 = get(0,'ScreenSize'); % get screen size

dataname                                                = fullfile(paths.results,'pop',...
    [paths.code 'PopData' suffix.save '.mat']);
monkey_data                                             = load(dataname);

% load entire study
popData                                                 = cell(length(monkeys),1);
for mo = 1:length(monkeys)
    popData{mo}                                         = load(fullfile(paths.results,'pop',...
        [monkeys{mo}(1) paths.code 'PopData' suffix.load '.mat']));
end
% savename = fullfile(paths.results,'pop','trtPopData');
% popData{1} = load(savename);
% savename = fullfile(paths.results,'pop','rrtPopData');
% popData{2} = load(savename);
% savename = fullfile(paths.results,'pop','mrtPopData');
% popData{3} = load(savename);

figure('Visible',visTog,'Color','w');
fig_pos                                                 = get(gcf,'Position').*[1 1 .5 .5];
set(gcf,'Position',fig_pos);

%% impose inclusion criteria
% monkey (subject) loop
for mo = 1:length(popData)
    
    popData{mo}.inclVec                                 = zeros(1,length(popData{mo}.popRadHeat));
    for i = 1:length(popData{mo}.popRadHeat)
        if any(any(popData{mo}.popRadHeat{i}>=2))
            popData{mo}.inclVec(i)                      = 1;
        end
    end
    
    popData{mo}.did_flip                                = 0;
    popData{mo}.didnt_flip                              = 0;
    popData{mo}.flipped_cell                            = nan(length(popData{mo}.popRadHeat),1);
    
    %% plot population summary
%     thumb                                               = 1;
    popData{mo}.popPlotRadHeat                          = cell(1,sum(popData{mo}.inclVec~=0));
    popData{mo}.popPlotTanHeat                          = cell(1,sum(popData{mo}.inclVec~=0));
    for i = 1:length(popData{mo}.inclVec) % step through included neurons
        
        if popData{mo}.inclVec(i)~=0
            popData{mo}.normPopRadHeat{i}               = nan(size(popData{mo}.popRadHeat{i}));
            popData{mo}.normPopTanHeat{i}               = nan(size(popData{mo}.popTanHeat{i}));
            
            gloo = length(popData{mo}.popRadHeat{i});
            for ii = 1:gloo % step through trajectories
                
                %                 if ii<=length(popData{mo}.popRadHeat{i})
                % invert tuning curve
                tuning_curve                            = popData{mo}.popRadHeat{i}(ii,:);
                
                %                    [~,temp_ind] = max(tuning_curve);
                %                    if temp_ind==find(monkey_data.radList==0) % &&mo==1
                %                       flipped = tuning_curve.*-1;
                %                       popData{mo}.normPopRadHeat(ii,:) = flipped;
                %                       disp(['TC flipped- Cell#: ' num2str(i)])
                
                %                 foo                                     = mean(mean([popData{mo}.popRadHeat{i}(:,monkey_data.radList==0) ...
                %                                 popData{mo}.popRadHeat{i}(:,monkey_data.radList==0)]));
                %                 snoo                                    = mean(mean([popData{mo}.popRadHeat{i}(:,monkey_data.radList==0) ...
                %                                 popData{mo}.popRadHeat{i}(:,end-2:end)]));

                % if the norm is less than the caricatures, invert
                if tuning_curve(monkey_data.radList==0)>mean(tuning_curve(end-3:end))
                    %                          mean(popData{mo}.popRadHeat{i}(ii,end-3:end))
                    popData{mo}.normPopRadHeat{i}(ii,:) = tuning_curve.*-1;
                    disp(['TC flipped- Cell#: ' num2str(i)])
                    popData{mo}.flipped_cell(i)         = 1;
                    popData{mo}.did_flip                = popData{mo}.did_flip+1;
                else
                    popData{mo}.normPopRadHeat{i}(ii,:) = tuning_curve;
                    disp(['Did not flip- Cell#: ' num2str(i)])
                    popData{mo}.didnt_flip              = popData{mo}.didnt_flip+1;
                end
                                
                %% normalize
                %             if sum(normPopRadHeat(ii,:))~=0
                popData{mo}.normPopRadHeat{i}(ii,:)     = popData{mo}.normPopRadHeat{i}(ii,:)...
                    +-min(popData{mo}.normPopRadHeat{i}(ii,:));
                popData{mo}.normPopRadHeat{i}(ii,:)     = popData{mo}.normPopRadHeat{i}(ii,:)...
                    ./ max(popData{mo}.normPopRadHeat{i}(ii,:));
                %             end
           
                if ii<=length(popData{mo}.normPopTanHeat{i})
                    popData{mo}.normPopTanHeat{i}(ii,:) = popData{mo}.popTanHeat{i}(ii,:);

                    % normalize
                    popData{mo}.normPopTanHeat{i}(ii,:) = popData{mo}.normPopTanHeat{i}(ii,:)...
                        +-(min(popData{mo}.normPopTanHeat{i}(ii,:)));
                    popData{mo}.normPopTanHeat{i}(ii,:) = popData{mo}.normPopTanHeat{i}(ii,:)...
                        ./max(popData{mo}.normPopTanHeat{i}(ii,:));
                end
            end
            
            % extract TC
            foo = nanmean(popData{mo}.normPopRadHeat{i});
            
            % shuffle TC
            foo_noise = foo(randperm(length(foo)));
            
            % define trajectory subsets
            x1 = 1:4;
            x2 = 4:length(foo);
            
            shuffSlopes1 = nan(1000,1);
            shuffSlopes2 = nan(1000,1);
            for w = 1:10000
                foo_noise = foo(randperm(length(foo)));
                baa = polyfit(x1,foo_noise(x1),1);
                shuffSlopes1(w) = baa(1);
                baa = polyfit(x2,foo_noise(x2),1);
            	shuffSlopes2(w) = baa(1);
            end
            shufSlopeMean_anticar = mean(shuffSlopes1);
            shufSlopeErr_anticar = std(shuffSlopes1);
            shufSlopeMean_rad = mean(shuffSlopes2);
            shufSlopeErr_rad = std(shuffSlopes2);
            
            
            p1 = polyfit(x1,foo(x1),1);
            
            keyboard
            
            
            p2 = polyfit(x2,foo(x2),1);
            p1n = polyfit(x1,foo_noise(x1),1);
            p2n = polyfit(x2,foo_noise(x2),1);
            
            sigSlopes = nan(1,2);
            inSigSlopes = nan(1,2);

            if p1(1)<(shufSlopeMean_anticar - 1.96*shufSlopeErr_anticar)
                sigSlopes(1) = p1(1);
            else
                inSigSlopes(1) = p1(1);
            end
            if p2(1)>(shufSlopeMean_rad + 1.96*shufSlopeErr_rad)
                sigSlopes(2) = p2(1);
            else
                inSigSlopes(2) = p2(1);
            end

            
            popData{mo}.normPopRadSlope{i} = sigSlopes;
            popData{mo}.normPopRadSlope_insig{i} = inSigSlopes;

            popData{mo}.normPopRadSlopeNoise{i} = [p1n(1) p2n(1)];
            
%             
%             popData{mo}.popPlotTanHeat{thumb}           = popData{mo}.normPopTanHeat;
%             
%             %% select species
%             %             if mo==1
%             %                 popData{mo}.popPlotRadHeat{thumb}        = popData{mo}.normPopRadHeat(monkey_data.humFacePairs,:);
%             %             else
%             popData{mo}.popPlotRadHeat{thumb}       = popData{mo}.normPopRadHeat;
%             %             end
%             thumb                                       = thumb+1;
        end
    end
    
    % plotMat = cell2mat(popPlotRadHeat');
    % figure
    % subplot(211)
    % for i = 1:length(plotMat)
    %     plot(plotMat(i,:)); hold on
    % end
    % plot(mean(plotMat),'r','LineWidth',2)
    %
    % subplot(212)
    % hist(plotMat(:,end)); hold on
    % yline(mean(plotMat(:,end)),'r');
    
    %% radial traj population average
    
    % define population data
    %     popData{mo}.popRadCurve = mean(cell2mat(popData{mo}.popRadAve'));
    popData{mo}.popRadCurve                             = nanmean(cell2mat(popData{mo}.normPopRadHeat'));
    popData{mo}.popLength                               = sum(~isnan(cell2mat(popData{mo}.normPopRadHeat')));
    popData{mo}.popRadErr                               = nanstd(cell2mat(popData{mo}.normPopRadHeat'))./sqrt(popData{mo}.popLength(1));
    
    if figsOn
        figure('Position',fig_pos,'Visible',visTog); %('Visible',visTog)
        
        % plot
        errorbar(popData{mo}.popRadCurve,popData{mo}.popRadErr,'LineWidth',ax_Lw); hold on
        tempY                                           = ylim;
        yline(find(popData{mo}.radList==0));
        set(gca,'XTick',1:length(monkey_data.radList),'XTickLabel',monkey_data.radList,'LineWidth',ax_Lw,...
            'XLim',[.5 length(monkey_data.radList)+.5],'YLim',tempY);
        xlabel('identity level')
        ylabel('ave normalized FR')
        title(['M' num2str(mo) ' - radial trajectory ave. (n = ' num2str(sum(popData{mo}.inclVec~=0)) ' cells)'])
        
        % save
        savename                                        = fullfile(paths.results,'pop',...
            [popData{mo}.fileList(1).name(1:3) '_popNormRadAve' suffix.save '.png']);
        export_fig(gcf,savename,'-nocrop',['-r' num2str(prs)])
        disp(['Saved: ' savename])
        
        close(gcf)
    end
    
    %% tangential traj population average
    
    %% figure out why this is empty
    popData{mo}.popTanCurve                              = nanmean(cell2mat(popData{mo}.normPopTanHeat'));
    popData{mo}.popLength                                = sum(~isnan(cell2mat(popData{mo}.normPopTanHeat')));
    popData{mo}.popTanErr                                = nanstd(cell2mat(popData{mo}.normPopTanHeat'))./sqrt(popData{mo}.popLength(1));
    
    if figsOn
        figure('Position',fig_pos,'Visible',visTog); %('Visible',visTog)
        
        % plot
        errorbar(popData{mo}.popTanCurve,popData{mo}.popTanErr,'LineWidth',ax_Lw); hold on
        set(gca,'XTick',(1:length([0 monkey_data.tanList 100]))-1,'XTickLabel',[0 monkey_data.tanList 100],'LineWidth',ax_Lw,...
            'XLim',[.5 length([0 monkey_data.tanList 100])-1.5],'YLim',tempY);
        xlabel('identity level')
        ylabel('ave normalized FR')
        title(['M' num2str(mo) ' - tan trajectory ave. (n = ' num2str(mo) ' cells)'])
        
        % save
        savename                                        = fullfile(paths.results,'pop',...
            [popData{mo}.fileList(1).name(1:3) '_popNormTanAve' suffix.save '.png']);
        export_fig(gcf,savename,'-nocrop',['-r' num2str(prs)])
        disp(['Saved: ' savename])
        close(gcf)
    end
end

keyboard

%% plot population slopes

%% paired line plots

popSlopes           = cell2mat(popData{mo}.normPopRadSlope');
popSlopesInsig      = cell2mat(popData{mo}.normPopRadSlope_insig');
popSlopesNoise      = cell2mat(popData{mo}.normPopRadSlopeNoise');

figure('Position',[250   250   1100   400],'Color','w');
subplot(121)
line([0.5 2.5],[0 0],'Color',[.8 .8 .87]); hold on
lH = line(repmat([1 2],length(popSlopes),1)',...
[popSlopesNoise(:,1) popSlopes(:,1)]','Color','b','Marker','o');
set(gca,'XTick',[1 2],'XTickLabel',{'shuffled slopes' 'TC slopes'},...
    'YTick',[-.25 -.125 0 .125 .25])
text(1.45,.175,'***','Fontsize',18);
xlim([0.5 2.5])
ylim([-.25 .25])
title('anti-face slopes')

subplot(122)
line([0.5 2.5],[0 0],'Color',[.8 .8 .8]); hold on
lH = line(repmat([1 2],length(popSlopes),1)',...
[popSlopesNoise(:,2) popSlopes(:,2)]','Color','b','Marker','o');
set(gca,'XTick',[1 2],'XTickLabel',{'shuffled slopes' 'TC slopes'},...
    'YTick',[-.1 -.05 0 .05 .1])
text(1.45,-.05,'***','Fontsize',18);
xlim([0.5 2.5])
ylim([-.1 .1])
title('ident traj slopes')


[p,h] = ranksum(popSlopesNoise(:,1),popSlopes(:,1)) 

[p,h] = ranksum(popSlopesNoise(:,2),popSlopes(:,2)) 

savename                    = fullfile(paths.results,'pop',['totPopRadSlopeLine' suffix.save]);
saveas(gcf,[savename '.eps'],'epsc');
export_fig(gcf,[savename '.png'],'-nocrop',['-r' num2str(prs)])
close(gcf)

%% histograms of slopes

yLims = [];

hF = figure('Position',[250   250   1100   400],'Color','w');
subplot(1,2,1)
% hAx1 = axes('Position',[0.1 0.1 0.3 .8]);
hist(popSlopesInsig(:,1),20)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.7 .7 .7],'EdgeColor','w','facealpha',0.75)
hold on;
hist(popSlopes(:,1),20)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75);
axis tight
xlim([-max(abs(xlim)) max(abs(xlim))])
line([0 0],ylim,'Color','k')
yLims = max([yLims ylim]);
ylabel('# of neurons')
title('anti-face slope')

subplot(1,2,2)
% hAx2 = axes('Position',[0.45 0.1 0.3 .8]);
hist(popSlopesInsig(:,2),20)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.7 .7 .7],'EdgeColor','w','facealpha',0.75)
hold on;
hist(popSlopes(:,2),20)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75);
axis tight
xlim([-.1 .1])
yLims = max([yLims ylim]);
line([0 0],[0 yLims],'Color','k')
title('ident traj slope')

% set(hAx1,'ylim',[0 yLims])

% get the axes handle
a = get(hF, 'CurrentAxes');
% get the handles for children of the axes -- these are the data series handles
c = get(a, 'Children');
% generate a legend command using these "children"
h = legend(c([2 3]),'significant','insignificant','Location','NorthWest');

% legend(hAx2, 'TC slopes', {'shuffled'; 'TC slopes'}); 
% set(h,'position',[0.85 0.4 .01 .3])

savename                    = fullfile(paths.results,'pop',['totPopRadSlope10kshuff' suffix.save]);
saveas(gcf,[savename '.eps'],'epsc');
export_fig(gcf,[savename '.png'],'-nocrop',['-r' num2str(prs)])
close(gcf)
%%

% num2str(popData{1}.did_flip/(popData{1}.did_flip+popData{1}.didnt_flip))
% num2str(popData{2}.did_flip/(popData{2}.did_flip+popData{2}.didnt_flip))
% num2str(popData{3}.did_flip/(popData{3}.did_flip+popData{3}.didnt_flip))

if figsOn
    
    % plot radial aves
    figure('Color','w','Visible',visTog); hold on % 'Visible',visTog,
    errorbar(popData{1}.popRadCurve,popData{1}.popRadErr,'Color',[.8 1 .8],'LineWidth',eB_Lw);
    errorbar(popData{2}.popRadCurve,popData{2}.popRadErr,'Color',[1 .8 .8],'LineWidth',eB_Lw);
    errorbar(popData{3}.popRadCurve,popData{3}.popRadErr,'Color',[.8 .8 1],'LineWidth',eB_Lw);
    errorbar(mean([popData{1}.popRadCurve; popData{2}.popRadCurve; popData{3}.popRadCurve]),...
        mean([popData{1}.popRadErr; popData{2}.popRadErr; popData{3}.popRadErr]),'Color',[.5 .5 .5],'LineWidth',eB_Lw)
    tempY                                                   = ylim;
    yline(find(monkey_data.radList==0));
    set(gca,'XTick',1:length(monkey_data.radList),'XTickLabel',monkey_data.radList,'LineWidth',ax_Lw,'YLim',tempY) ;
    xlabel('identity level')
    ylabel('ave normalized FR')
    
    tot_inclusion                                           = [popData{1}.inclVec popData{2}.inclVec popData{3}.inclVec];
    title(['radial trajectory ave. (n = ' num2str(sum(tot_inclusion)) ' cells)'])
    legend(['M1 (Toroid [' num2str(sum(popData{1}.inclVec)) ' cells - AF])'],...
        ['M2 (Rhombus [' num2str(sum(popData{2}.inclVec)) ' cells - AF])'],...
        ['M3 (Matcha [' num2str(sum(popData{3}.inclVec)) ' cells - AM])'],...
        'Ave.','Location','Southeast')
    
    savename                    = fullfile(paths.results,'pop',['totPopNormRadAve' suffix.save]);
    saveas(gcf,[savename '.eps'],'epsc');
    export_fig(gcf,[savename '.png'],'-nocrop',['-r' num2str(prs)])
    close(gcf)
    
end
%
% % plot tang. aves
% figure('Color','w'); hold on % 'Visible',visTog,
% errorbar(tanPopCurve,tanPopErr,'Color',[.8 1 .8],'LineWidth',eB_Lw);
% errorbar(rhombPop.tanPopCurve,rhombPop.tanPopErr,'Color',[1 .8 .8],'LineWidth',eB_Lw);
% errorbar(mean([tanPopCurve; rhombPop.tanPopCurve]),...
%     mean([tanPopErr; rhombPop.tanPopErr]),'Color',[.25 .7 1],'LineWidth',eB_Lw*2)
% errorbar(mean([tanPopCurve; rhombPop.tanPopCurve]),...
%     mean([tanPopErr; rhombPop.tanPopErr]),'Color',[.5 .5 .5],'LineWidth',eB_Lw/2)
% set(gca,'XTick',1:length([0 monkey_data.tanList 100]),'XTickLabel',[0 monkey_data.tanList 100],'LineWidth',ax_Lw,'YLim',tempY) ;
% xlabel('identity level')
% ylabel('ave normalized FR')
% title(['M1 & M2 - tan trajectory ave. (n = ' num2str(sum(inclVec~=0)) '+' ...
%     num2str(sum(rhombPop.inclVec~=0)) ' cells)'])
% legend('M1','M2','Ave.','Location','Southeast')
%
% savename = [paths.results 'pop/totPopNormTanAve.eps'];
% saveas(gcf,savename,'epsc');
% savename = [paths.results 'pop/totPopNormTanAve.png'];
% export_fig(gcf,savename,'-nocrop','-nocrop',['-r' num2str(prs)])
% close(gcf)
%
% humRadMat = [];
% monRadMat = [];
% humTanMat = [];
% monTanMat = [];
%
% slopeList.slopeAnti = [];
% slopeList.slopeCar = [];
% slopeList.slopeAntiShuff = [];
% slopeList.slopeCarShuff = [];
% slopeList.slopeTan = [];
%
% for i = 1:length(fileList)
%
%     % radial trajectories
%     radCellHeat = popRadHeat{i};
%     humRadMat = [humRadMat; radCellHeat(1:12,:)];
%     monRadMat = [monRadMat; radCellHeat(13:end,:)];
%
%     faceRadAve = mean(radCellHeat);
%     humRadAve = mean(radCellHeat(1:12,:));
%     monRadAve = mean(radCellHeat(13:end,:));
%
%     %     flippedRadAve = abs(faceRadAve - faceRadAve(monkey_data.radList==0));
%     %     flippedRadAve = flippedRadAve./max(flippedRadAve);
%     normFaceRadAve = faceRadAve+-(min(faceRadAve));
%     normFaceRadAve = normFaceRadAve./max(normFaceRadAve);
%
%     if normFaceRadAve(monkey_data.radList==0)>mean(normFaceRadAve(7:end))
%         normFaceRadAve = normFaceRadAve.*-1;
%     end
%
%     y1 = normFaceRadAve(monkey_data.radList<=0);
%     x1 = 1:length(y1);
%     coeff = polyfit(x1,y1,1);
%     slopeAnti = coeff(1);
%
%     y1 = normFaceRadAve(monkey_data.radList>=0);
%     x1 = 1:length(y1);
%     coeff = polyfit(x1,y1,1);
%     slopeCar = coeff(1);
%
%     figure(1414)
%     plot(slopeAnti,slopeCar,'b.'); hold on
%     xlabel('anti');
%     ylabel('cari');
%     %     axis equal
%     xline(0);
%     yline(0);
%
%     %     y1 = monRadAve(monkey_data.radList<=0);
%     %     x1 = 1:length(y1);
%     %     coeff = polyfit(x1,y1,1);
%     %     slopeAnti = coeff(1);
%     %
%     %     y1 = monRadAve(monkey_data.radList<=100);
%     %     x1 = 1:length(y1);
%     %     coeff = polyfit(x1,y1,1);
%     %     slopeCar = coeff(1);
%     %
%     %     plot(slopeAnti,slopeCar,'g.')
%
%     y1 = normFaceRadAve(monkey_data.radList<=0);
%     x1 = randperm(length(y1));
%     coeff = polyfit(x1,y1,1);
%     slopeAntiShuff = coeff(1);
%
%     y1 = normFaceRadAve(monkey_data.radList>=100);
%     x1 = randperm(length(y1));
%     coeff = polyfit(x1,y1,1);
%     slopeCarShuff = coeff(1);
%
%     plot(slopeAntiShuff,slopeCarShuff,'.','Color',[.7 .7 .7])
%
%     % tangential trajectories
%     tanCellHeat = popTanHeat{i};
%     humTanMat = [humTanMat; tanCellHeat(1:12,:)];
%     monTanMat = [monTanMat; tanCellHeat(13:end,:)];
%
%     faceTanAve = mean(tanCellHeat);
%     humTanAve = mean(tanCellHeat(1:12,:));
%     monTanAve = mean(tanCellHeat(13:end,:));
%
%     flippedTanAve = faceTanAve+-(min(faceTanAve));
%     flippedTanAve = flippedTanAve./max(flippedTanAve);
%
%     y1 = flippedTanAve;
%     x1 = 1:length(y1);
%     coeff = polyfit(x1,y1,1);
%     slopeTan = coeff(1);
%
%     slopeList.slopeAnti = [slopeList.slopeAnti; slopeAnti];
%     slopeList.slopeCar = [slopeList.slopeCar; slopeCar];
%     slopeList.slopeAntiShuff = [slopeList.slopeAntiShuff; slopeAntiShuff];
%     slopeList.slopeCarShuff = [slopeList.slopeCarShuff; slopeCarShuff];
%     slopeList.slopeTan = [slopeList.slopeTan; slopeTan];
%
% end
%
% savename = [paths.results 'pop/' exp_code 'SlopeScatter.png'];
% export_fig(1414,savename,'-nocrop',['-r' num2str(prs)])
%
% xAxis = -.25:.01:.25;
% antiHist = hist(slopeList.slopeAnti,xAxis);
% carHist = hist(slopeList.slopeCar,xAxis);
% shuffHist = hist(slopeList.slopeCarShuff,xAxis);
% tanHist = hist(slopeList.slopeTan,xAxis);
%
% figure
% stairs(xAxis,shuffHist,'Color',[.7 .7 .7]); hold on
% stairs(xAxis,tanHist,'Color',[.8 .8 0]);
% stairs(xAxis,antiHist,'Color',[.75 0 .75]);
% stairs(xAxis,carHist,'Color',[0 .8 0]);
% legend('shuffle slope','tang. slope','antiface slope','car. slope')
%
% savename = [paths.results 'pop/' exp_code 'SlopeBar.png'];
% export_fig(gcf,savename,'-nocrop',['-r' num2str(prs)]);
% close(gcf)
%
% prs = 150
%
% C = {'b','r','g','y','m'};
% plotList = [1 7 12 13 14];
% figure
% tempPos =  [0 425 1600 200];
% set(gcf,'Position',tempPos,'PaperUnits','inches','PaperPosition',tempPos/prs,'Visible',visTog)
% for a = 1:length(plotList)
%     % subplot(length(plotList),1,a)
%     tempPlot = popRadAve{plotList(a)}-min(popRadAve{plotList(a)});
%     tempPlot = tempPlot./max(tempPlot);
%     plot(tempPlot,'Color',C{a},'LineWidth',lw); hold on
%     ylim([0 1])
%     set(gca,'XTick',1:radLen+1,'XTickLabel',monkey_data.radList,'YTick',ylim,'LineWidth',ax_Lw) ;
%     xlabel('identity level')
%     ylabel('norm FR (Hz)')
%     title('human ave.(n = 12 faces)')
%     suplabel([monkey 'radial morphs'],'t')
% end
% xTemp = xlim;
% xlim([xTemp(1)+.5 xTemp(2)-.5])
% savename = [paths.results 'aves/' exp_code 'bestRadTrajAve.png'];
% export_fig(gcf,savename,'-nocrop',['-r' num2str(prs)]);
% disp(['saved: ' savename])
% close(gcf)
%
% figure;      % open figure for population heatmap
% set(gcf,'Position',fullscr,'PaperUnits','inches','PaperPosition',fullscr/prs,'Visible',visTog)
% for i = 1:length(popRadHeat)
%     tempCMapMax = max([max(popRadHeat{i}) max(popTanHeat{i})]);
%     [cmap] = buildcmap('kbwr');
%     colormap(cmap);   % set colormap
%     subplot(2,length(popRadHeat),i);
%     pH = pcolor(popRadHeat{i});
%     set(pH,'EdgeColor','none')
%     title(fileList(i).name(4:end-4))
%     % axis square
%     hold on
%     caxis([0 max(tempCMapMax)]) % min(tempCMapMin)
%     %             shading interp
%     if i ==1
%         set(gca,'XTick',1:radLen,'XTickLabel',monkey_data.radList,...
%             'YTick',1:length(faceList),'YTickLabel',faceList) ;
%         % rotateticklabel(gca,90)
%         xlabel('identity level (radial)','FontWeight','bold','FontSize',axFsz);
%         ylabel('face','FontWeight','bold','FontSize',axFsz);
%     else
%         set(gca,'XTick',[],'YTick',[])
%     end
% end
% for i = 1:length(popTanHeat)
%     tempCMapMax = max([max(popRadHeat{i}) max(popTanHeat{i})]);
%     [cmap] = buildcmap('kbwr');
%     colormap(cmap);   % set colormap
%     subplot(2,length(popRadHeat),i+length(popRadHeat));
%     pH = pcolor(popTanHeat{i});
%     set(pH,'EdgeColor','none')
%     hold on
%     caxis([0 max(tempCMapMax)]) % min(tempCMapMin)
%
%     if i ==1
%         set(gca,'XTick',1:tanLen,'XTickLabel',monkey_data.tanList,...
%             'YTick',1:length(faceList),'YTickLabel',faceList) ;
%         xlabel('identity level (tang)','FontWeight','bold','FontSize',axFsz);
%         ylabel('face','FontWeight','bold','FontSize',axFsz);
%     else
%         set(gca,'XTick',[],'YTick',[])
%     end
% end
% % do some fancy saving here
% savename = [paths.results  'popHeatRadTan.png'];
% export_fig(gcf,savename,'-nocrop',['-r' num2str(prs)]);
% disp(['saved: ' savename])
% close all
