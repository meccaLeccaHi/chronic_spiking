function plot_spheres(data,plot_limits,savename,prs,paths,faceToggle,varargin)
% playing around with sphere function to make nice facespace plot
%
% Usage: plot_spheres
% DATA is 2 cells of 12x9 matrices of response values
% PLOT_LIMITS are the color axis limits [min max]
% SAVENAME and PRS are self-explanatory
% PATHS is the project path structure containing at minimum paths.rare field
% FACETOGGLE is a 1 or a 0, indicating whether to plot faces or not
%
% last modified 6-22-16
% apj

% opengl software
% warning('off','MATLAB:prnRenderer:opengl')

if ~isempty(varargin{1})
    if regexp(varargin{1},'o')
        visTog              = varargin{1};
    else
        disp('varargin unrecognized')
        return
    end
else
    visTog                  = 'off';
end

Fsz                         = 10;           % tick label font size
axFsz                       = 12;           % axis label font size

if isunix
    stim_dir                = '/mnt/raid1/projects/stimuli/stimCode/FaceLearning/Face_On_Transparent';
else
    stim_dir                = '\\srvr1\LeopoldShare\projects\stimuli\stimCode\FaceLearning\Face_On_Transparent';
end


% faceLearningStim_dir_list = dir(fullfile(stim_dir,'*.png'));
addpath(fullfile(paths.physio,'ptb'))
gL                          = getFaceLearningList;

vecs =    [1.15     0       -.7;
    .25     -1      -.5;
    -1      0       -.5;
    0       1.3     -.8;
    .25     0       1;
    -.25    -.05    1.1;
    1       0       .5;
    0       1       .5;
    -1      0       .75;
    .25    -1       .5;
    .2      -.2     -1;
    -.2     .2      -1];

faceStim                    = 8:8:192;
[x,y,z]                     = sphere;
cmin                        = plot_limits(1);
cmax                        = plot_limits(2);
caxis([cmin cmax]);
% cc                          = jet(10);
ncol                        = 32;
cc                          = flipud(hot(ncol));
cutoffs                     = nan(ncol,2);
level                       = (cmax-cmin)/ncol;
for u = 1:ncol
    cutoffs(u,1)            = cmax-((cmax-cmin)*u/ncol);
    cutoffs(u,2)            = cmax-((cmax-cmin)*u/ncol)+level;
end
cutoffs(1,2)                = 1.01*cmax;
cutoffs(end,1)              = cmin-abs(cmin/100);
pos                         = [375 350 1200 550];

figure(37);
if regexp(visTog,'off')
    pos(1)                  = pos(1)*20;
end
clf(37);



%% ?-is z-buffer necessary?
set(37,'Position',pos,'Color',[1 1 1],'Visible',visTog); % ,'renderer','zbuffer'
hIm                         = cell(1,2);

% loop through each data species present
dataNum                     = length(data(:,1))/12;
for s = 1:dataNum
    %     if nargin>1 && s>1
    %         data                = data2;
    %     end
    loopData                = data([1:12]+(s-1)*12,:);
    morph_len               = length(loopData(1,:));
    subplot(1,dataNum,s)
    lStep                   = 8; % distance of each step away from center
    
    % plot 'halo'
    r                       = 10;
    h                       = mesh(x*r,y*r,z*8);
    set(h,'FaceAlpha',0,'EdgeAlpha',.35,'EdgeColor',[.7 .7 .7],...
        'DiffuseStrength',0,'AmbientStrength',0,'EdgeLighting','none','FaceLighting','none');
    hold on;
    
    temp                    = mean(loopData(:,1));
    normCol                 = temp>=cutoffs(:,1)&temp<=cutoffs(:,2);
    if sum(normCol)>=0
        normCol             = find(temp>=cutoffs(:,1),1,'first')-1;
    end
    % sphere centered at origin
    b                       = surf(x,y,z,'Facecolor',cc(normCol,:),'LineStyle','none');
    set(b,'DiffuseStrength',.9,'AmbientStrength',.6 ,'EdgeLighting','phong','FaceLighting','phong');
    
    % specify colormap based on range of data
    for i = 1:morph_len-1
        for ii = 1:length(loopData(:,1));
            col             = loopData(ii,i+1)>=cutoffs(:,1)&loopData(ii,i+1)<cutoffs(:,2);
            if sum(col)==0
                col         = find(loopData(ii,i+1)>=cutoffs(:,1),1,'first')-1;
            elseif sum(col)>1
                col         = find(loopData(ii,i+1)>=cutoffs(:,1),1,'first');
            end
            b               = surf(x+lStep*i*vecs(ii,1)/morph_len,y+lStep*i*vecs(ii,2)/morph_len,z+lStep*i*vecs(ii,3)/morph_len,...
                    'Facecolor',cc(col,:),'LineStyle','none');
            set(b,'DiffuseStrength',.9,'AmbientStrength',.6 ,'FaceLighting','phong');
            hold on
        end
    end
    
    light
    
    if faceToggle
        hIm{s}                  = nan(1,length(loopData(:,1)));
        % plot face images here
        i = morph_len+1;
        for ii = 1:length(loopData(:,1));
            if isempty(strfind(savename,'Monk'))
                imNum           = faceStim(ii+(s-1)*length(loopData(:,1)));
            else
                imNum           = faceStim(ii+12+(s-1)*length(loopData(:,1)));
            end
            
            [a,~,alpha]         = imread(fullfile(stim_dir,gL{imNum}));
            
            xLoc                = lStep*i*vecs(ii,1)/morph_len;
            yLoc                = lStep*i*vecs(ii,2)/morph_len;
            zLoc                = lStep*i*vecs(ii,3)/morph_len;
            
            foo                 = flip(cat(3,im2double(a),im2double(alpha)),1);  % flipdim
            
            hIm{s}(ii)          = imsurf(foo,[xLoc+.5 yLoc-1.6 zLoc-1.65],[-1 -1 0],[-1 1 0],.008);
            
            set(hIm{s}(ii),'LineStyle','none','CDataMapping','direct','EdgeColor','none',...
                'DiffuseStrength',0,'AmbientStrength',0,'EdgeLighting','none','FaceLighting','none');
        end
    end
    
    axis tight equal off;
    %     axLims = axis;
    
    if s==1
        %         hct = title('humans','FontWeight','bold'); hold on
        
        %         colPos = [.5 .165 .11 .015];
        %                 ha = axes('Position',colPos);
        %                 set(ha,'XTick',[],'YTick',[])
        %                 colormap(hot);
        %                 hcb = colorbar('Location','North','Position',colPos);
        %                 tempLabel = str2num(num2str([(minResp(1)*100)/100 (maxResp*100)/100],3));
        %                 set(get(hcb,'XLabel'),'String','{\itHz}','FontWeight','bold');
        %                 set(hcb,'XTick',get(hcb,'xlim'),'XTickLabel',tempLabel,...
        %                     'YTick',[],'FontWeight','bold');
        
        % create legend
        colormap(hot);
        hcb                 = colorbar('Location','South');
        pos                 = get(hcb,'Position');
        if dataNum==1
            pos             = pos.*[3.7 1.4 .125 .5];
        else
            pos             = pos.*[1.825 1.3 .25 .5];
        end
        
        set(hcb,'Position',pos);
        set(get(hcb,'YLabel'),'String','{\itHz}','FontWeight','bold','FontSize',axFsz,'Rot',0)
        tempPos1            = get(get(hcb,'YLabel'),'Position');
        set(get(hcb,'YLabel'),'Position',tempPos1.*[1 -1 1]);
        yAxCb               = get(hcb,'XLim');
        temp                = round([cmin cmax].*10)./10;
        set(hcb,'XTick',yAxCb,'XTickLabel',temp,'FontWeight','bold','FontSize',Fsz);
    else
        tempPos             = get(gca,'Position');
        set(gca,'Position',tempPos.*[.7 1 1 1])
        
        %         %% create movie
        %         %axLims = axis;
        %         winSize = get(37,'Position');
        %         winSize(1:2) = [0 0];
        %         numFrames = 180;
        %         A = moviein(numFrames,37,winSize);
        %         set(37,'NextPlot','replacechildren')
        %
        %         for i = 1:numFrames
        % %             if i==numFrames/4;
        % %                 % delete face images here
        % %                 for ii = 1:length(data(:,1));
        % %                     delete(hIm{1}(ii))
        % %                     delete(hIm{2}(ii))
        % %                 end
        % %                 hIm{1} = nan(1,length(data(:,1)));
        % %                 hIm{2} =hIm{1};
        % %
        % %                 % plot face images here
        % %                 i = morph_len+1;
        % %                 for ii = 1:length(data(:,1));
        % %                     for p = 1:2
        % %                         subplot(1,2,p)
        % %                         imNum = faceStim(ii+(p-1)*length(data(:,1)));
        % %                         temp = im2double(flipdim(img.img{imNum},1));
        % %
        % %                         xLoc = lStep*i*vecs(ii,1)/morph_len;
        % %                         yLoc = lStep*i*vecs(ii,2)/morph_len;
        % %                         zLoc = lStep*i*vecs(ii,3)/morph_len;
        % %
        % %                         hIm{p}(ii) = imsurf(temp,[xLoc yLoc-1 zLoc-1],[0 1 0],[1  0 0],.01);
        % %                         % imsurf(imageIn,upperLeftPoint3,normal,imXDirVec,scale,varargin)
        % %                         set(hIm{p}(ii),'LineStyle','none','CDataMapping','direct','EdgeColor','none',...
        % %                             'DiffuseStrength',0,'AmbientStrength',0,'EdgeLighting','none','FaceLighting','none');
        % %                     end
        % %                 end
        % %                 axis tight equal off
        % %                 %             elseif i==numFrames/4*3
        % %                 %                 % delete face images here
        % %                 %                 for ii = 1:length(data(:,1));
        % %                 %                     delete(hIm{1}(ii))
        % %                 %                     delete(hIm{2}(ii))
        % %                 %                 end
        % %                 %                 hIm{1} = nan(1,length(data(:,1)));
        % %                 %                 hIm{2} =hIm{1};
        % %                 %
        % %                 %                 % plot face images here
        % %                 %                 i = morph_len+1;
        % %                 %                 for ii = 1:length(data(:,1));
        % %                 %                     for p = 1:2
        % %                 %                         subplot(1,2,p)
        % %                 %                         imNum = faceStim(ii+(p-1)*length(data(:,1)));
        % %                 %                         temp = im2double(flipdim(img.img{imNum},1));
        % %                 %
        % %                 %                         xLoc = lStep*i*vecs(ii,1)/morph_len;
        % %                 %                         yLoc = lStep*i*vecs(ii,2)/morph_len;
        % %                 %                         zLoc = lStep*i*vecs(ii,3)/morph_len;
        % %                 %
        % %                 %                         hIm{p}(ii) = imsurf(temp,[xLoc yLoc-1 zLoc-1],[-1 -1 0],[-1 1 0],.01);
        % %                 %                         % imsurf(imageIn,upperLeftPoint3,normal,imXDirVec,scale,varargin)
        % %                 %                         set(hIm{p}(ii),'LineStyle','none','CDataMapping','direct','EdgeColor','none',...
        % %                 %                             'DiffuseStrength',0,'AmbientStrength',0,'EdgeLighting','none','FaceLighting','none');
        % %                 %                   view  end
        % %                 %                     axis tight equal off
        % %             end
        %
        %             subplot(1,2,1)
        %             view(i,15)
        %             %set(gca,'XLim',axLims(1:2),'YLim',axLims(3:4),'ZLim',axLims(5:6));
        %
        %             subplot(1,2,2)
        %             view(i,15)
        %             %set(gca,'XLim',axLims(1:2),'YLim',axLims(3:4),'ZLim',axLims(5:6));
        %             drawnow
        %
        %             % add to movie matrix
        %             A(:,i) = getframe(37,winSize);
        %             % add to gif
        %             im = frame2im(A(:,i));
        %             [imind,cm] = rgb2ind(im,256);
        %             if i==1;
        %                 imwrite(imind,cm,[movDir 'FSgif' spikeName 'Movie.gif'],'gif',...
        %                     'Loopcount',inf,'DelayTime',.05,'TransparentColor',0);
        %             else
        %                 imwrite(imind,cm,[movDir 'FSgif' spikeName 'Movie.gif'],'gif',...
        %                     'WriteMode','Append','DelayTime',.05,'TransparentColor',0);
        %             end
        %             % end
        %         end
        % %         movie(37,A,1,9,winSize)
        % %         save (fullfile(movDir,[spikeName 'mMovie.mat']) ,'A')
        % %         movie2avi(A,fullfile(movDir,[spikeName 'Movie.avi']))
        
    end
end


% do some fancy saving
tempPos                     = get(37,'Position');
if regexp(visTog,'off')
    tempPos(1)              = tempPos(1)*100;
end
set(37,'Position',tempPos/prs,'Position',tempPos,'Visible',visTog);
export_fig(37,savename,['-r' num2str(prs)])

cuts                        = strfind(savename,filesep);
disp(['saved = ' savename(cuts(end)+1:end)])
close(37)

end