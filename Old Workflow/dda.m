% dda(IM,m,CAX) shows the image sequence IM using the caxis CAX. IM must be a
% 3 dimensional array, and can be either raw LED data, converted
% hemodynamics or gcamp data. CAX is the data you are choosing to view,
% wither 'hbo', 'hbr', 'hbt', or 'gcamp'.
%
% The bottom slider allows you to move through the data set, and the text
% box allows for custom input of the data displayed.
%
% When 1 number is input, the frame of that run is shown. You can also use
% the slider to adjust what frame you are looking at.
%
% When 2 numbers are input, the average of that frame range is shown. The
% slider is not usable in this mode.
%
% When 3 numbers are input, the mean of the first two numbers range is
% subtracted from the third number (think of it as a pseudomovie). Moving
% the slider in this mode changes the frame of interest, while keeping the
% baseline frames unchanged.
%
% When 4 numbers are used, the mean of the 3rd and 4th number range is
% subtracted from the mean of the 1st and 2nd number range. Using the
% slider changes the frames used for the 3rd and 4th numbers, while
% keeping the range equal.

function m = dda(im,m,cax)
close all
sflag = 0;
a = [1 .5 .5];
if isstruct(im)
    sflag = 1;
    imchbo = im.chbo;
    imchbr = im.chbr;
    imchbt = im.chbt;
    imgcamp = im.gcamp;
    try
        imgR = im.R;
        imgX = im.C(:,:,1);
        imgY = im.C(:,:,2);
    catch
    end
    if ndims(imchbo) == 4
        imchbo = squeeze(nanmean(imchbo,4));
        imchbr = squeeze(nanmean(imchbr,4));
        imchbt = squeeze(nanmean(imchbt,4));
        imgcamp = squeeze(nanmean(imgcamp,4));
    end
end

if ndims(im) == 4
    
    im = squeeze(nanmean(im,4));
end

s = 1; caxscale = 1; range = 0; npoints = 1;
if isstr(cax)
    if strcmp(cax,'chbo')
        cax = [-1 1]*1e-5;
        if isstruct(im)
            im = imchbo;
        end
    elseif strcmp(cax,'chbr')
        cax = [-1 1]*1e-5;
        if isstruct(im)
            im = imchbr;
        end
    elseif strcmp(cax,'chbt')
        cax = [-1 1]*1e-5;
        if isstruct(im)
            im = imchbt;
        end
    elseif strcmp(cax,'gcamp')
        cax = [-.02 .02];
        if isstruct(im)
            im = imgcamp;
        end
    end
end
hFig = figure;
imagesc(im(:,:,s));
colorbar; colormap jet; caxis(cax); axis image;
set(gca,'XTick',[]); set(gca,'YTick',[]);
slider1 = uicontrol(hFig,'Style','slider','Callback',@slider1_Callback,'units','normalized','Position',[.15 0 .85 .05]);
slider2 = uicontrol(hFig,'Style','slider','Callback',@slider2_Callback,'units','normalized','Position',[.95 .1 .025 .82]);
edit1 = uicontrol(hFig,'Style','edit','Callback',@edit1_Callback,'units','normalized','Position',[0 0 .15 .05]);
edit2 = uicontrol(hFig,'Style','edit','Callback',@edit2_Callback,'units','normalized','Position',[0 .1 .15 .05]);
set(edit2,'String',1)
uicontrol(hFig,'Style','pushbutton','Callback',@button1_Callback,'units','normalized','Position',[0 .05 .15 .05],'String','Single TC');
uicontrol(hFig,'Style','pushbutton','Callback',@button2_Callback,'units','normalized','Position',[0 .15 .15 .05],'String','Multi TC');
uicontrol(hFig,'Style','pushbutton','Callback',@button3_Callback,'units','normalized','Position',[0 .95 .15 .05],'String','DONE');
set(slider1,'value',1);
set(slider1,'max',size(im,3)); %
set(slider1,'min',1);
set(slider2,'max',log10(10)); %
set(slider2,'min',log10(.1));
set(slider2,'value',log10(1))
    function slider1_Callback(slider1,~,~)
        if numel(s) == 1
            s = round(get(slider1,'Value'));
            set(edit1,'String',mat2str(s));
        elseif numel(s) == 2
            s(1) = round(get(slider1,'Value'));
            s(2) = round(get(slider1,'Value'))+range;
            set(edit1,'String',mat2str(s));
        elseif numel(s) == 3
            s(3) = round(get(slider1,'Value'));
            set(edit1,'String',mat2str(s));
        elseif numel(s) == 4
            s(3) = round(get(slider1,'Value'));
            s(4) = round(get(slider1,'Value'))+range;
            set(edit1,'String',mat2str(s));
        end
        update_fig
    end

    function slider2_Callback(slider2,~,~)
        caxscale = 10^get(slider2,'Value');
        update_fig
    end

    function edit1_Callback(edit1,~,~)
        try
            s = round(str2num(get(edit1,'String')));
        catch
        end
        if numel(s) == 4
            range = abs(s(4)-s(3));
        elseif numel(s) == 2
            range = abs(s(2)-s(1));
        end
        update_fig
    end

    function edit2_Callback(edit2,~,~)
        npoints = str2num(get(edit2,'String'));
    end

    function button1_Callback(~,~)
        m.x = []; m.y = [];
        [m.x,m.y] = ginput(2); m.x = m.x'; m.y = m.y';
        m.x = m.x + .5; m.y = m.y + .5;
        m.x = floor(m.x/(size(im,1)/32))*size(im,1)/32;
        m.y = floor(m.y/(size(im,1)/32))*size(im,1)/32;
        p = patch([m.x(1),m.x(2)+1,m.x(2)+1,m.x(1)]-.5,[m.y(1),m.y(1),m.y(2)+1,m.y(2)+1]-.5,'k');
        set(p,'FaceAlpha',.2);
        figure
        if sflag == 0
            TC = squeeze(squeeze(nanmean(nanmean(im(m.y(1):m.y(2),m.x(1):m.x(2),:),2),1)));
            try
                plot(m.t,TC)
            catch
                plot(TC)
            end
        else
            numsub = 4;
            try
                subplot(numsub,1,3)
                tcam = 1/30:1/30:30; % Needs work
                imgR = hampel(imgR,10);
                imgR(imgR == 0) = NaN;
                imgR(end,:) = NaN;
                for i = 1:size(imgR,2)
                    patch(tcam,imgR(:,i),'Black','EdgeAlpha',.1);
                    hold on
                end
                plot(tcam,nanmean(imgR,2),'LineWidth',3)
                ylim([10 25])
                
                subplot(numsub,1,3)
                plot(mean(imgX,2)-mean(imgX(1,:)))
                hold on
                plot(mean(imgY,2)-mean(imgY(1,:)))
                legend('X','Y')
            catch
                numsub = 2;
            end
            m.TCo = squeeze(squeeze(nanmean(nanmean(imchbo(m.y(1):m.y(2),m.x(1):m.x(2),:),2),1)));
            m.TCr = squeeze(squeeze(nanmean(nanmean(imchbr(m.y(1):m.y(2),m.x(1):m.x(2),:),2),1)));
            m.TCt = squeeze(squeeze(nanmean(nanmean(imchbt(m.y(1):m.y(2),m.x(1):m.x(2),:),2),1)));
            m.TCg = squeeze(squeeze(nanmean(nanmean(imgcamp(m.y(1):m.y(2),m.x(1):m.x(2),:),2),1)));
            subplot(numsub,1,1)
            plot(m.t,[m.TCo,m.TCr,m.TCt])
            subplot(numsub,1,2)
            plot(m.t,m.TCg)      
        end
    end

    function button2_Callback(~,~)
        m.x = []; m.y = [];
        for i = 1:npoints
            [x,y] = ginput(1);
            x = x + .5; y = y + .5;
            m.x(i,1) = floor(x*32/size(im,1))*16;
            m.y(i,1) = floor(y*32/size(im,1))*16;
            rectangle('Position',[m.x(i)-.5,m.y(i)-.5,size(im,1)/32,size(im,1)/32],'FaceColor',abs([a(1)-i/npoints a(2)-i/npoints a(3)-i/npoints]))
            text(m.x(i)-size(im,1)/32,m.y(i)-size(im,1)/32,mat2str(i),'Color','k');
            drawnow
        end
        m.x(:,2) = m.x(:,1)+size(im,1)/32-1;
        m.y(:,2) = m.y(:,1)+size(im,1)/32-1;
        
        figure
        if sflag == 0;
            clear legstr;
            for i = 1:npoints
                legendstring{i} = mat2str(i);
                plot(smooth(squeeze(squeeze(nanmean(nanmean(im(m.y(i,1):m.y(i,2),m.x(i,1):m.x(i,2),:),2),1)))),'Color',abs([a(1)-i/npoints a(2)-i/npoints a(3)-i/npoints]));
                hold on
                
            end
            legend(legendstring)
        else
            numsub = 6;
            try
                subplot(numsub,1,5)
                tcam = 1/30:1/30:30; % Needs work
                imgR = hampel(imgR,10);
                imgR(imgR == 0) = NaN;
                imgR(end,:) = NaN;
                for i = 1:size(imgR,2)
                    patch(tcam,imgR(:,i),'Black','EdgeAlpha',.1);
                    hold on
                end
                plot(tcam,nanmean(imgR,2),'LineWidth',3)
                ylim([10 25])
                
                subplot(numsub,1,6)
                plot(mean(imgX,2)-mean(imgX(1,:)))
                hold on
                plot(mean(imgY,2)-mean(imgY(1,:)))
                legend('X','Y')
            catch
                numsub = 4;
            end
            
            for i = 1:npoints
                subplot(numsub,1,1)
                plot(m.t,squeeze(squeeze(nanmean(nanmean(imchbo(m.y(i,1):m.y(i,2),m.x(i,1):m.x(i,2),:),2),1))),'Color',abs([a(1)-i/npoints a(2)-i/npoints a(3)-i/npoints]))
                hold on
                subplot(numsub,1,2)
                plot(m.t,squeeze(squeeze(nanmean(nanmean(imchbr(m.y(i,1):m.y(i,2),m.x(i,1):m.x(i,2),:),2),1))),'Color',abs([a(1)-i/npoints a(2)-i/npoints a(3)-i/npoints]))
                hold on
                subplot(numsub,1,3)
                plot(m.t,squeeze(squeeze(nanmean(nanmean(imchbt(m.y(i,1):m.y(i,2),m.x(i,1):m.x(i,2),:),2),1))),'Color',abs([a(1)-i/npoints a(2)-i/npoints a(3)-i/npoints]))
                hold on
                subplot(numsub,1,4)
                plot(m.t,squeeze(squeeze(nanmean(nanmean(imgcamp(m.y(i,1):m.y(i,2),m.x(i,1):m.x(i,2),:),2),1))),'Color',abs([a(1)-i/npoints a(2)-i/npoints a(3)-i/npoints]))
                hold on
                legstr{i} = mat2str(i);
            end
            subplot(numsub,1,1); title('CHBO'); legend(legstr); 
            subplot(numsub,1,2); title('CHBR');
            subplot(numsub,1,3); title('CHBT'); 
            subplot(numsub,1,4); title('GCaMP');             
        end       
    end

    function button3_Callback(~,~)
        assignin('base','m',m)
        close all
    end

    function update_fig
        figure(hFig)
        if numel(s) == 1
            imagesc(im(:,:,s));
            set(slider1,'value',s);
        elseif numel(s) == 2
            try
                imagesc(mean(im(:,:,s(1):s(2)),3));
            catch 
            end
        elseif numel(s) == 3
            imagesc(im(:,:,s(3))-mean(im(:,:,s(1):s(2)),3))
            set(slider1,'value',s(3));
        elseif numel(s) == 4
            try
                imagesc(mean(im(:,:,s(3):s(4)),3)-mean(im(:,:,s(1):s(2)),3));
            catch
            end
            set(slider1,'value',s(3));
        end
        colorbar; colormap jet; caxis(cax*caxscale); axis image
        set(gca,'XTick',[]); set(gca,'YTick',[]);
    end
end