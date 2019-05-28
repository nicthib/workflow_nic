function SCAPE_NNLS_seedpicker(im,window)
close all
ss = size(im);
hFig = figure('Position',[0 1920 1920 400]);
polyFig = figure('Position',[0 0 1920 400]);
win1 = 5;
win2 = 2;
s = 1+win1;
imtmp = [];
seedpix = [];
imforreggrow = [];
Jout = [];
a = 1;
binmask = ones(ss(1:3));
slider1 = uicontrol(hFig,'Style','slider','Callback',@slider1_Callback,'units','normalized','Position',[.15 0 .85 .05]);
button1 = uicontrol(hFig,'Style','pushbutton','Callback',@button1_Callback,'units','normalized','Position',[0 .1 .1 .05],'String','Pick TC');
button2 = uicontrol(hFig,'Style','pushbutton','Callback',@button2_Callback,'units','normalized','Position',[.1 .1 .1 .05],'String','Save');
text1 = uicontrol(hFig,'Style','edit','Callback',@text1_Callback,'units','normalized','Position',[.2 .1 .1 .05],'String','');

set(slider1,'value',1);
set(slider1,'max',size(im,4)); %
set(slider1,'min',1);
getim
    function getim
        if window == 1
            im2 = squeeze(mean(im(:,:,:,s:s+win2),4)).*binmask;
            im1 = squeeze(mean(im(:,:,:,s-win1:s-1),4)).*binmask;
            imside = squeeze(max(im2-im1,[],3)');
            imleft = squeeze(max(im2-im1,[],1));
            imtmp = [imside imleft];
            imforreggrow = im2-im1;
        else
            imside = squeeze(max(im(:,:,:,s),[],3)');
            imleft = squeeze(max(im(:,:,:,s),[],1));
            imtmp = [imside imleft];
            imforreggrow = squeeze(max(im(:,:,:,s),4));
        end
    end
    function slider1_Callback(slider1,eventdata,handles)
        s = round(get(slider1,'Value'));
        update_fig
    end
    function text1_Callback(text1,eventdata,handles)
        s = str2double(get(text1,'String'));
        update_fig
    end

    function button1_Callback(button1,~,~)
        [X,Y] = ginput(2);
        seedpix = round([X(1) Y(1) X(2)-ss(1)]);
        [P,J] = regionGrowing(imforreggrow, seedpix,[],20);
        if size(P,1) < 1000
            figure(polyFig)
            subplot(121)
            shp = alphaShape(P(:,1),P(:,2),P(:,3),1);
            plot(shp,'EdgeColor','none','FaceColor',rand(1,3))
            hold on
            axis image
            Jout{a} = J;
            a = a+1;
            binmask = binmask-J;
        else
            disp('TOO LARGE')
        end
    end

    function button2_Callback(button2,~,~)
        save('Output.mat','Jout');
    end

    function update_fig
        figure(hFig)
        getim
        imagesc(imtmp); axis image; colormap jet
        set(gca,'XTick',[]); set(gca,'YTick',[]);
        caxis([-50 50])
    end
end