function qthresh(im)
close all
im = im;
hFig = figure; thresh1 = 1; thresh2 = 255;

s = 2; imtmp = im(:,:,s);
hIm = imagesc(imtmp);
axis image; colormap gray;
set(gca,'XTick',[]); set(gca,'YTick',[]);
slider1 = uicontrol(hFig,'Style','slider','Callback',@slider1_Callback,'units','normalized','Position',[.15 0 .85 .05]);
button1 = uicontrol(hFig,'Style','pushbutton','Callback',@button1_Callback,'units','normalized','Position',[0 .1 .15 .05],'String','TC');
edit1 = uicontrol(hFig,'Style','edit','Callback',@edit1_Callback,'units','normalized','Position',[0 0.05 .15 .05]);
edit2 = uicontrol(hFig,'Style','edit','Callback',@edit2_Callback,'units','normalized','Position',[0 0 .15 .05]);

set(slider1,'value',1);
set(slider1,'max',size(im,3)); %
set(slider1,'min',1);
    function slider1_Callback(slider1,eventdata,handles)
        s = round(get(slider1,'Value'));
        update_fig
    end

    function button1_Callback(button1,~,~)
        im(im<thresh1) = 0;
        im(im>thresh2) = 0;
        im(thresh2 >= im >= thresh1) = 1;
%         im(im < thresh) = 0;
%         im(im >= thresh) = 1;
        im = logical(im);
        assignin('base','thresh1',thresh1)
        assignin('base','thresh2',thresh2)
        close(hFig);
    end

    function edit1_Callback(edit1,~,~)
        thresh1 = round(str2num(get(edit1,'String')));
        update_fig
    end
    function edit2_Callback(edit1,~,~)
        thresh2 = round(str2num(get(edit1,'String')));
        update_fig
    end


    function update_fig
        figure(hFig)
        imtmp = im(:,:,s);
        imtmpmask = imtmp;
        imtmpmask(imtmpmask<thresh1) = 0;
        imtmpmask(imtmpmask>thresh2) = 0;
        imtmpmask(thresh2 >= imtmpmask >= thresh1) = 1;
%         imtmpmask(imtmp<thresh(1)) = 0;
%         imtmpmask(imtmp>thresh(1)) = 255;
        imagesc(imfuse(imtmp,imtmpmask)); axis image; %colormap gray;
        set(gca,'XTick',[]); set(gca,'YTick',[]);
    end
end