function SuperSimpleBlock
% Setting up figures
close all
pFig = figure;
mFig = figure;
b = .025;

figure(pFig);
plot1 = axes('Parent',pFig,'units','normalized','Position',[0+b 0+b .5-b .5-b]);
colormap jet; axis square;
set(plot1,'XTick',[]); set(gca,'YTick',[]);
slider1 = uicontrol(pFig,'Style','slider','Callback',@slider1_Callback,'units','normalized','Position',[0 0 .45 b]);

plot2 = axes('Parent',pFig,'units','normalized','Position',[0+b .5+b .5-b .5-b]);
colormap jet; axis square;
set(plot2,'XTick',[]); set(gca,'YTick',[]);
slider2 = uicontrol(pFig,'Style','slider','Callback',@slider2_Callback,'units','normalized','Position',[0 .5 .45 b]);

plot3 = axes('Parent',pFig,'units','normalized','Position',[.5+b 0+b .5-b .5-b]);
colormap jet; axis square;
set(plot3,'XTick',[]); set(gca,'YTick',[]);
slider3 = uicontrol(pFig,'Style','slider','Callback',@slider3_Callback,'units','normalized','Position',[.5 0 .45 b]);

plot4 = axes('Parent',pFig,'units','normalized','Position',[.5+b .5+b .5-b .5-b]);
colormap jet; axis square;
set(plot4,'XTick',[]); set(gca,'YTick',[]);
slider4 = uicontrol(pFig,'Style','slider','Callback',@slider4_Callback,'units','normalized','Position',[.5 .5 .45 b]);

figure(mFig)
list1 = uicontrol(mFig,'Style','listbox','Callback',@list1_Callback,'units','normalized','Position',[0 1 .5 .5]);

function list1_Callback(button1,~,~)




end