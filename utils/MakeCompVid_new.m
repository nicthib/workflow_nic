function MakeCompVid_new(H,W,m,gcamp,resid,filename)
close all
fig1 = figure('Color','Black');
fig2 = figure('Color','Black');
fig = figure('Color','Black','Position',[0,0,1400,700]);

cmapgcamp = colormap('gray');
for i = 1:size(H,2)
    
    figure(fig1)
    imagesc(gcamp(:,:,i))
    axis image; caxis([-.01 .05]); axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); colormap(cmapgcamp);
    set(fig1,'Position',[100,100,700,700]);
    tmp = getframe(fig1); cdata1 = tmp.cdata;
    
    figure(fig2)
    imagesc(reshape(W*diag(H(:,i))*m.cmap,[m.sz m.sz 3])*m.brightness)
    axis image; axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    set(fig2,'Position',[100,100,700,700]);
    tmp = getframe(fig2); cdata2 = tmp.cdata;
    
    figure(fig)
    imagesc(cat(2,cdata1,cdata2));
    axis image; axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    text(700,-20,sprintf('t = %0.1f s',i*3/m.framerate),'Color','w','FontSize',14,'HorizontalAlignment','Center');
    M(i) = getframe(gcf);
end
close(fig);
vid = VideoWriter(filename);
vid.FrameRate = m.framerate/3; % if interp is used, framerate is 3x higher.
open(vid); writeVideo(vid,M); close(vid)