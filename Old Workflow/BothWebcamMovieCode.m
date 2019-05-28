j = 1;
cd([(m.STIMpath{m.ttu(j)}) '/' m.stimlist{m.ttu(j)}(find(m.stimlist{m.ttu(j)}~='_')) '_webcam']);
for i=0:length(dir)-4
    webcam_movie1(:,:,i+1)=rgb2gray(imread([num2str(i) '.jpg']));
end

stimnum = regexp(m.stimlist{m.ttu(j)},'\d','Match');
cd([(m.CCDpath{m.ttu(j)}) '/' m.stimlist{m.ttu(j)}(1:4) '_webcam' stimnum{1}]);
for i=0:length(dir)-4
    webcam_movie2(:,:,i+1)=rgb2gray(imread([m.stimlist{m.ttu(j)}(1:4) num2str(i) '.jpg']));
end

%% For running analysis


%%

close all
gcamp_f = smooth3(data.gcamp,'box',[3 3 3]);
chbt_f = smooth3(data.chbt,'box',[3 3 3]);
m.jetc = colormap(jet); m.jetc(1,:) = [0 0 0];
m.bonec = colormap(bone); m.bonec(1,:) = [0 0 0];
%%
figure('Color','black','Position',[100,100,800,800])

for i = 1:size(webcam_movie1,3)
    ax1 = subplot(221);
    imagesc(reshape(reshape(squeeze(comp(:,:,ceil(i*m.fr/30),:)*30),[512*512 m.ncomps])*bla,[512 512 3]))
    axis square; axis off;  set(gca,'xcolor','k','ycolor','k','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    text(10,20,sprintf('%0.1f sec',round(i*10/30)/10),'Color','white','FontSize',12)
%     imagesc(gcamp_f(:,:,ceil(i/2.8831)))
%     axis image; axis off; set(gca,'xcolor','k','ycolor','k','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
%     text(10,20,sprintf('%0.1f sec',round(i*10/30)/10),'Color','white','FontSize',12)
%     caxis([-.05 .05])

    ax2 = subplot(222);
    imagesc(chbt_f(:,:,ceil(i/2.8831)))
    axis square; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    text(10,20,sprintf('%0.1f sec',round(i*10/30)/10),'Color','white','FontSize',12)
    caxis([-10 10]*1e-6)

    ax3 = subplot(223);
    imagesc(webcam_movie1(:,:,i));
    text(10,20,sprintf('%0.1f sec',round(i*10/30)/10),'Color','black','FontSize',12)
    axis image
    set(gca,'xcolor','k','ycolor','k','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    
    ax4 = subplot(224);
    imagesc(webcam_movie2(:,:,i));
    text(10,20,sprintf('%0.1f sec',round(i*10/30)/10),'Color','black','FontSize',12)
    axis image
    set(gca,'xcolor','k','ycolor','k','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    
    colormap(ax1,m.bonec)
    colormap(ax2,m.jetc)
    colormap(ax3,'gray')
    colormap(ax4,'gray')
    
    
    M(i) = getframe(gcf);
    clf
end

%squeeze(nanmean(C(i*3,:,:),2))'