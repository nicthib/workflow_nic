%% Taking H, making a scrolling note animation.
close all
scrollrate = 1;
full_fig = zeros([100 500+3.*1800]);
fig_sz = [100 500];
line_thk = 1; staff = [20:20+line_thk 35:35+line_thk 50:50+line_thk 65:65+line_thk 80:80+line_thk];
note_y = 12.5:7.5:97.5;
full_fig(staff,:) = 1;

full_fig = insertShape(full_fig,'FilledCircle',[(round(MIDImat(1,:)*3*30/(31.22*100))+250)',note_y(mod(MIDImat(4,:),12)+1)',5*ones(size(MIDImat,2),1)],'Opacity',1);
imshow(full_fig); axis image;
%%
xlim([0 500])
fig = figure('Color','Black','Position',[100 100 600 800]);

for i = 1:5399
    subplot(4,1,1:3)
    if mod(i,3) == 1
        imagesc(squeeze(compvid(:,:,:,ceil(i*m.fr/90)))*brightness)
        axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
        text(95,256,sprintf('%0.1f seconds',i/30),'Color','w','FontSize',12)
    end
    title('LSQ Output','Color','w')
    
    subplot(4,1,4)
    imshow(full_fig(:,i:i+500))
    line([250 250],[0 100],'Color','w')
    M(i) = getframe(gcf);
end