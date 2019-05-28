function [x,y] = TCpicker(data,npoints,rectsize,c)
imagesc(squeeze(nanmean(data(:,:,1:size(data,3)/10),3))); axis square; caxis([-c c]);
c = colormap(jet(npoints));
for i = 1:npoints
    [x(i,1),y(i,1)] = ginput(1);
    rectangle('Position',[x(i),y(i),rectsize,rectsize],'FaceColor',c(i,:))
    drawnow
end
x(:,2) = x(:,1)+rectsize;
y(:,2) = y(:,1)+rectsize;
x = round(x); y = round(y);