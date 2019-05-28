close all
keep =  [2,3,4,8,10,11,15,16,18,21,22,23];
cmap = hsv(length(keep));
figure
mm=0
for i = keep
mm=mm+1;    
plot([1:length(H)]/10,H(i,:)+.3*mm,'color',cmap(mm,:));
    hold on; 
end
xlabel('time(s)');

figure
W = W/max(max(W));
imgfull = zeros([m.sz1 m.sz2 3]);
imgfulltop = zeros([m.sz1 m.sz2 3]);
%imgfullleft = zeros([m.sz1 m.sz2 3]);
imgfullcat = [];
mm=0;
for i = keep
    mm=mm+1;
    imgtmp = reshape(Wl(:,i)*cmap(mm,:),[m.sz1 m.sz2 3]);
    imgfullcat = [imgfullcat;rot90(imgtmp,3)];
end
imagesc([flipud(imgfullcat*4)])
hold on
plot([20:53],[250:283]*0+1280, 'w', 'linewidth', 2)
daspect([3.1,2.9,1])
axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);