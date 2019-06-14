tmp = m.BW_nnmf(:,:,1);
%imagesc(tmp)
B = bwboundaries(tmp);
imagesc(max(data.gcamp,[],3).*m.BW_nnmf(:,:,1))
caxis([-.3 .3])
hold on
x = []; y = [];
for i = 1:10
    
    
    [ytmp,xtmp] = ginput(1);
    y = [y ytmp]; x = [x xtmp];
    scatter(y(end),x(end),'k')
    
    
end
[V,C,XY]=VoronoiLimit(x,y,'bs_ext',B{1},'figure','on');
axis image