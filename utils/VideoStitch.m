function VideoStitch(H1,W1,H2,W2,webcam1,webcam2,m,filename)
close all
if ~isfield(m,'moviespeed')
    m.moviespeed = 1;
end
vidObj = VideoWriter([filename '.avi']); 
vidObj.FrameRate = m.framerate-2*m.framerate*(1-m.interp)/3;
open(vidObj)
fig = figure('Color','Black','Position',[0 0 1280 720]);
offset = 130; scale = .4;
movielength = size(H1,2)*m.nLEDs/m.framerate;
t = m.nLEDs/m.framerate:m.nLEDs/m.framerate:movielength;

if m.interp
    for i = 1:size(H1,1)
        H1_v(i,:) = interp(H1(i,:),3); 
        H2_v(i,:) = interp(H2(i,:),3);
    end
    H1_i(isnan(H1_i)) = 0; H2_i(isnan(H2_i)) = 0;

else
    H1_v = H1; H2_v = H2;
end

for i = 1:size(H1_v,2)
    cmap1 = m.cmap; cmap2 = m.cmap;
    cmap1(H1_v(:,i)<0,:) = 1; cmap2(H2_v(:,i)<0,:) = 1;
    H1_v(:,i) = abs(H1_v(:,i));
    H2_v(:,i) = abs(H2_v(:,i));
    frame1 = im2u8sc(reshape(W1*diag(H1_v(:,i))*cmap1,[m.sz m.sz 3]),m.caxis1).*uint8(imresize(m.BW,1/m.dsf));
    frame2 = im2u8sc(reshape(W2*diag(H2_v(:,i))*cmap2,[m.sz m.sz 3]),m.caxis2).*uint8(imresize(m.BW,1/m.dsf));
    if ~isempty(webcam1)
        frame3 = imresize(padarray(imcrop(imrotate(webcam1(:,:,i*(3-2*m.interp)),12),[offset offset offset+(640*scale) offset+(480*scale)]),[20 20],'both'),[480,640]);
        frame4 = imresize(padarray(imcrop(webcam2(:,:,i*(3-2*m.interp)),[65 20 540 380]),[20 20]),[480 640]);
        frame3 = repmat(frame3,[1 1 3]);
        frame4 = repmat(frame4,[1 1 3]);
    end
    if ~isempty(webcam1)
        fullframe = [frame3,frame4; padarray(frame1(11:490,:,:),[0 64],'both') padarray(frame2(11:490,:,:),[0 64],'both')];
    else
        fullframe = [frame1 frame2];
    end
    imagesc(fullframe)
    caxis([0 1])
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    if isfield(m,'movietext')
        text(.5,1,sprintf([m.movietext{i} ', t = %0.1f s'],t(i)*m.moviespeed),'Color','w','FontSize',16,'HorizontalAlignment','Center','Units','normalized');
    else
        text(.5,1,sprintf('t = %0.1f s',t(i)*m.moviespeed),'Color','w','FontSize',16,'HorizontalAlignment','Center','Units','normalized');
    end
    writeVideo(vidObj,getframe(fig));
end
close all
close(vidObj)
